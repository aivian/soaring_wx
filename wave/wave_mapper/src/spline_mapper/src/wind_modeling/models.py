# basic ros/numpy/etc imports
import rospy
import numpy as np
import scipy.io as sio
import scipy.sparse as sps
import copy
# import messages
from wave_mapping_msgs.msg import AircraftState
from wave_mapping_msgs.msg import EnergyState
from wave_mapping_msgs.msg import SoaringMap
from wave_mapping_msgs.msg import SplineModel
# import nat'e energy map
from soaring_maps.energy_map import EnergyMap
# import the splines
from geometry.bSpline import TensorProductSpline
# import some utilities
from geodedics.geodedics import lla2ned
# import required stuff to read parameters file
import xml.etree.ElementTree as et
import os

class ThermalBSplineModel:
    """ Models thermals using tensor product bsplines

    """
    mission_center = np.zeros((3,))
    last_sector = (2,2)
    sector_size = 0
    n_sector_knots = 0
    rate = 1
    model = TensorProductSpline()
    ac_lla = np.zeros((3,))
    boundary_conditions = None
    P = np.identity(2)
    #energy_map = EnergyMap()

    def __init__(self, param_file=None, rate=1, pub_model=None):
        """ constructor for the thermal model
        
        Starts things up and gives it knowledge of the configuration
        
        Args:
            param_file: xml paramter file
            
        Returns:
            the thermal model object
        """
        
        if param_file is None:
            self.mission_center = np.zeros((3,))
            self.energy_map = None
        else:
            dir_name =  os.path.dirname(__file__)
            map_file = os.path.join(
                dir_name, '../../../soaring_parameters/map.mat')
            matlab_map_data = sio.loadmat( map_file )
            self.mission_center = np.asarray([
                np.squeeze(matlab_map_data['lat_center']),
                np.squeeze(matlab_map_data['lon_center']),
                np.mean(matlab_map_data['elevation'])])
            self.energy_map = EnergyMap(1.0, matlab_map_data)
                
        self.rate = rate
        
        self.pub = pub_model
        
        k_i = np.arange(-75,80,10)
        k_b = np.vstack((np.arange(-105,-75,10), np.arange(85,115,10)))
        order = 3
        coeffs = np.zeros((pow(len(k_i) + order - 1,2),))
        
        self.sector_size = 50
        self.n_sector_knots = 5
        self.model = TensorProductSpline(
            (k_i, copy.deepcopy(k_i)),
            (order, copy.deepcopy(order)),
            (k_b, copy.deepcopy(k_b)),
            coeffs)
            
        self.P = sps.identity(len(coeffs))*1
        self.Q = sps.identity(len(coeffs))*0.1/pow(self.rate,2.0)
        self.R = 1.0        
        
    def time_update(self, stuff=None):
        """ Time update for kalman filter
        
        Propagate model forward one timestep
        
        Args: None
        
        Returns: None
        
        Shrinks lift with an assumed e-folding time of 60 seconds, 
        covariance in the model is generated at the rate of Q (a model
        property) Boundary conditions not yet implemented.
        """
        
        self.model.coords *= np.exp(-1/(60*self.rate))
        self.P = self.P + self.Q
        self.apply_boundary_conditions()
        
    def measurement_update(self, energy_state):
        """ Measurement update from energy state
        
        A callback, uses energy state message to update the model with
        the latest observations.
        
        Args:
            energy_state - an energy state message containing data
                            about the aircraft position and energy
        
        Returns: None
        
        Right now the measurement uncertainty is assumed to be R, a
        model parameter that must be specified.
        """
        
        if not energy_state.is_valid:
            rospy.loginfo('invalid energy')
            return
        self.ac_lla = np.zeros((3,))
        self.ac_lla[0] = np.squeeze(np.asarray(energy_state.meas_lla.x))
        self.ac_lla[1] = np.squeeze(np.asarray(energy_state.meas_lla.y))
        self.ac_lla[2] = np.squeeze(np.asarray(energy_state.meas_lla.z))
        self.ac_ned = np.squeeze(np.asarray(
            lla2ned([self.ac_lla], self.mission_center)))
        if not self.energy_map.is_in_geofence(
            self.ac_lla[0], self.ac_lla[1]):
            rospy.loginfo('out of map')
            return
        transition = self.check_sector()
        if (transition != (0,0)):
            self.repanel(transition)
        H = self.model.basis((self.ac_ned[0], self.ac_ned[1]))
        HS = sps.csr_matrix(H)
        HP = HS*self.P
        HPH = (HP*H.T)[0]
        K = HP.T / (HPH + self.R)
        Y = energy_state.ed - energy_state.ed_polar
        Z = Y - self.model.eval((self.ac_ned[0], self.ac_ned[1]))
        self.model.coords += np.squeeze(np.asarray(K.todense()))*Z
        M = sps.identity(len(self.model.coords)) - K*HS
        self.P = M*self.P
        if self.pub is not None:
            spline_msg = SplineModel()
            spline_msg.header.stamp = rospy.Time.now()
            spline_msg.knots_x = self.model.splines[0].knots
            spline_msg.knots_y = self.model.splines[1].knots
            spline_msg.order_x = self.model.splines[0].order
            spline_msg.order_y = self.model.splines[1].order
            spline_msg.c = np.ravel(self.model.coords).tolist()
            self.pub.publish(spline_msg)
            
        
    def check_sector(self):
        """ check where in the map we are
        
        Right now the map is composed of nine "panels," the model tries
        to keep the aircraft in the center panel and will repanel if
        the aircraft flies into one of the outer panels. Current panel
        size is hard-coded as 50m, with 10ish knots per panel.
        
        Args: None (aircraft position information is held in the object
            and updated in the measurement update step)
            
        Returns:
            sector, a tuple specifying the location of the aircraft.
            there is one entry for each dimension of the model. For the
            current, 2d model this means that the center panel is (0,0)
            
        Sector Arrangement:
                 _-1___0___1__
                1|___|___|___|
                0|___|___|___|
               -1|___|___|___|
        """
        # figure out which sector we're in
        map_center = np.asarray(
            [x.interior_knots[7] for x in self.model.splines])
        sector_width = np.asarray(
            [(x.interior_knots[-1]-x.interior_knots[0])/3.0 
            for x in self.model.splines])
        sector = tuple(np.floor((self.ac_ned[0:2] - map_center
            + np.asarray([25,25]))/sector_width))
        return sector
        
    def repanel(self, transition):
        """ move the knots to keep the aircraft in the central sector
        
        If the aircraft flies from the central sector to another one
        we want to repanel, creating a new sector ahead of us and
        forgetting the most distant sector behind us. This will keep
        the aircraft in the middle of our map.
        
        Args: None (aircraft position is stored in the object)
        
        Returns: None
        """
        rospy.loginfo('repanel')
        rospy.loginfo(transition)
        delta = (np.cumprod(np.asarray(self.model.shape))
            /self.model.shape[0])
        for i,sector,spl in zip(
            range(0,len(transition)), transition, self.model.splines):
            if sector > 0:
                spl.interior_knots += self.sector_size
                spl.boundary_knots += self.sector_size
                spl.knots += self.sector_size
                if i == 0:
                    for j in range(0, self.model.shape[i]):
                        for k in range(0,12):
                            self.model.coords[j+k*18] = (
                                self.model.coords[j+(k+6)*18])
                        for k in range(12,18):
                            self.model.coords[j+k*6] = 0.0
                else:
                    for j in range(0,self.model.shape[i]):
                        for k in range(0,12):
                            self.model.coords[j*delta[i]+k] = (
                                self.model.coords[j*delta[i]+k+6])
                        for k in range(12,18):
                            self.model.coords[j*delta[i]+k] = 0.0
            elif sector < 0:
                spl.interior_knots -= self.sector_size
                spl.boundary_knots -= self.sector_size
                spl.knots -= self.sector_size
                #for j in range(0,self.model.shape[i]):
                #    for k in range(0,22):
                #        self.model.coords[j*delta[i]+k+11] = (
                #            self.model.coords[j*delta[i]+k])
                #    for k in range(0,11):
                #        self.model.coords[j*delta[i]+k] = 0.0
                if i == 0:
                    for j in range(0, self.model.shape[i]):
                        for k in range(17,5,-1):
                            self.model.coords[j+(k)*18] = (
                                self.model.coords[j+(k-6)*18])
                        for k in range(0,6):
                            self.model.coords[j+k*18] = 0.0
                else:
                    for j in range(0,self.model.shape[i]):
                        for k in range(17,5,-1):
                            self.model.coords[j*delta[i]+k] = (
                                self.model.coords[j*delta[i]+k-6])
                        for k in range(0,6):
                            self.model.coords[j*delta[i]+k] = 0.0
                    
        
    def update_boundary_conditions(self, energy_map):
        """ update BC info
        
        A callback, fires when energy map is published. We'll save the
        energy map so that we'll have it as boundary conditions on our
        model.
        
        Args:
            energy_map, an energy map message containing the strength
                of the lift in each cell of the map
        
        Returns: None
        """
        self.energy_map.wz_mean = energy_map.data
        
    def apply_boundary_conditions(self):
        """ enforce BC's
        
        When we run the update step, we want to ensure that our map 
        boundary matches that of the lower-resolution map that Nate is
        building. This will apply the saved energy map as BC's for this
        model.
        
        Args: None
        
        Returns: None
        """
        
        do_stuff = None
