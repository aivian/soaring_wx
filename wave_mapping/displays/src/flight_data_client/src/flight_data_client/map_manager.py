import sys
import os
from PyQt4 import QtGui, QtCore
import time
import copy

import rospy
import rospkg
import scipy.io as sio
import numpy as np
from math import sqrt

#from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import pyqtgraph.opengl as pgl
from pyqtgraph.dockarea import *

# import protocol
from communications.AVIa_tcp import AVIaTCPParser as AVIaParser
from communications.socket_connections import pysock

# import serial packets
from communications.AVIa_packet import AVIaPacket
from wave_mapping_packets.AircraftStatePacket import AircraftStatePacket
from wave_mapping_packets.SplineModelPacket import SplineModelPacket

# import widgets
import graphics.plot_widgets as plot_widgets
from graphics.glider3d import glider3d

# more stuff
from geometry.quaternion import Quaternion
from geometry.bSpline import TensorProductSpline

# map things
from soaring_maps.energy_map import EnergyMap
from soaring_maps.lift_map import LiftMap

# things I've filed away in the src folder of this package for manageability
from flight_data_client.connection_manager import SocketServicer

import pdb

class SplineManager(QtGui.QWidget):
    def __init__(self):
        super(SplineManager, self).__init__()
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
        
        self.N = self.model.basis((np.linspace(-75,75,50),np.linspace(-75,75,50)))
        
        self.spl_viz = None
    
    def register_visualization(self, spl_viz=None):
        if spl_viz is None:
            return
        
        self.spl_viz = spl_viz
        self.spl_viz.update_data(x_new = np.linspace(-75,75,50),
            y_new = np.linspace(-75,75,50), z_new = np.zeros((50,50)))
            
    def new_data(self, c):
        w = np.dot(self.N, c);
        w = np.reshape(w, (50,50))
        if self.spl_viz is None:
            return
        self.spl_viz.update_sfc_data(w)

class SplineImgManager(QtGui.QWidget):
    def __init__(self):
        super(SplineImgManager, self).__init__()
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
        
        self.N = self.model.basis((np.linspace(-75,75,50),np.linspace(-75,75,50)))
        
        self.spl_img_display = None
    
    def register_visualization(self, spl_img_display=None):
        if spl_img_display is None:
            return
        
        self.spl_img_display = spl_img_display
        self.spl_img_display.update_data(np.zeros((50,50)))
            
    def update_data(self, c):
        w = np.dot(self.N, c);
        w = np.reshape(w, (50,50))
        if self.spl_img_display is None:
            return
        self.spl_img_display.update_data(np.flipud(w))
        
class MapManager(QtGui.QWidget):
    
    def __init__(self, dt):
        super(MapManager, self).__init__()
        """ set up the soaring maps
        """
        # rospack stuff for getting paths
        rospack = rospkg.RosPack()
        base_path = rospack.get_path('soaring_maps')
        map_file = os.path.join(base_path, '../soaring_parameters/map.mat')
        matlab_map_data = sio.loadmat(map_file)
        # initialize map
        self.map = EnergyMap(dt, matlab_map_data)
        # save time of last map update..
        self.last_update = rospy.Time.now()
        self.img_viewer = None
        self.img_display = None
        self.cmin = 0.
        self.cmax = 255.
            
    def time_update(self):
        # simple method calls the time update of the energy map object
        self.map.time_update()
    
    def measurement_update(self, now, wind, lla):
        # do a measurement update of the map
        self.map.measurement_update(lla[0], lla[1], wind[2])
        # save time of last map update..
        self.last_update = rospy.Time.now()
        img = self.colorize(np.array(self.map.output_energy_img()))
        if self.img_viewer is not None:
            self.img_viewer.setImage(img)
            
    def colorize(self, img):
        nimg = sqrt(img.shape[0])
        img = np.reshape(img,(nimg,nimg))
        c = np.zeros((nimg, nimg, 3))
        crange = self.cmax - self.cmin
        c[:,:,0] = (img - self.cmin)/crange
        c[:,:,1] = (1 - abs(img/crange))
        c[:,:,2] = -(img - self.cmax)/crange
        return c
