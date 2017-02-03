# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 13:47:05 2014

@author: Nate Depenbusch, AVIA Lab
"""

# to add to map data:
# 3(n_latxn_lon) arrays giving components of the unit normal vector in each cell
# The lat,lon,alt of the center of the map

import numpy as np
import matplotlib.path as mplpath
import scipy.ndimage as ndimage
from scipy.stats import norm
import fastkml as kml
from soaring_maps.map_classes import Map

class EnergyMap(Map):
    def __init__( self, dt, data ):
        # call the base map class
        super( EnergyMap, self ).__init__(data)
        # declare some things only used in this map
        self.C = np.array([1])
        self.var_max = 4.0 # sigma_wz squared
        self.meas_noise = 0.5
        self.meas_rad = 50.0
        self.expected_lifespan = 15*60
        # note that edot and wz should be equivalent if polar, motor, and external contributions are correctly removed
        self.wz_mean = np.zeros((self.lat_n,self.lon_n), dtype=float )
        self.wz_var = self.var_max*np.ones((self.lat_n,self.lon_n), dtype=float )
        # lastly, tune the filter parameters
        self._tune_filter(dt)
        
    def measurement_update( self, lat, lon, edot ):
        # method does a measurement update to the map
        # first find the current cell
        cur_idx = self.current_cell( lat, lon )
        # we will assign the measurement in edot to neighboring cells, but increase the measurement noise
        spread = 2.0*np.ceil(self.meas_rad/self.resolution)     
        idxs = self._index_range( cur_idx, spread )
        for i in range(idxs[0], idxs[1]):
            for j in range(idxs[2], idxs[3]):
                # find distance between provided lat, lon and cell i,j
                delta_lon = lon - self.lon_mid[0,j]
                R = 6371000.0
                d = np.arccos(np.sin(self.lat_mid[0,i])*np.sin(lat) + np.cos(self.lat_mid[0,i])*np.cos(lat)*np.cos(delta_lon))*R
                # check that the distance is reasonable
                if (d < 1.0*self.meas_rad):
                    # find the measurement noise associated with this distance
                    Q = self._distance_2_noise(d)
                    # compute the Kalman gain
                    K = self.wz_var[i,j]*self.C.T/(self.C*self.wz_var[i,j]*self.C.T + Q)
                    # update the estimate in cell i,j
                    self.wz_mean[i,j] = self.wz_mean[i,j] + K*(edot - self.C*self.wz_mean[i,j])
                    # update the varience in cell i,j
                    self.wz_var[i,j] = (1 - K*self.C)*self.wz_var[i,j]
        
    def _distance_2_noise( self, d ):
        # simple private method to convert a distance to a measurement noise
        if (d <= 1.0*self.meas_rad):
            # if we can conceiveably measure the wind
            decay = 2.0
            return ( self.var_max*(d/self.meas_rad)**decay + self.meas_noise)
        else:
            # if we can't conceivably measure the wind
            return self.var_max
    
    def time_update( self ):
        # method for doing the time update on the map, current just a decay
        self.wz_mean = self.A * self.wz_mean
        self.wz_var = self.A*self.wz_var*self.A.T + self.R
        self.wz_var = np.minimum(self.wz_var, self.var_max)
        
    def _tune_filter( self, dt ):
        # private method for approximately tuning the Kalman filter
        # much of this is "good enough"
        # dt is the time interval at which the time-update will be called (sec)
        self.A = np.array([(self.expected_lifespan/1.4)/((self.expected_lifespan/1.4) + dt)])
        self.R = np.array([3.0 * (self.var_max * dt) / self.expected_lifespan])
                
    def output_uncertainty_flat( self ):
        # method that outputs a flatted and scaled uncertainty map
        return self._output_flat( (self.wz_var/self.var_max) )   
        
    def output_variance_flat( self ):
        # method that outputs a flatted and unscaled variance map
        return self._output_flat( self.wz_var ) 
        
    def output_uncertainty_img( self ):
        # method that outputs a flatted and scaled uncertainty map intended to be shown as an 8bit image
        # the map is flipped such that the image is displayed as expected
        scaled_map = np.uint8(-255.0*(np.flipud(self.wz_var/self.var_max)-1.0))
        return self._output_flat( scaled_map )  

    def output_energy_flat( self ):
        # method that spits out a flattened, unscaled energy map
        return self._output_flat( self.wz_mean )
        
    def output_energy_img( self ):
        # method that outputs a flatted and scaled uncertainty map intended to be shown as an 8bit image
        # the map is flipped such that the image is displayed as expected
        # limit the range first
        fsr = 4.0 # full scale range to be plotted
        energy_map = (0.5*fsr + np.maximum(np.minimum(self.wz_mean, 0.5*fsr), -0.5*fsr))/fsr
        scaled_map = np.uint8(255.0*np.flipud(energy_map))
        return self._output_flat( scaled_map )  
        
