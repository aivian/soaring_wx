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

class WindMap(Map):
    def __init__( self, data ): 
        # call the base map class
        super( WindMap, self ).__init__(data)
        self.wind_speed = 0.0
        self.wind_direction = 0.0   # direction the wind is COMING FROM, this would be the value read out by a 
        self.upwind = np.zeros((self.lat_n,self.lon_n), dtype=float )
        # find the bearing to each cell from the center
        self.bearing_array = np.zeros((self.lat_n,self.lon_n), dtype=float )
        self.distance_array = np.zeros((self.lat_n,self.lon_n), dtype=float )
        self._set_bearing_array()
        # counter for wind map update
        self.update_cnt = 0
    
    def _set_bearing_array( self ):
        for i in range(0, self.lat_n):
            for j in range(0, self.lon_n):
                self.bearing_array[i,j] =  self._get_bearing( [self.lat_center, self.lon_center], [self.lat_mid[0,i], self.lon_mid[0,j]] )
                self.distance_array[i,j] =  self._get_distance_fast( [self.lat_center, self.lon_center], [self.lat_mid[0,i], self.lon_mid[0,j]] )
        
    def set_wind( self, wind_from ):
        # method takes current measured wind vecotr and incorporates it with prior knowledge
        wind_from = np.array( wind_from )
        self.wind_direction = np.arctan2( -1.0*wind_from[1], -1.0*wind_from[0] )
        self.wind_speed = np.sqrt( wind_from.dot( wind_from ) ) # faster than linalg.norm
        # this doesn't need to be done on every message, maybe every 5th? (2 Hz)
        if (self.update_cnt == 50):
            self.update_cnt = 0
            for i in range(0, self.lat_n):
                for j in range(0, self.lon_n):
                    bearing_vec = np.array([ np.cos( self.bearing_array[i,j] ), np.sin( self.bearing_array[i,j] ) ]) * self.distance_array[i,j]
                    self.upwind[i,j] = np.inner( -1.0*wind_from, bearing_vec ) 
        self.update_cnt = self.update_cnt + 1
    
    def output_priority_flat( self ):
        # compute the priority due to a cell based on the current wind
        upwind_min =  np.amin(self.upwind)
        upwind_range = np.amax(self.upwind) - upwind_min
        if upwind_range > 0.0:
            output = ((self.upwind - upwind_min) / upwind_range )**( self.wind_speed / 10.0 )
        else:
            output = np.ones((self.lat_n,self.lon_n))
        return self._output_flat( output )
        
