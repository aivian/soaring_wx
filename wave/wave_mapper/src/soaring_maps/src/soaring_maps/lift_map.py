# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 13:47:05 2014

@author: Nate Depenbusch, AVIA Lab
"""

import numpy as np
import matplotlib.path as mplpath
import scipy.ndimage as ndimage
from scipy.stats import norm
import fastkml as kml
from soaring_maps.map_classes import Map

# map analysis class, is handed data from the energy map and does some feature recognition
class LiftMap(Map):
    def __init__( self, data ): 
        # call the base map class
        super( LiftMap, self ).__init__(data)
        # we need a mean and variation here too
        self.wz_mean = np.zeros((self.lat_n,self.lon_n), dtype=float )
        self.wz_var = 1.0*np.ones((self.lat_n,self.lon_n), dtype=float )
        self.P = np.zeros((self.lat_n,self.lon_n), dtype=float )
        # empty feature arrays
        self.num_features = 0
        self.feature_loc = np.zeros((self.num_features,2), dtype=float )
        self.feature_mean = np.zeros((1,self.num_features), dtype=float )
        self.feature_var = 1.0*np.ones((1,self.num_features), dtype=float )
        # some data for the thermal lift template (or kernel)
        avg_rad =  75.0 # minimum radius of a lift sources needed to circle
        min_rad = 20.0
        res = 25.0
        self._set_kernel( avg_rad, min_rad, res)
        self.labeled_array = np.zeros((self.lat_n,self.lon_n))
        # exponential scale factor for distance, 0 implies distance isn't a factor in selecting a lift source
        self.distance_scale = 0.25
        
    def set_variance( self, data ):
        # method takes a flattened variance array and expands and saves it
        data = np.array(data)
        self.wz_var = data.reshape((self.lat_n,self.lon_n))
        
    def set_mean( self, data ):
        # method takes a flattened mean array and expands and saves it
        data = np.array(data)
        self.wz_mean = data.reshape((self.lat_n,self.lon_n))
    
    def _set_kernel( self, out_rad, in_rad, res ):
        # method for setting the kernel (mask) for a minimum size thermal # make a kernel or template for a usable thermal
        outer_radius = max(np.round( out_rad / res ), 1.0)
        inner_radius = np.round( 20.0 / res )
        k = np.zeros((2*outer_radius+1, 2*outer_radius+1))
        y,x = np.ogrid[-outer_radius:outer_radius+1, -outer_radius:outer_radius+1]
        mask = x**2 + y**2 <= outer_radius**2
        mask2 = x**2 + y**2 < inner_radius**2
        k[mask] = 1.0
        k[mask2] = 0.0
        self.kernel = k/np.sum(k)
        
    def set_threshold( self, V_threshold, P_threshold_min = 0.5, P_threshold_max = 0.97 ):
        # method returns a thresholded map
        # convlove the expected wind speed as well
        convolved_mean = ndimage.convolve( self.wz_mean, self.kernel, mode='constant', cval=0.0)
        # find the varience correlated with the convolved value in each cell
        max_unc = np.amax(self.wz_var)
        convolved_var = ndimage.convolve( self.wz_var, self.kernel**2, mode='constant', cval=max_unc)
        # compute the probability that the wz in each cell exceeds V_threshold
        self.P = 1 - norm.cdf( V_threshold, convolved_mean, np.sqrt(convolved_var) )
        # set probability threshold to 2 std's above mean
        P_threshold = min(max( self.P.mean() + 2.0*self.P.std(), P_threshold_min ), P_threshold_max)
        # check the probability versus a threshold
        binary_image = ( self.P >= P_threshold )
        # do two erosions and one dilation to clean up noise
        #binary_image = ndimage.binary_erosion(binary_image, None, 1)
        binary_image = ndimage.binary_dilation(binary_image, None, 3) # this may be excessive, but keep dilation first so thatsingle points aren't lost
        binary_image = ndimage.binary_erosion(binary_image, None, 3)  # this may also be excessive
        # find features and label
        self.labeled_array, self.num_features = ndimage.label( binary_image )
        # allocate correctly sized arrays
        self.feature_loc = np.zeros( (self.num_features, 2) )
        self.feature_var = np.zeros( self.num_features )
        self.feature_mean = np.zeros( self.num_features )
        if (self.num_features > 0):
            # find the grid (pixel) locations of each feature (rounded to the nearest pixel)
            grid_location = np.around(ndimage.maximum_position( convolved_mean, self.labeled_array, index=np.arange(1, self.num_features+1) ))
            for k in range(0, self.num_features):
                # find the lat, lon of each
                 loc = self.lat_lon( grid_location[k][0], grid_location[k][1] )
                 self.feature_loc[k][0] = loc[0]
                 self.feature_loc[k][1] = loc[1]
                 # 
                 self.feature_var[k] = convolved_var[grid_location[k][0], grid_location[k][1]]
                 self.feature_mean[k] = convolved_mean[grid_location[k][0], grid_location[k][1]]
        # return the number of features
        return np.float(self.num_features)
        
    def find_best_lift( self, lat, lon ):
        # method picks the "best" lift given the current lat/lon of the
        if ( self.num_features > 0 ):
            # if there are detected lift sources, allocate a payoff array indicating the attractiveness of each feature
            payoff = np.zeros( self.num_features )
            for k in range( 0, self.num_features ):
                # check if feature k is within the geofence, should be map fence here
                if (self.is_in_geofence( self.feature_loc[k,0], self.feature_loc[k,1] )):
                    # get distance from current aircraft lat/lon to feature k
                    distance = self._get_distance_fast( [lat, lon], [self.feature_loc[k,0], self.feature_loc[k,1]] )
                    # compute payoff of feature k, though anything within 100m has the same divisor
                    payoff[k] = self.feature_mean[k]/(max(distance,100.0)**self.distance_scale)           
            # sort the computed payoffs in order of best to worst
            idx = np.argsort( -1.0*payoff )
            # return the first value in the sorted array indicating the feature with the highest payoff
            return [ self.feature_loc[idx[0],0], self.feature_loc[idx[0],1] ]
        else:
            # return the lat/lon of the center of the map..
            return [ self.lat_center, self.lon_center ]
            
    def output_probability_img( self ):
        # method that outputs a flatted and scaled uncertainty map intended to be shown as an 8bit image
        # the map is flipped such that the image is displayed as expected
        scaled_map = np.uint8( 255.0*np.flipud(self.P) )
        return self._output_flat( scaled_map ) 
        
    def output_lift_flat( self ):
        # method that spits out a flattened, unscaled lift probability map
        return self._output_flat( self.P )
    
