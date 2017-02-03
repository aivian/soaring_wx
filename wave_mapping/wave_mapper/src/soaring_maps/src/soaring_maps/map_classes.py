# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 13:47:05 2014

@author: Nate Depenbusch, AVIA Lab
"""

# to add to map data:
# 3(n_latxn_lon) arrays giving components of the unit normal vector in each cell
# The lat,lon,alt of the center of the map

# needs documentation on creating the map initialization info

import numpy as np
#import matplotlib as mpl
import matplotlib.path as mplpath
import Pysolar as solar
import fastkml as kml
import scipy.ndimage as ndimage
from scipy.stats import norm

class Map(object):
    # base map class
    def __init__( self, data ):
        self.lat_mid = data['lat_mid']
        self.lon_mid = data['lon_mid']
        self.lat_n = np.size(self.lat_mid)
        self.lon_n = np.size(self.lon_mid)
        self.resolution = data['dx']
        self.elevation = data['elevation']
        # start with a simple geofence around the whole region [lon, lat]
        posts = [[self.lon_mid[0,0], self.lat_mid[0,0]],
                 [self.lon_mid[0,0], self.lat_mid[0,-1]],
                [self.lon_mid[0,-1], self.lat_mid[0,-1]],
                [self.lon_mid[0,-1], self.lat_mid[0,0]]]
        # use the set method to set the fence   
        self._set_geofence( posts )
        # save coordinates of the map center
        self.lat_center = data['lat_center'] #self.lat_mid[0, np.round(self.lat_n / 2.0)]
        self.lon_center = data['lon_center'] #self.lon_mid[0, np.round(self.lon_n / 2.0)]
        
    def is_in_geofence( self, lat, lon ):
        # method to check if provided lat, lon is inside geofence
        return self.bounds.contains_point( (lat, lon) )
    
    def import_geofence( self, filename ):
        # method takes a filename and reads the first feature in as a geofence, then sets the fence-posts
        doc = file( filename ).read()
        k = kml.KML()
        k.from_string(doc)
        features = list(k.features())
        f = list(features[0].features())
        geo = f[0].geometry
        self._set_geofence( np.array(geo.exterior.coords) )
    
    def _set_geofence( self, fence_posts ):
        # method for setting the geo-fence boundary given a 2d list of posts
        # [[lon0,lat0,alt0],[lon1,lat1,alt1],...[lonf,latf,altf]]
        # this is the order Google Earth kml files come in, altitude of posts is ignored
        # the posts can be closed or not, radians or degrees, either as a list or np.array (in above format)
        #convert the supplied list to a np.array (this doesn't do anything if it's already an array)
        fence_posts = np.array(fence_posts)
        # check that the region is closed
        if (not np.array_equal( fence_posts[0], fence_posts[-1] )):
            # first, close the boundary by adding the first post again to the end
            fence_posts = np.vstack( (fence_posts,fence_posts[0,:]) )
        # check that the supplied posts are in radians
        if (( np.amax(fence_posts) > np.pi ) or ( np.amin(fence_posts) < -1.0*np.pi )):
            fence_posts = np.deg2rad(fence_posts)
        # find the number of posts (including the doubled first post)
        n_posts = fence_posts.shape[0]
        # declare an empty list for the codes
        fence_codes = [mplpath.Path.LINETO]*n_posts
        # make the first code a "start" and the last code a "close"
        fence_codes[0] = mplpath.Path.MOVETO
        fence_codes[-1] = mplpath.Path.CLOSEPOLY
        # set the Path object that defines the lat, lon bounds of the region
        # the bounds are [lat, lon] order and do not include altitude, thus the flipping and limiting
        self.bounds = mplpath.Path( np.fliplr(fence_posts[:,0:2]), fence_codes )
        
    def current_cell( self, lat, lon ):
        # NOTE, THIS SHOULD RETURN THE NEAREST VALID CELL
        # returns the indices of the cell closes to the provided lat and lon
        # this uses the spherical law of cosines rather than the Haversine formula
        # http://www.movable-type.co.uk/scripts/latlong.html
        lon_index = np.tile(self.lon_mid, (self.lat_n,1))
        lat_index = np.tile(self.lat_mid.T, (1,self.lon_n))
        delta_lon = (lon - lon_index)
        R = 6371000.0
        d = np.arccos(np.sin(lat_index)*np.sin(lat) + np.cos(lat_index)*np.cos(lat)*np.cos(delta_lon))*R
        # find minimum distance
        i_min = d.argmin(0)[0]
        j_min = d.argmin(1)[0]
        return np.array([i_min,j_min])
        
    def _index_range( self, cur_idx, spread ):
        # internal method for safely finding the range of indices from idx
        min_i = int(max( cur_idx[0] - spread , 0 )) 
        max_i = int(min( cur_idx[0] + spread , self.lat_n-1 ))
        min_j = int(max( cur_idx[1] - spread , 0 ))  
        max_j = int(min( cur_idx[1] + spread , self.lon_n-1 ))
        return [min_i, max_i, min_j, max_j]
        
    def _output_flat( self, map_data ):
        # method for flattening the map before publishing
        # convert the numpy array to a row-major flat array, then un-numpy it
        return map_data.ravel().tolist()
        
    def _get_distance_fast( self, X0, X1 ):
        lat0 = X0[0]
        lon0 = X0[1]
        lat1 = X1[0]
        lon1 = X1[1]
        delta_lon = lon1 - lon0
        R = 6371000.0
        return np.arccos(np.sin(lat1)*np.sin(lat0) + np.cos(lat1)*np.cos(lat0)*np.cos(delta_lon))*R
        
    def _get_bearing( self, X0, X1 ):
        lat0 = X0[0]
        lon0 = X0[1]
        lat1 = X1[0]
        lon1 = X1[1]
        delta_lon = lon1 - lon0
        # compute bearing to center of cell i,j
        return np.arctan2( (np.sin(delta_lon)*np.cos(lat1)), (np.cos(lat0)*np.sin(lat1) - np.sin(lat0)*np.cos(lat1)*np.cos(delta_lon)) )
        
    def lat_lon( self, idx_i, idx_j):
        # method retuns the lat lon of a set of indices
        return [ self.lat_mid[0,idx_i], self.lon_mid[0,idx_j] ]  
        
    def get_elevation( self, lat, lon ):
        # method returns the current agl height of the aircraft given a lat, lon
        idx = self.current_cell( lat, lon )
        return self.elevation[idx[0],idx[1]]

        
