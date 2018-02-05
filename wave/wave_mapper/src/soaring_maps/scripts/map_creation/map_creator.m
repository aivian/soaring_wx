% this script expects that map data has already been imported, it will
% create and save the relevant data to initialize the ROS mapping system

classdef map_creator < handle
    properties
        lat_0 % rad
        lon_0 % rad
        Lx
        Ly
        dx
        ned
        ned_lon_vec
        ned_lat_vec
        nlcd_code
        nlcd_lon_vec
        nlcd_lat_vec
        Nx
        Ny
        x_lim
        y_lim
        x_center
        y_center
        x_corner
        y_corner
        lat_center
        lon_center
        lat_corner
        lon_corner
        elevation_centers
        elevation_corners
        land_cover 
        terrain_normal
        land_cover_img
        orthoimage_img
        lat_lim
        lon_lim
    end
    
    methods
        function obj = map_creator(lat0,lon0,Lx,Ly,dx,filename)
            % specify a map center (lat0, lon0) (decimal deg), a map size (Lx, Ly) (m), 
            % and a cell size (dx) (m)
            obj.lat_0 = lat0*pi/180;    % decimal degrees of map center  
            obj.lon_0 = lon0*pi/180;   % decimal degrees of map center
            obj.Lx = Lx;
            obj.Ly = Ly;
            obj.dx = dx;

            % from lat0, lon0 find latlim and lonlim

            % use latlim and lonlim to download map data

            if ~exist(filename, 'file')
                % geodata region to download
                scale_f = 2.0;
                [latmin, lonmin] = obj.loc2geo( scale_f*-obj.Lx, scale_f*-obj.Ly );
                [latmax, lonmax] = obj.loc2geo( scale_f*obj.Lx, scale_f*obj.Ly );
                Latlim = 180/pi*[latmin, latmax];
                Lonlim = 180/pi*[lonmin, lonmax];
                % obtain and save data
                download_geo_data(Latlim, Lonlim, (2*scale_f*obj.Lx),(2*scale_f*obj.Ly),filename)
            end
            % parse saved map_data file
            load(filename);
            obj.ned = double(A_ned);
            obj.ned_lon_vec = lonlim(1)+[0.5:1:size(A_ned,2)-0.5]*R_ned(2,1); % lon vec is increasing from lon_min
            obj.ned_lat_vec = latlim(2)+[0.5:1:size(A_ned,1)-0.5]*R_ned(1,2); % lat vec is decreasing from lat_max
            obj.nlcd_code = nlcd_parse( A_lc );
            obj.nlcd_lon_vec = lonlim(1)+[0.5:1:size(A_lc,2)-0.5]*R_lc(2,1); % lon vec is increasing from lon_min
            obj.nlcd_lat_vec = latlim(2)+[0.5:1:size(A_lc,1)-0.5]*R_lc(1,2); % lat vec is decreasing from lat_max
            % hold onto images for plotting
            obj.land_cover_img = A2_lc;
            obj.orthoimage_img = A2_ort;
            % hold onto lat an lon lim of the geo data
            obj.lat_lim = latlim;
            obj.lon_lim = lonlim;

            % find the number of cells in each direction
            obj.Nx = floor(obj.Lx/obj.dx);     %number of cells in x-dir (north)
            obj.Ny = floor(obj.Ly/obj.dx);     %number of cells in y-dir (east)

            % find map limits in ned coordinate system
            obj.x_lim = [ -(obj.Nx*obj.dx)/2, (obj.Nx*obj.dx)/2 ];
            obj.y_lim = [ -(obj.Ny*obj.dx)/2, (obj.Ny*obj.dx)/2 ];

            % determine the locations of cell centers in NED
            obj.x_center = (obj.x_lim(1)+obj.dx/2):obj.dx:(obj.x_lim(2)-obj.dx/2);
            obj.y_center = (obj.y_lim(1)+obj.dx/2):obj.dx:(obj.y_lim(2)-obj.dx/2);

            % determine the locations of cell border intersections in NED
            obj.x_corner = obj.x_lim(1):obj.dx:obj.x_lim(2);
            obj.y_corner = obj.y_lim(1):obj.dx:obj.y_lim(2);

            % determine corresponding lat-lon arrays
            [obj.lat_center, obj.lon_center] = obj.loc2geo(obj.x_center,obj.y_center);
            [obj.lat_corner, obj.lon_corner] = obj.loc2geo(obj.x_corner,obj.y_corner);

            % make an n-d grid of the lat-lon arrays
            [lat_center_index, lon_center_index] = ndgrid(obj.lat_center, obj.lon_center);
            [lat_corner_index, lon_corner_index] = ndgrid(obj.lat_corner, obj.lon_corner);

            % interp2 (linear) the n-d grids to find elevations of cell centers and of
            % cell border intersections
            obj.elevation_centers = obj.elevation(lat_center_index, lon_center_index);
            obj.elevation_corners = obj.elevation(lat_corner_index, lon_corner_index);

            % interp2 (nearest) the n-d grid of centers to find land-cover
            obj.land_cover = obj.land_cover_code( lat_center_index, lon_center_index );

            % find array of unit normals
            obj.terrain_normal = zeros(obj.Nx,obj.Ny,3);
            for i = 1:obj.Nx
                for j = 1:obj.Ny
                    % find points A,B,C,D as defined in K.Cheng 2013 fig
                    % 2.4, note that elevation is a negative z coordinate
                    A = [obj.x_corner(i), obj.y_corner(j), -1.0*obj.elevation_corners(i,j)];
                    B = [obj.x_corner(i+1), obj.y_corner(j+1), -1.0*obj.elevation_corners(i+1,j+1)];
                    C = [obj.x_corner(i), obj.y_corner(j+1), -1.0*obj.elevation_corners(i,j+1)];
                    D = [obj.x_corner(i+1), obj.y_corner(j), -1.0*obj.elevation_corners(i+1,j)];
                    AB = A-B;
                    CD = C-D;
                    normal_vec = cross(AB,CD);
                    obj.terrain_normal(i,j,:) = normal_vec./(norm(normal_vec));
                end
            end
        end
        
        function [] = save_map_parameters( obj, filename )
            % method used to save map aparameters to a file for the ROS
            % system
            % elevation map
            elevation = obj.elevation_centers;
            save(filename, 'elevation');
            % nx*xy*3 array containing NED terrain normal vectors
            terrain_normal = obj.terrain_normal;
            save(filename, 'terrain_normal','-append');
            % array containing nlcd land cover codes
            land_cover = obj.land_cover;
            save(filename, 'land_cover','-append');
            % arrays containing the lat/lon corrdinated of the cell centers
            lat_mid = obj.lat_center; lon_mid = obj.lon_center; dx = obj.dx;
            save(filename, 'lat_mid','lon_mid','dx','-append');
            % the lat/lon coordinated of the center of the map
            lon_center = obj.lon_0; lat_center = obj.lat_0;
            save(filename, 'lat_center','lon_center','-append');
        end
        
        function [geoLat, geoLon] = loc2geo( obj, locX, locY )
            % LOC2GEO takes x and y points in the local (map) coordinate
            % frame and returns the corresponding latitudes and longitudes
            
            [dlat,dlon] = lg2dg( locX, locY, obj.lat_0, 0.0 );
            
            geoLat = obj.lat_0 + dlat; % (rad)
            geoLon = obj.lon_0 + dlon; % (rad)
        end
        
        function [locX, locY] = geo2loc( obj, geoLat, geoLon )
            % GEO2LOC takes lat and lon geodetic coordinates and returns
            % their corresponding x and y coordinates in the local (map)
            % system
            
            % geoLat and geoLon should be in radians!!!!!
            
            dlat = geoLat - obj.lat_0;
            dlon = geoLon - obj.lon_0;
            
            [dx,dy] = dg2lg( dlat, dlon, obj.lat_0, 0.0 );
            
            locX = 0 + dx;
            locY = 0 + dy;
        end

        function [zg] = elevation( obj, lat, lon, varargin )
            % lat is 
            lat = lat*180/pi;
            lon = lon*180/pi;

            zg = interp2( obj.ned_lon_vec, obj.ned_lat_vec, obj.ned, lon, lat, 'spline');          
        end % ground elevation

        function [zg] = land_cover_code( obj, lat, lon, varargin)
            % lat is 
            lat = lat*180/pi;
            lon = lon*180/pi;

            zg = interp2( obj.nlcd_lon_vec, obj.nlcd_lat_vec, obj.nlcd_code, lon,lat, 'nearest' );          
        end % nation land cover database code
    end
end

function [dlat,dlon]=lg2dg(dx,dy,lat,h,a,e2)
    % LG2DG  Converts local geodetic coordinates to Ælat,Ælon,Æh.
    %   Local origin at lat,lon.  If astronomic lat,h input,
    %   then output is in local astronomic system.  Vectorized.
    %   See also DG2LG.
    % Version: 2011-02-19
    % Useage:  [dlat,dlon]=lg2dg(dx,dy,lat,h,a,e2)
    %          [dlat,dlon]=lg2dg(dx,dy,lat,h)
    % Input:   dx   - vector of x (N) coordinates in local system
    %          dy   - vector of y (E) coordinates in local system
    %          lat  - vector of lats of local system origins (rad)
    %          h    - vector of hts of local system origins
    %          a    - ref. ellipsoid major semi-axis (m); default GRS80
    %          e2   - ref. ellipsoid eccentricity squared; default GRS80
    % Output:  dlat - vector of latitude differences (rad)
    %          dlon - vector of longitude differences (rad)

    % Copyright (c) 2011, Michael R. Craymer
    % All rights reserved.
    % Email: mike@craymer.com

    if nargin ~= 4 & nargin ~= 6
      warning('Incorrect number of input arguments');
      return
    end
    if nargin == 4
      [a,b,e2]=refell('wgs84');
    end

    v=a./sqrt(1-e2.*sin(lat).^2);
    r=v.*(1-e2)./(1-e2.*sin(lat).^2);
    dlat=dx./(r+h);
    dlon=dy./cos(lat)./(v+h);
end

function [a,b,e2,finv]=refell(type)
    % REFELL  Computes reference ellispoid parameters.
    %   Reference: Department of Defence World Geodetic
    %   System 1984, DMA TR 8350.2, 1987.
    %   TOPEX Reference: <http://topex-www.jpl.nasa.gov/aviso/text/general/news/hdbk311.htm#CH3.3>
    % Version: 1 Jun 04
    % Useage:  [a,b,e2,finv]=refell(type)
    % Input:   type - reference ellipsoid type (char)
    %                 CLK66 = Clarke 1866
    %                 GRS67 = Geodetic Reference System 1967
    %                 GRS80 = Geodetic Reference System 1980
    %                 WGS72 = World Geodetic System 1972
    %                 WGS84 = World Geodetic System 1984
    %                 ATS77 = Quasi-earth centred ellipsoid for ATS77
    %                 NAD27 = North American Datum 1927 (=CLK66)
    %                 NAD83 = North American Datum 1927 (=GRS80)
    %                 INTER = International
    %                 KRASS = Krassovsky (USSR)
    %                 MAIRY = Modified Airy (Ireland 1965/1975)
    %                 TOPEX = TOPEX/POSEIDON ellipsoid
    % Output:  a    - major semi-axis (m)
    %          b    - minor semi-axis (m)
    %          e2   - eccentricity squared
    %          finv - inverse of flattening

    % Copyright (c) 2011, Michael R. Craymer
    % All rights reserved.
    % Email: mike@craymer.com

    type=upper(type);
    if (type=='CLK66' | type=='NAD27')
      a=6378206.4;
      finv=294.9786982;
    elseif type=='GRS67'
      a=6378160.0;
      finv=298.247167427;
    elseif (type=='GRS80' | type=='NAD83')
      a=6378137.0;
      finv=298.257222101;
    elseif (type=='WGS72')
      a=6378135.0;
      finv=298.26;
    elseif (type=='WGS84')
      a=6378137.0;
      finv=298.257223563;
    elseif type=='ATS77'
      a=6378135.0;
      finv=298.257;
    elseif type=='KRASS'
      a=6378245.0;
      finv=298.3;
    elseif type=='INTER'
      a=6378388.0;
      finv=297.0;
    elseif type=='MAIRY'
      a=6377340.189;
      finv=299.3249646;
    elseif type=='TOPEX'
      a=6378136.3;
      finv=298.257;
    end
    f=1/finv;
    b=a*(1-f);
    e2=1-(1-f)^2;
end

function [codes] = nlcd_parse( image )
    % funciton intended to parse National Land Cover Dataset image into
    % codes, the legend can be found here: http://www.mrlc.gov/nlcd11_leg.php
    codes = zeros(size(image(:,:,1)));
    for i = 1:size(image,1)
        for j = 1:size(image,2)
            color = squeeze(image(i,j,:));
            if color == [71;107;160]
                % open water
                codes(i,j) = 11; 
            elseif color == [221;201;201]
                % developed, open space
                codes(i,j) = 21;
            elseif color == [216;147;130]
                % developed, low intensity
                codes(i,j) = 22;
            elseif color == [237;0;0]
                % developed, medium intensity
                codes(i,j) = 23;
            elseif color == [170;0;0]
                % developed, high intensity
                codes(i,j) = 24;
            elseif color == [178;173;163]
                % barren land (rock/sand/clay
                codes(i,j) = 31;
            elseif color == [104;170;99]
                % deciduous forest
                codes(i,j) = 41;
            elseif color == [28;99;48]
                % evergreen forest
                codes(i,j) = 42;
            elseif color == [181;201;142]
                % mixed forest
                codes(i,j) = 43;
            elseif color == [170;112;40]
                % dwarf scrub
                codes(i,j) = 51;
            elseif color == [219;216;60]
                % pasture/hay
                codes(i,j) = 81;
            elseif color == [112;163;186]
                % emergent herbaceous wetlands
                codes(i,j) = 95; 
            end
        end
    end
end