function [] = download_geo_data( latlim, lonlim, lx,ly, filename )
    % retreive info on wms server, for land cover, these two seem to be
    % equivalent
    %[cap,url]=wmsinfo('http://isse.cr.usgs.gov/arcgis/services/LandCover//USGS_EROS_LandCover_NLCD/MapServer/WmsServer?')
    [cap,url]=wmsinfo('http://raster.nationalmap.gov/arcgis/services/LandCover//USGS_EROS_LandCover_NLCD/MapServer/WmsServer?');
    % look at cap.Layer to find the relevant LayerName, here '30'
    land_cover_layer = WMSLayer('ServerURL',url,'LayerName','30');
    land_cover_layer = wmsupdate(land_cover_layer);
    % read in info, for ground cover we want a geotiff and not something
    % rastered
    width = max(512,ceil(ly/30));
    height = max(512,ceil(lx/30));
    [A_lc,R_lc]=wmsread(land_cover_layer,'LatLim',latlim,'LonLim',lonlim,'ImageFormat','image/tiff','ImageWidth',width,'ImageHeight',height);
    % plot to confirm
%     figure
%     geoshow(A_lc,R_lc)

    % NASA has the best elevation data as of the writing of this script
    layers=wmsfind('nasa.network*elev','SearchField','serverurl');
    % peruse layers to find the one we want, I like 'SRTM30 with Bathymetry (900m) merged with global ASTER (30m) and USGS NED (10m)'
    ned_layer = layers(6);
    ned_layer = wmsupdate(ned_layer);
    % here the 'image type' is a bil, which gives raw elevation data which is
    % really nice
    width = max(512,ceil(ly/10));
    height = max(512,ceil(lx/10));
    [A_ned,R_ned]=wmsread(ned_layer,'LatLim',latlim,'LonLim',lonlim,'ImageFormat','image/bil','ImageWidth',width,'ImageHeight',height);
    % plot to confirm
%     figure
%     geoshow(double(A_ned),R_ned,'displaytype','surface')

    % now for pictures! both of these work, the second seems to be black and
    % white
    [cap,url]=wmsinfo('http://raster.nationalmap.gov/arcgis/services/Orthoimagery/USGS_EROS_Ortho_SCALE/ImageServer/WMSServer?');
    % [cap,url]=wmsinfo('http://raster.nationalmap.gov/arcgis/services/Orthoimagery/USGS_EROS_Ortho_1Foot/ImageServer/WMSServer?')
    ortho_layer = WMSLayer('ServerURL',url,'LayerName','0');
    ortho_layer = wmsupdate(ortho_layer);
%     [A_ort,R_ort]=wmsread(ortho_layer,'LatLim',latlim,'LonLim',lonlim,'ImageWidth',3000,'ImageHeight',2000);
%     figure
%     geoshow(A_ort,R_ort)

    % and a little practice draping data over elevation data
    % grab a new land cover map
    [A2_lc,R2_lc]=wmsread(land_cover_layer,'LatLim',latlim,'LonLim',lonlim,'ImageFormat','image/tiff',...
        'imageheight',size(A_ned,1),'imagewidth',size(A_ned,2));
    % drape land cover over elvation
    figure('Renderer','opengl')
    usamap(latlim,lonlim)
    framem off; mlabel off; plabel off; gridm off;
    geoshow(double(A_ned),R_ned,'displaytype','surface','cdata',A2_lc)
    daspectm('m',5) % exaggerate the vertical axis by 5x
    axis vis3d
    % grab some ortho imagery
    [A2_ort,R2_ort]=wmsread(ortho_layer,'LatLim',latlim,'LonLim',lonlim,...
        'imageheight',size(A_ned,1),'imagewidth',size(A_ned,2));
    % drape orto imagery cover over elvation
    figure('Renderer','opengl')
    usamap(latlim,lonlim)
    framem off; mlabel off; plabel off; gridm off;
    geoshow(double(A_ned),R_ned,'displaytype','surface','cdata',A2_ort)
    daspectm('m',5) % exaggerate the vertical axis by 5x
    axis vis3d

    % note, can use geoshow to plot data on the map (flight paths, waypoints,
    % etc.). Can also export as a KML file for plotting in Google Earth.   
    save(filename,'A*','R*','latlim','lonlim')
end
