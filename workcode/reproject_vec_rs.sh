
ogr2ogr -f "ESRI Shapefile" milano_LST_poly_32n.shp milano_LST_poly.shp -s_srs EPSG:4326 -t_srs EPSG:32632 
ogr2ogr -f "ESRI Shapefile" firenze_LST_poly_32n.shp firenze_LST_poly.shp -s_srs EPSG:4326 -t_srs EPSG:32632 
ogr2ogr -f "ESRI Shapefile" bologna_LST_poly_32n.shp bologna_LST_polyn.shp -s_srs EPSG:4326 -t_srs EPSG:32632 
ogr2ogr -f "ESRI Shapefile" firenze_LST_poly_32n.shp firenze_LST_poly.shp -s_srs EPSG:4326 -t_srs EPSG:32632 
ogr2ogr -f "ESRI Shapefile" roma_LST_poly_32n.shp roma_LST_poly.shp -s_srs EPSG:4326 -t_srs EPSG:32632 



gdalwarp   -t_srs EPSG:32632 -r near milano_ispra_csuolo_clip.tif milano_ispra_csuolo_clip_32n.tif
gdalwarp   -t_srs EPSG:32632 -r near bologna_ispra_csuolo_clip.tif bologna_ispra_csuolo_clip_32n.tif
gdalwarp   -t_srs EPSG:32632 -r near firenze_ispra_csuolo_clip.tif firenze_ispra_csuolo_clip_32n.tif
gdalwarp   -t_srs EPSG:32632 -r near palermo_ispra_csuolo_clip.tif palermo_ispra_csuolo_clip_32n.tif
gdalwarp   -t_srs EPSG:32632 -r near Roma_ispra_csuolo_clip.tif roma_ispra_csuolo_clip_32n.tif



