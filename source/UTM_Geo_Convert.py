"""
Python module 'UTM_Geo_Convert.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia, heavily modified from python source code found at
http://stackoverflow.com/questions/343865/how-to-convert-from-utm-to-latlng-in-python-or-javascript/10239676

DISTRIBUTION and USE subject to 'LICENSE_GnuGPLv3.txt' and 'DISCLAIMER.txt' that accompany 
this file. This file is provided FREE OF COST and WITHOUT the author's WARRANTY or LIABILITY 
in any manner whatsoever. This file may be redistributed by the USER as long as the above 
AUTHORSHIP and COPYRIGHT, this statement, and the accompanying LICENSE and DISCLAIMER files 
remain intact. Treat others as you would be treated. Pay it forward. Valar dohaeris.

Send questions, bug reports, any related requests to matt.e.garcia@gmail.com

REFERENCE: If you use this software, please reference the following in your work products:
               Garcia, M., and P.A. Townsend, in review: "Climatological influences on the 
               forest growing season around western Lake Superior, USA." Submitted to J. 
               Geophys. Res. Atmos. on 5 April 2016.
           See also 'README.md', 'CITATION.txt', and 'ACKNOWLEDGEMENTS.txt' for more information.

USAGE: insert 'from UTM_Geo_Convert import *' line near head of script, then (for example)
           UTM_coords = geographic_to_utm(stn_lon, stn_lat, UTMzone)
           stn_easting = int(round(UTM_coords[1],0))
           stn_northing = int(round(UTM_coords[2],0))

PURPOSE: Conversion between UTM and geographic coordinates, or from one UTM zone to another

DEPENDENCIES: The 'osgeo.osr' package, which itself uses 'gdal'

INPUT: Provided by calling script

OUTPUT: Returned to calling script

RUN TIME: negligible
"""

from osgeo import osr  # uses GDAL via osr-to-python bindings

def get_utm_zone(long):    
    zone = int(1 + (long + 180.0) / 6.0)
    return zone

def is_northern(lat):
    if (lat < 0.0):
        return 0
    else:
        return 1

def utm_to_geographic(easting, northing, zone):
    utm_coordinate_system = osr.SpatialReference()
    utm_coordinate_system.SetWellKnownGeogCS("WGS84") # Set unprojected geographic coordinate system
    utm_coordinate_system.SetUTM(int(zone), int(is_northern(northing)))
    geog_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the unprojected geographic coordinate system 
    # create transform component
    utm_to_geog_transform = osr.CoordinateTransformation(utm_coordinate_system, geog_coordinate_system) # (<from>, <to>)
    lon, lat, alt = utm_to_geog_transform.TransformPoint(easting, northing, 0)
    return lon, lat

def geographic_to_utm(lon, lat, zonepref=-1):    
    utm_coordinate_system = osr.SpatialReference()
    utm_coordinate_system.SetWellKnownGeogCS("WGS84") # Set unprojected geographic coordinate system 
    if zonepref != -1:
        zone = zonepref
    else:
        zone = get_utm_zone(lon)
    utm_coordinate_system.SetUTM(int(zone), int(is_northern(lat)))
    geog_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the unprojected geographic coordinate system 
    # create transform component
    geog_to_utm_transform = osr.CoordinateTransformation(geog_coordinate_system, utm_coordinate_system) # (<from>, <to>)
    easting, northing, alt = geog_to_utm_transform.TransformPoint(lon, lat, 0)
    return zone, easting, northing

def utm_to_utm(easting_in, northing_in, zone_in, zone_out):
    utm_coordinate_system_in = osr.SpatialReference()
    utm_coordinate_system_in.SetWellKnownGeogCS("WGS84") # Set unprojected geographic coordinate system
    utm_coordinate_system_in.SetUTM(int(zone_in), int(is_northern(northing_in)))
    utm_coordinate_system_out = osr.SpatialReference()
    utm_coordinate_system_out.SetWellKnownGeogCS("WGS84") # Set unprojected geographic coordinate system
    utm_coordinate_system_out.SetUTM(int(zone_out), int(is_northern(northing_in)))
    # create transform component
    utm_to_utm_transform = osr.CoordinateTransformation(utm_coordinate_system_in, utm_coordinate_system_out) # (<from>, <to>)
    easting_out, northing_out, alt = utm_to_utm_transform.TransformPoint(easting_in, northing_in, 0)
    return zone_out, easting_out, northing_out

# end UTM_Geo_Convert.py
