"""
Python module 'UTM_Geo_Convert.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2014-2016 by Matthew Garcia
heavily modified from python source code found at
http://stackoverflow.com/questions/343865/how-to-convert-from-utm-to-latlng-in-python-or-javascript/10239676
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Conversion between UTM and geographic coordinates, or from one UTM
         zone to another

DEPENDENCIES: osgeo.osr (uses gdal)

USAGE: insert 'from UTM_Geo_Convert import *' near head of script, then
       (for example)
        UTM_coords = geographic_to_utm(stn_lon, stn_lat, UTMzone)
        stn_easting = int(round(UTM_coords[1],0))
        stn_northing = int(round(UTM_coords[2],0))

INPUT: coordinates provided by calling script

OUTPUT: coordinates returned to calling script
"""


from osgeo import osr  # uses GDAL via osr-to-python bindings


def get_utm_zone(lon):
    """ get proper UTM zone based on longitude """
    zone = int(1 + (lon + 180.0) / 6.0)
    return zone


def is_northern(lat):
    """ check if latitude is in Northern Hemisphere (thus + northing) """
    if lat < 0.0:
        return 0
    else:
        return 1


def utm_to_geographic(easting, northing, zone):
    """
    convert from UTM coordinates (zone, easting, northing)
            to geographic coordinates (longitude, latitude)
    """
    utm_coordinate_system = osr.SpatialReference()
    # Set unprojected geographic coordinate system
    utm_coordinate_system.SetWellKnownGeogCS("WGS84")
    utm_coordinate_system.SetUTM(int(zone), int(is_northern(northing)))
    # Clone ONLY the unprojected geographic coordinate system
    geog_coordinate_system = utm_coordinate_system.CloneGeogCS()
    # Create transform component with (<from>, <to>)
    utm_to_geog_transform = \
        osr.CoordinateTransformation(utm_coordinate_system,
                                     geog_coordinate_system)
    # Note returned 'alt' (altitude) is currently unused, thus '_'
    lon, lat, _ = utm_to_geog_transform.TransformPoint(easting, northing, 0)
    return lon, lat


def geographic_to_utm(lon, lat, zonepref=-1):
    """
    convert from geographic coordinates (longitude, latitude)
            to UTM coordinates (zone, easting, northing)
    """
    utm_coordinate_system = osr.SpatialReference()
    # Set unprojected geographic coordinate system
    utm_coordinate_system.SetWellKnownGeogCS("WGS84")
    if zonepref != -1:
        zone = zonepref
    else:
        zone = get_utm_zone(lon)
    utm_coordinate_system.SetUTM(int(zone), int(is_northern(lat)))
    # Clone ONLY the unprojected geographic coordinate system
    geog_coordinate_system = utm_coordinate_system.CloneGeogCS()
    # Create transform component with (<from>, <to>)
    geog_to_utm_transform = \
        osr.CoordinateTransformation(geog_coordinate_system,
                                     utm_coordinate_system)
    # Note returned 'alt' (altitude) is currently unused, thus '_'
    easting, northing, _ = geog_to_utm_transform.TransformPoint(lon, lat, 0)
    return zone, easting, northing


def utm_to_utm(easting_in, northing_in, zone_in, zone_out):
    """
    convert from UTM coordinates (zone_in, easting, northing)
            to UTM coordinates (zone_out, easting, northing)
    generally only useful for coordinate translation near edges
        of neighboring zones
    """
    utm_coordinate_system_in = osr.SpatialReference()
    # Set unprojected geographic coordinate system
    utm_coordinate_system_in.SetWellKnownGeogCS("WGS84")
    utm_coordinate_system_in.SetUTM(int(zone_in),
                                    int(is_northern(northing_in)))
    utm_coordinate_system_out = osr.SpatialReference()
    # Set unprojected geographic coordinate system
    utm_coordinate_system_out.SetWellKnownGeogCS("WGS84")
    utm_coordinate_system_out.SetUTM(int(zone_out),
                                     int(is_northern(northing_in)))
    # Create transform component with (<from>, <to>)
    utm_to_utm_transform = \
        osr.CoordinateTransformation(utm_coordinate_system_in,
                                     utm_coordinate_system_out)
    # Note returned 'alt' (altitude) is currently unused, thus '_'
    easting_out, northing_out, _ = \
        utm_to_utm_transform.TransformPoint(easting_in, northing_in, 0)
    return zone_out, easting_out, northing_out

# end UTM_Geo_Convert.py
