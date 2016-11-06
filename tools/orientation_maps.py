"""
Python script "orientation_maps.py"
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Create projected maps of study area and meteorological stations

DEPENDENCIES: h5py, numpy, pandas, matplotlib, mpl_toolkits
              'UTM_Geo_Convert' module (which has its own requirements)

USAGE: '$ python orientation_maps.py NCEI_WLS_19830101-20151031_stnmeta.csv
          clipped_ecoregions.h5'

INPUT: <NCEI_WLS_19830101-20151031_stnmeta.csv> is an output of
       process_NCEI_01.py
       <clipped_ecoregions.h5> is an output of process_NCEI_08.py

NOTE: run this script *after* you've run process_NCEI_01.py and _08.py so that
      the proper input files are in place

OUTPUT: 2 '.png' map figures
"""


import sys
import datetime
import h5py as hdf
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from UTM_Geo_Convert import utm_to_geographic


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


message(' ')
message('orientation_maps.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 3:
    message('input error: need h5 file containing study area boundary \
            information')
    sys.exit(1)
else:
    maskfname = sys.argv[2]
#
if len(sys.argv) < 2:
    message('input error: need csv file containing processed station metadata \
            information')
    sys.exit(1)
else:
    stnfname = sys.argv[1]
#
message('reading input data file %s' % stnfname)
stndata_df = pd.read_csv(stnfname)
nstns, nstncols = np.shape(stndata_df)
message('- found %d total stations with location and data information' % nstns)
lats_vals = list(stndata_df['LATITUDE'])
lons_vals = list(stndata_df['LONGITUDE'])
ndates_vals = list(np.array(stndata_df['NDATES']) / 365.25)
nprcp_vals = list(stndata_df['NPRCP'])
ntmax_vals = list(stndata_df['NTMAX'])
ntmin_vals = list(stndata_df['NTMIN'])
#
# get bounding coordinates for station collection
min_lat = np.floor(np.min(lats_vals))
max_lat = np.ceil(np.max(lats_vals))
min_lon = np.floor(np.min(lons_vals))
max_lon = np.ceil(np.max(lons_vals))
mid_lat = (min_lat + max_lat) / 2.0
mid_lon = (min_lon + max_lon) / 2.0
#
# create polar stereographic reference Basemap
message('plotting reference orientation map')
SWcorner = np.array([25.0, -110.0], dtype=np.float32)
NEcorner = np.array([55.0, -50.0], dtype=np.float32)
r = Basemap(projection='stere', lon_0=mid_lon, lat_0=90.0, lat_ts=45.0,
            llcrnrlat=SWcorner[0], llcrnrlon=SWcorner[1],
            urcrnrlat=NEcorner[0], urcrnrlon=NEcorner[1],
            resolution='i', area_thresh=10000)
fig1 = plt.figure(figsize=(4, 4))
ax = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
r.drawcoastlines()
r.drawcountries()
r.drawstates()
parallels = np.arange(0., 90., 10.)
r.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10)
meridians = np.arange(180., 360., 10.)
r.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=10)
x = [min_lon, min_lon, max_lon, max_lon, min_lon]
y = [min_lat, max_lat, max_lat, min_lat, min_lat]
r.plot(x, y, latlon=True, linewidth=2, color='r')
fname = 'map_orientation.png'
message('saving figure as %s' % fname)
plt.savefig(fname, dpi=300, bbox_inches='tight')
plt.clf()
#
# get actual study area boundaries
message('extracting study area boundaries from %s' % maskfname)
with hdf.File(maskfname, 'r') as h5infile:
    UTM_zone = np.copy(h5infile['grid/UTM_zone'])
    clip_bounds = np.copy(h5infile['grid/clip_bounds'])
min_east = clip_bounds[0]
max_north = clip_bounds[1]
max_east = clip_bounds[2]
min_north = clip_bounds[3]
sa_min_lon, sa_min_lat = utm_to_geographic(min_east, min_north, UTM_zone)
sa_max_lon, sa_max_lat = utm_to_geographic(max_east, max_north, UTM_zone)
sa_max_lon = sa_max_lon - 0.25
#
# sort stations by data provided (P only, T only, P + T)
lonsP = []
latsP = []
nyP = []
lonsT = []
latsT = []
nyT = []
lonsPT = []
latsPT = []
nyPT = []
for i in range(nstns):
    if (nprcp_vals[i] != 0) and (ntmax_vals[i] != 0):
        lonsPT.append(lons_vals[i])
        latsPT.append(lats_vals[i])
        nyPT.append(ndates_vals[i])
    if (nprcp_vals[i] == 0) and (ntmax_vals[i] != 0):
        lonsT.append(lons_vals[i])
        latsT.append(lats_vals[i])
        nyT.append(ndates_vals[i])
    if (nprcp_vals[i] != 0) and (ntmax_vals[i] == 0):
        lonsP.append(lons_vals[i])
        latsP.append(lats_vals[i])
        nyP.append(ndates_vals[i])
#
# create UTM station coverage Basemap
message('plotting station coverage map')
message('- %d P + T stations' % len(nyPT))
message('- %d P only stations' % len(nyP))
message('- %d T only stations' % len(nyT))
SWcorner = np.array([np.min([min_lat, sa_min_lat]),
                     np.min([min_lon, sa_min_lon])], dtype=np.float32)
NEcorner = np.array([np.max([max_lat, sa_max_lat]),
                     np.max([max_lon, sa_max_lon]) + 0.75], dtype=np.float32)
d1 = Basemap(projection='tmerc', lon_0=mid_lon, lat_0=90.0, lat_ts=mid_lat,
             llcrnrlat=SWcorner[0], llcrnrlon=SWcorner[1],
             urcrnrlat=NEcorner[0], urcrnrlon=NEcorner[1],
             resolution='h', area_thresh=500)
fig2 = plt.figure(figsize=(8, 8))
ax = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
d1.drawcoastlines()
d1.drawcountries()
d1.drawstates()
d1.drawmapscale(lon=-95.2, lat=44.4, lon0=mid_lon, lat0=mid_lat,
                length=100.0, barstyle='fancy')
parallels = np.arange(0., 90., 2.)
d1.drawparallels(parallels, labels=[1, 1, 0, 0], fontsize=10)
meridians = np.arange(180., 360., 2.)
d1.drawmeridians(meridians, labels=[0, 0, 1, 1], fontsize=10)
sa_x = [sa_min_lon, sa_min_lon, sa_max_lon, sa_max_lon, sa_min_lon]
sa_y = [sa_min_lat, sa_max_lat, sa_max_lat, sa_min_lat, sa_min_lat]
d1.plot(sa_x, sa_y, latlon=True, linewidth=3, linestyle='--', color='k')
cmap = plt.cm.YlGnBu
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.linspace(0, 32, 9)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
dPT_label = 'P + T (%d)' % len(nyPT)
dPT = d1.scatter(lonsPT, latsPT, 40, latlon=True, marker='o', c=nyPT,
                 cmap=cmap, label=dPT_label)
dP_label = 'P only (%d)' % len(nyP)
dP = d1.scatter(lonsP, latsP, 30, latlon=True, marker='s', c=nyP,
                cmap=cmap, label=dP_label)
dT_label = 'T only (%d)' % len(nyT)
dT = d1.scatter(lonsT, latsT, 30, latlon=True, marker='^', c=nyT,
                cmap=cmap, label=dT_label)
cbar = d1.colorbar(dPT, location='bottom', pad="5%", cmap=cmap, norm=norm,
                   boundaries=bounds, format='%1i')
cbar.set_label('data years')
legend = plt.legend(loc='upper left', fontsize='12')
legend.get_frame().set_facecolor('#EEEEEE')
fname = 'map_stations.png'
message('saving figure as %s' % fname)
plt.savefig(fname, dpi=300, bbox_inches='tight')
plt.clf()
#
message(' ')
message('orientation_maps.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end orientation_maps.py
