"""
Python script 'process_NCEI_13.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia

DISTRIBUTION and USE subject to 'LICENSE_GnuGPLv3.txt' and 'DISCLAIMER.txt' that accompany 
this file. This file is provided FREE OF COST and WITHOUT the author's WARRANTY or LIABILITY 
in any manner whatsoever. This file may be redistributed by the USER as long as the above 
AUTHORSHIP and COPYRIGHT, this statement, and the accompanying LICENSE and DISCLAIMER files 
remain intact. Treat others as you would be treated. Pay it forward. Valar dohaeris.

Send questions, bug reports, any related requests to matt.e.garcia@gmail.com

REFERENCE: If you use this software, please reference the following in your work products:
               Garcia, M., and P.A. Townsend, in review: "Climatological trends influencing 
               forest phenology around western Lake Superior, USA." Submitted to J. Geophys. 
               Res. Atmos. on 5 April 2016.
           See also 'README.md', 'CITATION.txt', and 'ACKNOWLEDGEMENTS.txt' for more information.

USAGE: 'python process_NCEI_13.py 1984 2013 ./analyses'

PURPOSE: Calculate climatological statistics over time on ecoregion clusters (as opposed to 
         individual ecoregions as in 'process_NCEI_10.py') 
        
DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O
              The 'Stats' and 'Plots' modules have their own requirements

INPUT: Aggregated grid datacubes and statistics grids in '.h5' file from 'process_NCEI_04.py'
       (with the naming convention 'analyses/[YYYY]-[YYYY]_derived_clim_grids.h5')
       Ecoregion maps/masks/grid info in 'clipped_ecoregions.h5' from 'process_NCEI_08.py'
       Aggregated ecoregion-based time series in '.h5' file from 'process_NCEI_10.py'
       Ecoregion clusters in '.txt' file from 'process_NCEI_12.py'

OUTPUT: Aggregated ecoregion cluster-based statistics in '.h5' files
        Maps of ecoregion clusters in '.png' format in your 'analyses' folder

RUN TIME: ~4 sec per variable per cluster

TO DO: Parallelize indicated loop using pp
"""

import os, sys, datetime
import scipy.ndimage.interpolation
import h5py as hdf
import numpy as np
from Stats import *
from Plots import *

def message(str):
    print str
    sys.stdout.flush()
    return    

def zoom(ingrids,zoomfac):
    outgrids = scipy.ndimage.interpolation.zoom(ingrids, (1,1.0/zoomfac,1.0/zoomfac), order=1)
    return outgrids

message(' ')
message('process_NCEI_13.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no input analyses directory path indicated, using ./analyses')
    path = './analyses'
else:
    path = sys.argv[3]
# 
if len(sys.argv) < 3:
    year_begin = 1984
    year_end = 2013 
else:
    year_begin = int(sys.argv[1])      
    year_end = int(sys.argv[2])      
years = np.arange(year_begin,year_end+1).astype(int)
nyears = len(years)
#
# get working area size/shape/location from clipped ecoregion maps file
infile = '%s/../data/clipped_ecoregions.h5' % path
message('extracting grid information and land mask from ecoregion maps file %s' % infile)
h5infile = hdf.File(infile,'r') 
UTM_zone = np.copy(h5infile['grid/UTM_zone'])
clip_bounds = np.copy(h5infile['grid/clip_bounds'])
UTM_bounds = clip_bounds[:4]
dx = dy = np.copy(h5infile['grid/pixelsize_orig'])
zoom_factor = np.copy(h5infile['grid/zoom_factor'])
ecoregions_map = np.copy(h5infile['eco_clip'])
landmask = np.copy(h5infile['landmask'])
nrows, ncols = np.shape(landmask)
landmask_reduced = np.copy(h5infile['landmask_reduced'])
h5infile.close()
message(' ')
#
tsfname = '%s/%d-%d_ecoregion_timeseries.h5' % (path,year_begin,year_end)
message('getting variable names, ecoregion designations, and polygon info from %s' % tsfname)
h5infile = hdf.File(tsfname,'r')
varnames = np.copy(h5infile['varnames'])
nvars = len(varnames)
ecoregion_IDs = np.copy(h5infile['ecoregion_IDs'])
necoregions = len(ecoregion_IDs)
poly_IDs = []
for i in range(0,necoregions):
    varname = 'ecoregion_%s_polygon_IDs' % ecoregion_IDs[i]
    poly_IDs.append(np.copy(h5infile[varname]))
h5infile.close()
message(' ')
#
clustersfname = '%s/%d-%d_ecoregion_clusters.txt' % (path,year_begin,year_end)
message('getting ecoregion cluster designations from %s' % clustersfname)
clusters = []
nclusters = 0
clustersf = open(clustersfname,'r')
for line in clustersf:
    line = line.rstrip()
    parts = line.split(': ')
    if ',' in parts[1]:
        ecoregions = []
        items = parts[1].split(', ')
        for item in items:
            ecoregions.append(item)
        clusters.append(ecoregions)
        nclusters += 1
    else:
        clusters.append([parts[1]])
        nclusters += 1
clustersf.close()
for i in range(0,nclusters):
    message('cluster %d contains %d ecoregions: %s' % (i+1,len(clusters[i]),str(sorted(clusters[i]))))
message(' ')
#
clusters_poly_IDs = []
for i in range(0,nclusters):
    clusters_poly_IDs.append([])
    for j in range(0,len(clusters[i])):
        for k in range(0,necoregions):
            if ecoregion_IDs[k] == clusters[i][j]:
                if len(clusters_poly_IDs[i]) == 0:
                    clusters_poly_IDs[i] = poly_IDs[k]
                else:
                    clusters_poly_IDs[i] = list(set(clusters_poly_IDs[i]) | set(poly_IDs[k]))
                break
for i in range(0,nclusters):
    message('cluster %d contains %d polygons: %s' % (i+1,len(clusters_poly_IDs[i]),str(sorted(clusters_poly_IDs[i]))))
message(' ')
#
tsfname = '%s/%d-%d_cluster_timeseries.h5' % (path,year_begin,year_end)
h5outfile = hdf.File(tsfname,'w')
h5outfile.create_dataset('meta', data='time series of all variables listed by ecoregion cluster')
h5outfile.create_dataset('ecoregion_IDs', data=ecoregion_IDs)
h5outfile.create_dataset('varnames', data=varnames)
h5outfile.create_dataset('nclusters', data=nclusters)
for i in range(0,nclusters):
    varname = 'cluster_%d_polygon_IDs' % (i+1)
    h5outfile.create_dataset(varname, data=clusters_poly_IDs[i])
h5outfile.close()
#
statsfname = '%s/%d-%d_cluster_stats.h5' % (path,year_begin,year_end)
h5outfile = hdf.File(statsfname,'w')
h5outfile.create_dataset('meta', data='statistics of all variables listed by ecoregion cluster')
h5outfile.create_dataset('ecoregion_IDs', data=ecoregion_IDs)
h5outfile.create_dataset('varnames', data=varnames)
h5outfile.create_dataset('nclusters', data=nclusters)
for i in range(0,nclusters):
    varname = 'cluster_%d_polygon_IDs' % (i+1)
    h5outfile.create_dataset(varname, data=clusters_poly_IDs[i])
h5outfile.close()
#
infile = '%s/%d-%d_derived_clim_grids.h5' % (path,year_begin,year_end)
message('getting gridded annual fields from %s' % infile)
vargrids = np.zeros((nvars-1,nyears,nrows,ncols))
h5infile = hdf.File(infile,'r') 
for i in range(1,nvars):
    var = 'grids_%s' % varnames[i]
    vargrids[i-1,:,:,:] = zoom(np.copy(h5infile[var]),zoom_factor)
    message('- %s' % varnames[i])
h5infile.close()
message(' ')
#
# TO DO: parallelize this loop at cluster-var level (getstats call)! 
# NOTE: h5py (and hdf5 libraries) typically *not* compiled for parallel h5 I/O operations 
#       -- hdf5 open command puts a lock on the h5 file until closed
#       --> save masks and stats after return from pp (could be heavy on RAM?)
nstats = 7
cluster_map = landmask
for i in range(0,nclusters):
    message('compiling individual variable statistics for cluster %d' % (i+1))
    npolys = len(clusters_poly_IDs[i])
    cluster_mask = landmask
    for j in range(0,npolys):
        cluster_mask = np.where(ecoregions_map == clusters_poly_IDs[i][j], 2, cluster_mask)
    titlestr = 'Ecoregion cluster %d' % (i+1)
    fname = '%s/cluster_maps/%d-%d_cluster_%d_mask.png' % (path,year_begin,year_end,(i+1))
    masked_map_plot_geo(cluster_mask,landmask,UTM_zone,UTM_bounds,'rainbow','none',titlestr,fname)
    cluster_mask = np.where(cluster_mask == 2, 1, 0)
    cluster_map = np.where(cluster_mask == 1, (i+1), cluster_map)
    area = np.sum(cluster_mask)*(float(dx)/1000.0)*(float(dy)/1000.0)
    cluster_series = np.zeros((nvars,nyears))
    cluster_series[0,:] = years
    cluster_stats = np.zeros((nvars,nstats))
    cluster_stats[0,0] = area
    for j in range(1,nvars):
        message('- %s' % varnames[j])
        cluster_series[j,:], cluster_stats[j,:] = getstats(vargrids[j-1,:,:,:],cluster_mask,nyears)
    #
    cluster_name = 'cluster_%d' % (i+1)
    mask_name = 'cluster_%d_mask' % (i+1)
    message('saving cluster %d mask and time series to %s' % (i+1,tsfname))
    h5outfile = hdf.File(tsfname,'r+')
    if mask_name in h5outfile.keys():
        del h5outfile[mask_name]
    h5outfile.create_dataset(mask_name, data=cluster_mask, dtype=np.int8, compression='gzip')
    if cluster_name in h5outfile.keys():
        del h5outfile[cluster_name]
    h5outfile.create_dataset(cluster_name, data=cluster_series, dtype=np.float32, compression='gzip')
    h5outfile.close()
    #
    message('saving cluster %d mask and stats to %s' % (i+1,statsfname))
    h5outfile = hdf.File(statsfname,'r+')
    if mask_name in h5outfile.keys():
        del h5outfile[mask_name]
    h5outfile.create_dataset(mask_name, data=cluster_mask, dtype=np.int8, compression='gzip')
    if cluster_name in h5outfile.keys():
        del h5outfile[cluster_name]
    h5outfile.create_dataset(cluster_name, data=cluster_stats, dtype=np.float32, compression='gzip')
    h5outfile.close()
    message(' ')
#
titlestr = 'Ecoregion cluster map'
fname = '%s/cluster_maps/%d-%d_cluster_map.png' % (path,year_begin,year_end)
masked_map_plot_geo(cluster_map,landmask,UTM_zone,UTM_bounds,'rainbow','none',titlestr,fname,0,1)
message(' ')
#
message('process_NCEI_13.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_13.py
