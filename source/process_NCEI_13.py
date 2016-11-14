"""
Python script 'process_NCEI_13.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Calculate climatological statistics over time on ecoregion clusters
            (as opposed to individual ecoregions as in 'process_NCEI_10.py')

DEPENDENCIES: h5py, numpy, scipy
              The 'Stats' and 'Plots' modules have their own requirements

USAGE: '$ python process_NCEI_13.py 1984 2013 ./analyses'

INPUT: '.h5' datacubes and statistics grids in file from 'process_NCEI_04.py'
       (with the naming convention
        'analyses/[YYYY]-[YYYY]_derived_clim_grids.h5')
       'clipped_ecoregions.h5' from 'process_NCEI_08.py'
       '.h5' ecoregion-based time series from 'process_NCEI_10.py'
       '.txt' ecoregion clusters from 'process_NCEI_12.py'

OUTPUT: '.h5' ecoregion cluster-based statistics
        '.png' maps of ecoregion clusters
"""


import sys
import datetime
import scipy.ndimage.interpolation
import h5py as hdf
import numpy as np
from Stats import getstats
from Plots import masked_map_plot_geo


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


def zoom(ingrids, zoomfac):
    outgrids = \
        scipy.ndimage.interpolation.zoom(ingrids,
                                         (1, 1.0 / zoomfac, 1.0 / zoomfac),
                                         order=1)
    return outgrids


message(' ')
message('process_NCEI_13.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no directory indicated, using ./analyses')
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
#
years = np.arange(year_begin, year_end + 1).astype(int)
nyears = len(years)
#
# get working area size/shape/location from clipped ecoregion maps file
infile = '%s/../data/clipped_ecoregions.h5' % path
message('extracting grid information and land mask from %s' % infile)
with hdf.File(infile, 'r') as h5infile:
    UTM_zone = np.copy(h5infile['grid/UTM_zone'])
    clip_bounds = np.copy(h5infile['grid/clip_bounds'])
    UTM_bounds = clip_bounds[:4]
    dx = dy = np.copy(h5infile['grid/pixelsize_orig'])
    zoom_factor = np.copy(h5infile['grid/zoom_factor'])
    ecoregions_map = np.copy(h5infile['eco_clip'])
    landmask = np.copy(h5infile['landmask'])
    nrows, ncols = np.shape(landmask)
    landmask_reduced = np.copy(h5infile['landmask_reduced'])
message(' ')
#
tsfname = '%s/%d-%d_ecoregion_timeseries.h5' % (path, year_begin, year_end)
message('getting variable names, ecoregion designations, and polygon info \
    from %s' % tsfname)
poly_IDs = []
with hdf.File(tsfname, 'r') as h5infile:
    varnames = np.copy(h5infile['varnames'])
    ecoregion_IDs = np.copy(h5infile['ecoregion_IDs'])
    for ecoregion_ID in ecoregion_IDs:
        varname = 'ecoregion_%s_polygon_IDs' % ecoregion_ID
        poly_IDs.append(np.copy(h5infile[varname]))
nvars = len(varnames)
necoregions = len(ecoregion_IDs)
message(' ')
#
clustersfname = '%s/%d-%d_ecoregion_clusters.txt' % \
    (path, year_begin, year_end)
message('getting ecoregion cluster designations from %s' % clustersfname)
clusters = []
nclusters = 0
with open(clustersfname, 'r') as clustersf:
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
for i, c in enumerate(clusters):
    message('cluster %d contains %d ecoregions: %s' %
            (i + 1, len(c), str(sorted(c))))
message(' ')
#
clusters_poly_IDs = []
for i in range(nclusters):
    clusters_poly_IDs.append([])
    for j in range(len(clusters[i])):
        for k in range(necoregions):
            if ecoregion_IDs[k] == clusters[i][j]:
                if len(clusters_poly_IDs[i]) == 0:
                    clusters_poly_IDs[i] = poly_IDs[k]
                else:
                    clusters_poly_IDs[i] = \
                        list(set(clusters_poly_IDs[i]) | set(poly_IDs[k]))
                break
for i, c in enumerate(clusters_poly_IDs):
    message('cluster %d contains %d polygons: %s' %
            (i + 1, len(c), str(sorted(c))))
message(' ')
#
tsfname = '%s/%d-%d_cluster_timeseries.h5' % (path, year_begin, year_end)
with hdf.File(tsfname, 'w') as h5outfile:
    h5outfile.create_dataset('meta',
                             data='time series of all variables listed by \
                                ecoregion cluster')
    h5outfile.create_dataset('ecoregion_IDs', data=ecoregion_IDs)
    h5outfile.create_dataset('varnames', data=varnames)
    h5outfile.create_dataset('nclusters', data=nclusters)
    for i in range(nclusters):
        varname = 'cluster_%d_polygon_IDs' % (i + 1)
        h5outfile.create_dataset(varname, data=clusters_poly_IDs[i])
#
statsfname = '%s/%d-%d_cluster_stats.h5' % (path, year_begin, year_end)
with hdf.File(statsfname, 'w') as h5outfile:
    h5outfile.create_dataset('meta',
                             data='statistics of all variables listed by \
                                ecoregion cluster')
    h5outfile.create_dataset('ecoregion_IDs', data=ecoregion_IDs)
    h5outfile.create_dataset('varnames', data=varnames)
    h5outfile.create_dataset('nclusters', data=nclusters)
    for i in range(nclusters):
        varname = 'cluster_%d_polygon_IDs' % (i + 1)
        h5outfile.create_dataset(varname, data=clusters_poly_IDs[i])
#
infile = '%s/%d-%d_derived_clim_grids.h5' % (path, year_begin, year_end)
message('getting gridded annual fields from %s' % infile)
vargrids = np.zeros((nvars - 1, nyears, nrows, ncols))
with hdf.File(infile, 'r') as h5infile:
    for i in range(1, nvars):
        var = 'grids_%s' % varnames[i]
        vargrids[i - 1, :, :, :] = zoom(np.copy(h5infile[var]), zoom_factor)
        message('- %s' % varnames[i])
message(' ')
#
nstats = 7
cluster_map = landmask
for i in range(nclusters):
    message('compiling individual variable statistics for cluster %d' %
            (i + 1))
    npolys = len(clusters_poly_IDs[i])
    cluster_mask = landmask
    for j in range(npolys):
        cluster_mask = np.where(ecoregions_map == clusters_poly_IDs[i][j],
                                2, cluster_mask)
    titlestr = 'Ecoregion cluster %d' % (i + 1)
    fname = '%s/cluster_maps/%d-%d_cluster_%d_mask.png' % \
        (path, year_begin, year_end, i + 1)
    masked_map_plot_geo(cluster_mask, landmask, UTM_zone, UTM_bounds,
                        'rainbow', 'none', titlestr, fname)
    cluster_mask = np.where(cluster_mask == 2, 1, 0)
    cluster_map = np.where(cluster_mask == 1, i + 1, cluster_map)
    area = np.sum(cluster_mask) * (float(dx) / 1000.0) * (float(dy) / 1000.0)
    cluster_series = np.zeros((nvars, nyears))
    cluster_series[0, :] = years
    cluster_stats = np.zeros((nvars, nstats))
    cluster_stats[0, 0] = area
    for j in range(1, nvars):
        message('- %s' % varnames[j])
        cluster_series[j, :], cluster_stats[j, :] = \
            getstats(vargrids[j - 1, :, :, :], cluster_mask, nyears)
    #
    cluster_name = 'cluster_%d' % (i + 1)
    mask_name = 'cluster_%d_mask' % (i + 1)
    message('saving cluster %d mask and time series to %s' % (i + 1, tsfname))
    with hdf.File(tsfname, 'r+') as h5outfile:
        if mask_name in h5outfile.keys():
            del h5outfile[mask_name]
        h5outfile.create_dataset(mask_name, data=cluster_mask,
                                 dtype=np.int8, compression='gzip')
        if cluster_name in h5outfile.keys():
            del h5outfile[cluster_name]
        h5outfile.create_dataset(cluster_name, data=cluster_series,
                                 dtype=np.float32, compression='gzip')
    #
    message('saving cluster %d mask and stats to %s' % (i + 1, statsfname))
    with hdf.File(statsfname, 'r+') as h5outfile:
        if mask_name in h5outfile.keys():
            del h5outfile[mask_name]
        h5outfile.create_dataset(mask_name, data=cluster_mask,
                                 dtype=np.int8, compression='gzip')
        if cluster_name in h5outfile.keys():
            del h5outfile[cluster_name]
        h5outfile.create_dataset(cluster_name, data=cluster_stats,
                                 dtype=np.float32, compression='gzip')
    message(' ')
#
titlestr = 'Ecoregion cluster map'
fname = '%s/cluster_maps/%d-%d_cluster_map.png' % (path, year_begin, year_end)
masked_map_plot_geo(cluster_map, landmask, UTM_zone, UTM_bounds, 'rainbow',
                    'none', titlestr, fname, 0, 1)
message(' ')
#
message('process_NCEI_13.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_13.py
