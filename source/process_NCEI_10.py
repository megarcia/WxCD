"""
Python script 'process_NCEI_10.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Calculate climatological statistics over time on grid-wide and
         ecoregion areas

DEPENDENCIES: h5py, numpy, scipy
              'Stats' and 'Plots' modules have their own requirements

USAGE: '$ python process_NCEI_10.py ecoregion_polygonIDs.txt 1984 2013
          ./analyses'

INPUT: '.h5' datacubes and statistics grids from 'process_NCEI_04b.py'
       (with the naming convention
        'analyses/[YYYY]-[YYYY]_derived_clim_grids.h5')
       'clipped_ecoregions.h5' from 'process_NCEI_08.py'

OUTPUT: '.h5' and '.csv' aggregated full-grid and ecoregion-based statistics
        '.png' ecoregion maps
"""


import sys
import datetime
import h5py as hdf
import numpy as np
import scipy.ndimage.interpolation
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
message('process_NCEI_10.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 5:
    message('input warning: no directory path indicated, using ./analyses')
    path = './analyses'
else:
    path = sys.argv[4]
#
if len(sys.argv) < 4:
    year_begin = 1984
    year_end = 2013
else:
    year_begin = int(sys.argv[2])
    year_end = int(sys.argv[3])
#
if len(sys.argv) < 2:
    message('input error: need polygon index/ID txt file')
    sys.exit(1)
else:
    ecopolysfname = '%s/../data/%s' % (path, sys.argv[1])
#
# get working area size/shape/location from clipped ecoregion maps file
infile = '%s/../data/clipped_ecoregions.h5' % path
message('extracting grid info and land mask from %s' % infile)
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
message('getting ecoregion polygon designations from %s' % ecopolysfname)
ecoregion_IDs = []
poly_IDs = []
with open(ecopolysfname, 'r') as ecopolysf:
    for line in ecopolysf:
        line = line.rstrip()
        parts = line.split(': ')
        ecoregion_IDs.append(str(parts[0]))
        if ',' in parts[1]:
            polys = []
            items = parts[1].split(',')
            for item in items:
                polys.append(int(item))
            poly_IDs.append(polys)
        else:
            poly_IDs.append([int(parts[1])])
necoregions = len(ecoregion_IDs)
#
varnames = ['area', 'tavg_90d_avg_at_veq', 'tavg_90d_var_at_veq',
            'prcp_90d_sum_at_veq', 'prcp_90d_nd0_sum_at_veq',
            'prcp_90d_nd10_sum_at_veq', 'prcp_90d_nd25_sum_at_veq',
            'tavg_90d_avg_at_ssol', 'tavg_90d_var_at_ssol',
            'prcp_90d_sum_at_ssol', 'prcp_90d_nd0_sum_at_ssol',
            'prcp_90d_nd10_sum_at_ssol', 'prcp_90d_nd25_sum_at_ssol',
            'tavg_90d_avg_at_aeq', 'tavg_90d_var_at_aeq',
            'prcp_90d_sum_at_aeq', 'prcp_90d_nd0_sum_at_aeq',
            'prcp_90d_nd10_sum_at_aeq', 'prcp_90d_nd25_sum_at_aeq',
            'tavg_90d_avg_at_wsol', 'tavg_90d_var_at_wsol',
            'prcp_90d_sum_at_wsol', 'prcp_90d_nd0_sum_at_wsol',
            'prcp_90d_nd10_sum_at_wsol', 'prcp_90d_nd25_sum_at_wsol',
            'chill_d_at_ssol', 'tmin_frz_days_at_ssol', 'intensity_winter',
            'doy_last_spring_tmin_frz', 'gdd_last_spring_tmin_frz',
            'prcp_365d_at_eoy', 'doy_plateau_begin', 'gdd_plateau_begin',
            'gdd_veq_to_aeq', 'doy_first_autumn_tmin_frz', 'doy_plateau_end',
            'gdd_plateau_end', 'frost_free_season_days', 'days_plateau_length',
            'gdd_plateau_length', 'intensity_plateau']
nvars = len(varnames)
#
years = np.arange(year_begin, year_end + 1).astype(int)
nyears = len(years)
#
infile = '%s/%d-%d_derived_clim_grids.h5' % (path, year_begin, year_end)
message('getting gridded annual fields from %s' % infile)
vargrids = np.zeros((nvars, nyears, nrows, ncols))
with hdf.File(infile, 'r') as h5infile:
    for i in range(1, nvars):
        var = 'grids_%s' % varnames[i]
        vargrids[i, :, :, :] = zoom(np.copy(h5infile[var]), zoom_factor)
        message('- %s' % varnames[i])
message(' ')
#
nstats = 7
message('compiling individual variable statistics for full grid area')
# study area in sq km
area = np.sum(landmask) * (float(dx) / 1000.0) * (float(dy) / 1000.0)
grid_series = np.zeros((nvars, nyears))
grid_series[0, :] = years
grid_stats = np.zeros((nvars, nstats))
grid_stats[0, 0] = area
for i in range(1, nvars):
    grid_series[i, :], grid_stats[i, :] = \
        getstats(vargrids[i, :, :, :], landmask, nyears)
    message('- %s' % varnames[i])
#
tsfname = '%s/%d-%d_full_grid_timeseries.h5' % (path, year_begin, year_end)
message('saving mask and time series to %s' % tsfname)
with hdf.File(tsfname, 'w') as h5outfile:
    h5outfile.create_dataset('meta',
                             data='time series of all variables over entire \
                                  grid area')
    h5outfile.create_dataset('varnames', data=varnames)
    h5outfile.create_dataset('full_grid', data=grid_series,
                             dtype=np.float32, compression='gzip')
csvname = '%s/%d-%d_full_grid_timeseries.csv' % (path, year_begin, year_end)
np.savetxt(csvname, grid_series, delimiter=',')
#
statsfname = '%s/%d-%d_full_grid_stats.h5' % (path, year_begin, year_end)
message('saving mask and stats to %s' % statsfname)
with hdf.File(statsfname, 'w') as h5outfile:
    h5outfile.create_dataset('meta',
                             data='statistics of all variables over entire \
                                  grid area')
    h5outfile.create_dataset('varnames', data=varnames)
    h5outfile.create_dataset('full_grid', data=grid_stats,
                             dtype=np.float32, compression='gzip')
csvname = '%s/%d-%d_full_grid_stats.csv' % (path, year_begin, year_end)
np.savetxt(csvname, grid_stats, delimiter=',')
message(' ')
#
#
tsfname = '%s/%d-%d_ecoregion_timeseries.h5' % (path, year_begin, year_end)
with hdf.File(tsfname, 'w') as h5outfile:
    h5outfile.create_dataset('meta',
                             data='variable time series listed by ecoregion')
    h5outfile.create_dataset('ecoregion_IDs', data=ecoregion_IDs)
    h5outfile.create_dataset('varnames', data=varnames)
    for i in range(necoregions):
        varname = 'ecoregion_%s_polygon_IDs' % ecoregion_IDs[i]
        h5outfile.create_dataset(varname, data=poly_IDs[i])
#
statsfname = '%s/%d-%d_ecoregion_stats.h5' % (path, year_begin, year_end)
with hdf.File(statsfname, 'w') as h5outfile:
    h5outfile.create_dataset('meta',
                             data='variable statistics listed by ecoregion')
    h5outfile.create_dataset('ecoregion_IDs', data=ecoregion_IDs)
    h5outfile.create_dataset('varnames', data=varnames)
    for i in range(necoregions):
        varname = 'ecoregion_%s_polygon_IDs' % ecoregion_IDs[i]
        h5outfile.create_dataset(varname, data=poly_IDs[i])
#
for i in range(necoregions):
    message('compiling individual variable statistics for ecoregion %s' %
            ecoregion_IDs[i])
    npolys = len(poly_IDs[i])
    ecoregion_mask = landmask
    for j in range(npolys):
        ecoregion_mask = \
            np.where(ecoregions_map == poly_IDs[i][j], 2, ecoregion_mask)
    titlestr = 'Ecoregion %s' % ecoregion_IDs[i]
    fname = '%s/ecoregion_maps/ecoregion_%s_mask.png' % \
        (path, ecoregion_IDs[i])
    masked_map_plot_geo(ecoregion_mask, landmask, UTM_zone, UTM_bounds,
                        'rainbow', 'none', titlestr, fname)
    ecoregion_mask = np.where(ecoregion_mask == 2, 1, 0)
    # ecoregion area in sq km
    area = np.sum(ecoregion_mask) * (float(dx) / 1000.0) * (float(dy) / 1000.0)
    ecoregion_series = np.zeros((nvars, nyears))
    ecoregion_series[0, :] = years
    ecoregion_stats = np.zeros((nvars, nstats))
    ecoregion_stats[0, 0] = area
    for j in range(1, nvars):
        message('- %s' % varnames[j])
        ecoregion_series[j, :], ecoregion_stats[j, :] = \
            getstats(vargrids[j, :, :, :], ecoregion_mask, nyears)
    mask_name = ecoregion_IDs[i] + '_mask'
    message('saving mask and time series to %s' % tsfname)
    with hdf.File(tsfname, 'r+') as h5outfile:
        h5outfile.create_dataset(mask_name, data=ecoregion_mask,
                                 dtype=np.int8, compression='gzip')
        h5outfile.create_dataset(ecoregion_IDs[i], data=ecoregion_series,
                                 dtype=np.float32, compression='gzip')
    message('saving mask and stats to %s' % statsfname)
    with hdf.File(statsfname, 'r+') as h5outfile:
        h5outfile.create_dataset(mask_name, data=ecoregion_mask,
                                 dtype=np.int8, compression='gzip')
        h5outfile.create_dataset(ecoregion_IDs[i], data=ecoregion_stats,
                                 dtype=np.float32, compression='gzip')
    message(' ')
#
message('process_NCEI_10.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_10.py
