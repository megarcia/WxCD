"""
Python script 'process_NCEI_10.py'
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
               Garcia, M., and P.A. Townsend (in review): "Recent climatological trends and 
               potential influences on forest phenology around western Lake Superior, USA." 
               Submitted to J. Geophys. Res. Atmos. on 5 April 2016, revised 20 August 2016.
           See also 'README.md', 'CITATION.txt', and 'ACKNOWLEDGEMENTS.txt' for more information.

USAGE: 'python process_NCEI_10.py ecoregion_polygonIDs.txt 1984 2013 ./analyses'

PURPOSE: Calculate climatological statistics over time on grid-wide and ecoregion areas

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O
              The 'Stats' and 'Plots' modules have their own requirements

INPUT: Aggregated grid datacubes and statistics grids in '.h5' file from 'process_NCEI_04.py'
       (with the naming convention 'analyses/[YYYY]-[YYYY]_derived_clim_grids.h5')
       Ecoregion maps/masks/grid info in 'clipped_ecoregions.h5' from 'process_NCEI_08.py'

OUTPUT: Aggregated full-grid and ecoregion-based statistics in '.h5' and '.csv' files
        Maps of individual ecoregions in '.png' format in your 'analyses' folder

RUN TIME: ~5 sec per ecoregion per data year

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
message('process_NCEI_10.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
# 
if len(sys.argv) < 5:
    message('input warning: no input analyses directory path indicated, using ./analyses')
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
message('getting ecoregion polygon designations from %s' % ecopolysfname)
ecoregion_IDs = []
poly_IDs = []
ecopolysf = open(ecopolysfname,'r')
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
ecopolysf.close()
necoregions = len(ecoregion_IDs)
#
varnames = []
varnames.append('area')
varnames.append('tavg_90d_avg_at_veq')
varnames.append('tavg_90d_var_at_veq')
varnames.append('prcp_90d_sum_at_veq')
varnames.append('prcp_90d_nd0_sum_at_veq')
varnames.append('prcp_90d_nd10_sum_at_veq')
varnames.append('prcp_90d_nd25_sum_at_veq')
varnames.append('tavg_90d_avg_at_ssol')
varnames.append('tavg_90d_var_at_ssol')
varnames.append('prcp_90d_sum_at_ssol')
varnames.append('prcp_90d_nd0_sum_at_ssol')
varnames.append('prcp_90d_nd10_sum_at_ssol')
varnames.append('prcp_90d_nd25_sum_at_ssol')
varnames.append('tavg_90d_avg_at_aeq')
varnames.append('tavg_90d_var_at_aeq')
varnames.append('prcp_90d_sum_at_aeq')
varnames.append('prcp_90d_nd0_sum_at_aeq')
varnames.append('prcp_90d_nd10_sum_at_aeq')
varnames.append('prcp_90d_nd25_sum_at_aeq')
varnames.append('tavg_90d_avg_at_wsol')
varnames.append('tavg_90d_var_at_wsol')
varnames.append('prcp_90d_sum_at_wsol')
varnames.append('prcp_90d_nd0_sum_at_wsol')
varnames.append('prcp_90d_nd10_sum_at_wsol')
varnames.append('prcp_90d_nd25_sum_at_wsol')
varnames.append('chill_d_at_ssol')
varnames.append('tmin_frz_days_at_ssol')
varnames.append('intensity_winter')
varnames.append('doy_last_spring_tmin_frz')
varnames.append('gdd_last_spring_tmin_frz')
varnames.append('prcp_365d_at_eoy')
varnames.append('doy_plateau_begin')
varnames.append('gdd_plateau_begin')
varnames.append('gdd_veq_to_aeq')
varnames.append('doy_first_autumn_tmin_frz')
varnames.append('doy_plateau_end')
varnames.append('gdd_plateau_end')
varnames.append('frost_free_season_days')
varnames.append('days_plateau_length')
varnames.append('gdd_plateau_length')
varnames.append('intensity_plateau')
nvars = len(varnames)
#
years = np.arange(year_begin,year_end+1).astype(int)
nyears = len(years)
#
infile = '%s/%d-%d_derived_clim_grids.h5' % (path,year_begin,year_end)
message('getting gridded annual fields from %s' % infile)
vargrids = np.zeros((nvars,nyears,nrows,ncols))
h5infile = hdf.File(infile,'r') 
for i in range(1,nvars):
    var = 'grids_%s' % varnames[i]
    vargrids[i,:,:,:] = zoom(np.copy(h5infile[var]),zoom_factor)
    message('- %s' % varnames[i])
h5infile.close()
message(' ')
#
nstats = 7
message('compiling individual variable statistics for full grid area')
area = np.sum(landmask)*(float(dx)/1000.0)*(float(dy)/1000.0)  # study area in sq km
grid_series = np.zeros((nvars,nyears))
grid_series[0,:] = years
grid_stats = np.zeros((nvars,nstats))
grid_stats[0,0] = area
for i in range(1,nvars):
    grid_series[i,:], grid_stats[i,:] = getstats(vargrids[i,:,:,:],landmask,nyears)
    message('- %s' % varnames[i])
#
tsfname = '%s/%d-%d_full_grid_timeseries.h5' % (path,year_begin,year_end)
message('saving mask and time series to %s' % tsfname)
h5outfile = hdf.File(tsfname,'w')
h5outfile.create_dataset('meta', data='time series of all variables over entire grid area')
h5outfile.create_dataset('varnames', data=varnames)
h5outfile.create_dataset('full_grid', data=grid_series, dtype=np.float32, compression='gzip')
h5outfile.close()
csvname = '%s/%d-%d_full_grid_timeseries.csv' % (path,year_begin,year_end)
np.savetxt(csvname, grid_series, delimiter=',')
#
statsfname = '%s/%d-%d_full_grid_stats.h5' % (path,year_begin,year_end)
message('saving mask and stats to %s' % statsfname)
h5outfile = hdf.File(statsfname,'w')
h5outfile.create_dataset('meta', data='statistics of all variables over entire grid area')
h5outfile.create_dataset('varnames', data=varnames)
h5outfile.create_dataset('full_grid', data=grid_stats, dtype=np.float32, compression='gzip')
h5outfile.close()
csvname = '%s/%d-%d_full_grid_stats.csv' % (path,year_begin,year_end)
np.savetxt(csvname, grid_stats, delimiter=',')
message(' ')
#
#
tsfname = '%s/%d-%d_ecoregion_timeseries.h5' % (path,year_begin,year_end)
h5outfile = hdf.File(tsfname,'w')
h5outfile.create_dataset('meta', data='variable time series listed by ecoregion')
h5outfile.create_dataset('ecoregion_IDs', data=ecoregion_IDs)
h5outfile.create_dataset('varnames', data=varnames)
for i in range(0,necoregions):
    varname = 'ecoregion_%s_polygon_IDs' % ecoregion_IDs[i]
    h5outfile.create_dataset(varname, data=poly_IDs[i])
h5outfile.close()
#
statsfname = '%s/%d-%d_ecoregion_stats.h5' % (path,year_begin,year_end)
h5outfile = hdf.File(statsfname,'w')
h5outfile.create_dataset('meta', data='variable statistics listed by ecoregion')
h5outfile.create_dataset('ecoregion_IDs', data=ecoregion_IDs)
h5outfile.create_dataset('varnames', data=varnames)
for i in range(0,necoregions):
    varname = 'ecoregion_%s_polygon_IDs' % ecoregion_IDs[i]
    h5outfile.create_dataset(varname, data=poly_IDs[i])
h5outfile.close()
#
# TO DO: parallelize this loop! 
# NOTE: h5py (and hdf5 libraries) typically *not* compiled for parallel h5 I/O operations 
#       -- hdf5 open command puts a lock on the h5 file until closed
#       --> save masks and stats after return from pp (could be heavy on RAM?)
for i in range(0,necoregions):
    message('compiling individual variable statistics for ecoregion %s' % ecoregion_IDs[i])
    npolys = len(poly_IDs[i])
    ecoregion_mask = landmask
    for j in range(0,npolys):
        ecoregion_mask = np.where(ecoregions_map == poly_IDs[i][j], 2, ecoregion_mask)
    titlestr = 'Ecoregion %s' % ecoregion_IDs[i]
    fname = '%s/ecoregion_maps/ecoregion_%s_mask.png' % (path,ecoregion_IDs[i])
    masked_map_plot_geo(ecoregion_mask,landmask,UTM_zone,UTM_bounds,'rainbow','none',titlestr,fname)
    ecoregion_mask = np.where(ecoregion_mask == 2, 1, 0)
    area = np.sum(ecoregion_mask)*(float(dx)/1000.0)*(float(dy)/1000.0)  # ecoregion area in sq km
    ecoregion_series = np.zeros((nvars, nyears))
    ecoregion_series[0,:] = years
    ecoregion_stats = np.zeros((nvars, nstats))
    ecoregion_stats[0,0] = area
    for j in range(1,nvars):
        message('- %s' % varnames[j])
        ecoregion_series[j,:], ecoregion_stats[j,:] = getstats(vargrids[j,:,:,:],ecoregion_mask,nyears)
    mask_name = ecoregion_IDs[i] + '_mask'
    message('saving mask and time series to %s' % tsfname)
    h5outfile = hdf.File(tsfname,'r+')
    h5outfile.create_dataset(mask_name, data=ecoregion_mask, dtype=np.int8, compression='gzip')
    h5outfile.create_dataset(ecoregion_IDs[i], data=ecoregion_series, dtype=np.float32, compression='gzip')
    h5outfile.close()
    message('saving mask and stats to %s' % statsfname)
    h5outfile = hdf.File(statsfname,'r+')
    h5outfile.create_dataset(mask_name, data=ecoregion_mask, dtype=np.int8, compression='gzip')
    h5outfile.create_dataset(ecoregion_IDs[i], data=ecoregion_stats, dtype=np.float32, compression='gzip')
    h5outfile.close()
    message(' ')
#
message('process_NCEI_10.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_10.py
