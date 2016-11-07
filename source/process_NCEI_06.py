"""
Python script 'process_NCEI_06.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Summarizing key meteorological variables on a daily basis

DEPENDENCIES: h5py, numpy
              'Date_Convert' module has no external requirements

USAGE: '$ python process_NCEI_06.py 1984 2013 ./grids'

INPUT: '.h5' output from process_NCEI_03.py
       (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_2.h5')

OUTPUT: New '.csv' files with daily time series (7 variables)
        (not sure if we need these any longer)
        '.h5' file with daily time series and aggregated mean daily values
        (with the naming convention
         'analyses/[YYYY]-[YYYY]_clim_values_by_doy.h5')
"""


import sys
import datetime
import glob
import h5py as hdf
import numpy as np
from Date_Convert import doy_to_date


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


def write_to_file(h5file, gvar, gdata):
    if gvar in h5file.keys():
        del h5file[gvar]
    h5file.create_dataset(gvar, data=gdata, dtype=np.float32,
                          compression='gzip')
    message('- %s %s' % (gvar, str(gdata.shape)))
    return


message(' ')
message('process_NCEI_06.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no input grids directory indicated, using ./grids')
    path = './grids'
else:
    path = sys.argv[3]
#
if len(sys.argv) < 3:
    year_begin = 1984
    year_end = 2013
else:
    year_begin = int(sys.argv[1])
    year_end = int(sys.argv[2])
years = np.arange(year_begin, year_end + 1).astype(int)
nyears = len(years)
#
doy_begin = 1
doy_end = 365
doys = np.arange(doy_begin, doy_end + 1).astype(int)
ndoys = len(doys)
#
wxlist = glob.glob('%s/*_NCEI_grids_2.h5' % path)
message('found %d weather derivative grid files' % len(wxlist))
message(' ')
#
infile = wxlist[0]
message('extracting grid info from weather derivatives file %s' % infile)
with hdf.File(infile, 'r') as h5infile:
    ncols = np.copy(h5infile['grid/ncols'])
    nrows = np.copy(h5infile['grid/nrows'])
message(' ')
#
prcp_by_year_doy = np.zeros((nyears, ndoys))
tmin_by_year_doy = np.zeros((nyears, ndoys))
tmax_by_year_doy = np.zeros((nyears, ndoys))
tavg_by_year_doy = np.zeros((nyears, ndoys))
cd_by_year_doy = np.zeros((nyears, ndoys))
cdd_by_year_doy = np.zeros((nyears, ndoys))
gdd_by_year_doy = np.zeros((nyears, ndoys))
#
for j, year in enumerate(years):
    message('processing grids for %d' % year)
    for i, doy in enumerate(doys):
        mmdd = doy_to_date(year, doy)
        filename = '%s/%d%d_NCEI_grids_2.h5' % (path, year, str(mmdd).zfill(4))
        with hdf.File(filename, 'r') as h5infile:
            grid_prcp = np.copy(h5infile['grid_prcp'])
            prcp_by_year_doy[j, i] = np.mean(grid_prcp)
            grid_tmin = np.copy(h5infile['grid_tmin'])
            tmin_by_year_doy[j, i] = np.mean(grid_tmin)
            grid_tmax = np.copy(h5infile['grid_tmax'])
            tmax_by_year_doy[j, i] = np.mean(grid_tmax)
            grid_tavg = np.copy(h5infile['grid_tavg'])
            tavg_by_year_doy[j, i] = np.mean(grid_tavg)
            grid_cd = np.copy(h5infile['chill_d'])
            cd_by_year_doy[j, i] = np.mean(grid_cd)
            grid_cdd = np.copy(h5infile['chill_dd'])
            cdd_by_year_doy[j, i] = np.mean(grid_cdd)
            grid_gdd = np.copy(h5infile['grow_dd'])
            gdd_by_year_doy[j, i] = np.mean(grid_gdd)
message(' ')
#
outfile = '%s/../analyses/%d-%d_clim_values_by_doy.h5' % \
    (path, year_begin, year_end)
message('writing metadata and series information to %s' % outfile)
with hdf.File(outfile, 'w') as h5outfile:
    h5outfile.create_dataset('meta/filename', data=outfile)
    h5outfile.create_dataset('meta/created',
                             data=datetime.datetime.now().isoformat())
    h5outfile.create_dataset('meta/by',
                             data='M. Garcia, UWisconsin-Madison FWE')
    h5outfile.create_dataset('meta/last_updated',
                             data=datetime.datetime.now().isoformat())
    h5outfile.create_dataset('meta/at',
                             data='climatological and statistics values by \
                                DOY')
    message('- 5 metadata items saved')
    write_to_file(h5outfile, 'prcp_by_year_doy', prcp_by_year_doy)
    write_to_file(h5outfile, 'tmin_by_year_doy', tmin_by_year_doy)
    write_to_file(h5outfile, 'tmax_by_year_doy', tmax_by_year_doy)
    write_to_file(h5outfile, 'tavg_by_year_doy', tavg_by_year_doy)
    write_to_file(h5outfile, 'cd_by_year_doy', cd_by_year_doy)
    write_to_file(h5outfile, 'cdd_by_year_doy', cd_by_year_doy)
    write_to_file(h5outfile, 'gdd_by_year_doy', gdd_by_year_doy)
    message('- 6 series arrays saved')
message(' ')
#
prcp_mean_by_doy = np.mean(prcp_by_year_doy, axis=0)
prcp_min_by_doy = np.min(prcp_by_year_doy, axis=0)
prcp_max_by_doy = np.max(prcp_by_year_doy, axis=0)
#
tmin_mean_by_doy = np.mean(tmin_by_year_doy, axis=0)
tmin_min_by_doy = np.min(tmin_by_year_doy, axis=0)
tmin_max_by_doy = np.max(tmin_by_year_doy, axis=0)
#
tmax_mean_by_doy = np.mean(tmax_by_year_doy, axis=0)
tmax_min_by_doy = np.min(tmax_by_year_doy, axis=0)
tmax_max_by_doy = np.max(tmax_by_year_doy, axis=0)
#
tavg_mean_by_doy = np.mean(tavg_by_year_doy, axis=0)
tavg_min_by_doy = np.min(tavg_by_year_doy, axis=0)
tavg_max_by_doy = np.max(tavg_by_year_doy, axis=0)
#
cd_mean_by_doy = np.mean(cd_by_year_doy, axis=0)
cd_min_by_doy = np.min(cd_by_year_doy, axis=0)
cd_max_by_doy = np.max(cd_by_year_doy, axis=0)
#
cdd_mean_by_doy = np.mean(cdd_by_year_doy, axis=0)
cdd_min_by_doy = np.min(cdd_by_year_doy, axis=0)
cdd_max_by_doy = np.max(cdd_by_year_doy, axis=0)
#
gdd_mean_by_doy = np.mean(gdd_by_year_doy, axis=0)
gdd_min_by_doy = np.min(gdd_by_year_doy, axis=0)
gdd_max_by_doy = np.max(gdd_by_year_doy, axis=0)
#
message('writing summary values by DOY to %s' % outfile)
with hdf.File(outfile, 'r+') as h5outfile:
    write_to_file(h5outfile, 'prcp_mean_by_doy', prcp_mean_by_doy)
    write_to_file(h5outfile, 'prcp_min_by_doy', prcp_min_by_doy)
    write_to_file(h5outfile, 'prcp_max_by_doy', prcp_max_by_doy)
    write_to_file(h5outfile, 'tmin_mean_by_doy', tmin_mean_by_doy)
    write_to_file(h5outfile, 'tmin_min_by_doy', tmin_min_by_doy)
    write_to_file(h5outfile, 'tmin_max_by_doy', tmin_max_by_doy)
    write_to_file(h5outfile, 'tmax_mean_by_doy', tmax_mean_by_doy)
    write_to_file(h5outfile, 'tmax_min_by_doy', tmax_min_by_doy)
    write_to_file(h5outfile, 'tmax_max_by_doy', tmax_max_by_doy)
    write_to_file(h5outfile, 'tavg_mean_by_doy', tavg_mean_by_doy)
    write_to_file(h5outfile, 'tavg_min_by_doy', tavg_min_by_doy)
    write_to_file(h5outfile, 'tavg_max_by_doy', tavg_max_by_doy)
    write_to_file(h5outfile, 'cd_mean_by_doy', cd_mean_by_doy)
    write_to_file(h5outfile, 'cd_min_by_doy', cd_min_by_doy)
    write_to_file(h5outfile, 'cd_max_by_doy', cd_max_by_doy)
    write_to_file(h5outfile, 'cdd_mean_by_doy', cdd_mean_by_doy)
    write_to_file(h5outfile, 'cdd_min_by_doy', cdd_min_by_doy)
    write_to_file(h5outfile, 'cdd_max_by_doy', cdd_max_by_doy)
    write_to_file(h5outfile, 'gdd_mean_by_doy', gdd_mean_by_doy)
    write_to_file(h5outfile, 'gdd_min_by_doy', gdd_min_by_doy)
    write_to_file(h5outfile, 'gdd_max_by_doy', gdd_max_by_doy)
    message('- 21 series saved')
message(' ')
#
message('process_NCEI_06.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_06.py
