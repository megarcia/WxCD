"""
Python script 'process_NCEI_03_prcp_180d.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Temporal calculation of PRCP 180-day accumulation

DEPENDENCIES: h5py, numpy
              'process_NCEI_03_aux' module has its own requirements

USAGE: '$ python process_NCEI_03_prcp_180d.py NCEI_WLS_1983 1983 ./grids'

INPUT: copied '.h5' file from process_NCEI_03_preprocess.py
       (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_2.h5')

OUTPUT: updated daily '.h5' file with new accumulation grid
        (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_2.h5')
        year-end '.h5' and '.pickle' files with rolling accounted variable
"""


import sys
import datetime
import glob
import h5py as hdf
import numpy as np
from process_NCEI_03_aux import get_stn_lists, write_stn_lists, \
    write_to_file, cube_sum


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


message(' ')
message('process_NCEI_03_prcp_180d.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no input directory indicated,, using ./grids')
    path = './grids'
else:
    path = sys.argv[3]
#
if len(sys.argv) < 3:
    message('input error: need year to process')
    sys.exit(1)
else:
    this_year = int(sys.argv[2])
#
if len(sys.argv) < 2:
    message('input error: need prefix for weather data h5 file')
    sys.exit(1)
else:
    NCEIfname = sys.argv[1]
h5infname = '%s/../data/%s_processed.h5' % (path, NCEIfname)
#
message('reading dates information from %s' % h5infname)
with hdf.File(h5infname, 'r') as h5infile:
    all_dates = np.copy(h5infile['dates'])
message('- information for %d total dates found' % len(all_dates))
dates = sorted([j for j in all_dates if int(j // 1E4) == this_year])
message('- processing %d dates in %d' % (len(dates), this_year))
message(' ')
#
prev_year = this_year - 1
vars_files = sorted(glob.glob('%s/*_year_end_prcp_180d.h5' % path))
use_vars_file = False
if len(vars_files) > 0:
    for vars_file in vars_files:
        if str(prev_year) in vars_file:
            use_vars_file = True
            varfname = vars_file
            break
#
# if rolling accounting variable files exist to be carried over
# from previous year
if use_vars_file:
    message('extracting prcp_180d datacube from %s' % varfname)
    with hdf.File(varfname, 'r') as h5infile:
        nrows = np.copy(h5infile['nrows'])
        ncols = np.copy(h5infile['ncols'])
        prcp_180d = np.copy(h5infile['prcp_180d'])
    message('extracting station lists')
    prcp_180d_stns = get_stn_lists(path, prev_year, 'prcp_180d_stns')
else:  # otherwise, initialize the variable space(s)
    h5infname = '%s/%d_NCEI_grids_2.h5' % (path, dates[0])
    message('extracting grid information from %s' % h5infname)
    with hdf.File(h5infname, 'r') as h5infile:
        nrows = np.copy(h5infile['grid/nrows'])
        ncols = np.copy(h5infile['grid/ncols'])
    message('establishing prcp_180d datacube')
    prcp_180d = np.zeros((180, nrows, ncols))
    prcp_180d_stns = []
message(' ')
#
for date in dates:
    h5infname = '%s/%d_NCEI_grids_2.h5' % (path, date)
    message('extracting PRCP grid from %s' % h5infname)
    with hdf.File(h5infname, 'r') as h5infile:
        prcp_stns = np.copy(h5infile['stns/prcp_stns'])
        prcp = np.copy(h5infile['grid_prcp'])
    #
    year = date // 10000
    month = (date - (year * 10000)) // 100
    day = date - (year * 10000) - (month * 100)
    #
    grid_prcp_180d, prcp_180d_stns_all, prcp_180d, prcp_180d_stns = \
        cube_sum(180, prcp_180d, prcp, prcp_180d_stns, prcp_stns)
    message('- calculated updated 180-day running precipitation total, \
            mean %.1f' % np.mean(grid_prcp_180d))
    #
    h5outfname = '%s/%d_NCEI_grids_2.h5' % (path, date)
    message('saving grids to %s' % h5outfname)
    with hdf.File(h5outfname, 'r+') as h5outfile:
        del h5outfile['meta/last_updated']
        h5outfile.create_dataset('meta/last_updated',
                                 data=datetime.datetime.now().isoformat())
        del h5outfile['meta/at']
        outstr = 'prcp_180d'
        h5outfile.create_dataset('meta/at', data=outstr)
        write_to_file(h5outfile, 'prcp_180d_sum', grid_prcp_180d,
                      'prcp_180d_stns', prcp_180d_stns_all)
    message(' ')
#
# save rolling accounting variable for next year's run
varfname = '%s/%d_year_end_prcp_180d.h5' % (path, this_year)
message('saving variable datacube to %s' % varfname)
with hdf.File(varfname, 'w') as h5outfile:
    h5outfile.create_dataset('nrows', data=nrows)
    h5outfile.create_dataset('ncols', data=ncols)
    h5outfile.create_dataset('prcp_180d', data=prcp_180d,
                             dtype=np.float32, compression='gzip')
message('saving station lists')
write_stn_lists(path, this_year, 'prcp_180d_stns', prcp_180d_stns)
#
message('process_NCEI_03_prcp_180d.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_03_prcp_180d.py
