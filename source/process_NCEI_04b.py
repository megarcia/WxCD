"""
Python script 'process_NCEI_04b.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Calculating statistics (mean/stdev/trend/p-value over analysis period)
         for aggregated climatological grids. Numerous variables are addressed.

DEPENDENCIES: h5py, numpy, scipy.stats

USAGE: '$ python process_NCEI_04b.py 1984 2013 ./analyses'

INPUT: '.h5' output file from process_NCEI_04a.py
       (with the naming convention
        'analyses/[YYYY]-[YYYY]_derived_clim_grids.h5')

OUTPUT: Same '.h5' file with new calculated statistics datacubes
"""


import sys
import datetime
import h5py as hdf
import numpy as np
from scipy.stats import pearsonr


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


def regress(ny, yvals):
    xvals = np.arange(0, ny)
    coeffs = np.polyfit(xvals, yvals, 1)
    corr, sig = pearsonr(xvals, yvals)
    return coeffs[0], corr, sig


def write_to_file(h5file, gvar, gdata):
    if gvar in h5file.keys():
        del h5file[gvar]
    h5file.create_dataset(gvar, data=gdata, dtype=np.float32,
                          compression='gzip')
    message('- saved %s %s' % (gvar, str(gdata.shape)))
    return


message(' ')
message('process_NCEI_04b.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no input directory indicated, using ./analyses')
    path = './analyses'
else:
    path = sys.argv[3]
#
if len(sys.argv) < 3:
    message('no dates specified, analyzing 1984-2013 period')
    year_begin = 1984
    year_end = 2013
else:
    year_begin = int(sys.argv[1])
    year_end = int(sys.argv[2])
#
h5fname = '%s/%d-%d_derived_clim_grids.h5' % (path, year_begin, year_end)
message('extracting variable information from %s' % h5fname)
with hdf.File(h5fname, 'r') as h5infile:
    varnames = h5infile.keys()
gridvarnames = [var for var in varnames if var[:6] == 'grids_']
ngridvars = len(gridvarnames)
message('- found %d collections of climatological grids' % ngridvars)
message(' ')
#
message('calculating statistics grids and writing to %s' % h5fname)
with hdf.File(h5fname, 'r+') as h5outfile:
    if 'stats_order' not in h5outfile['meta'].keys():
        h5outfile.create_dataset('meta/stats_order',
                                 data='mean, stdev, min, max, trend, \
                                 pearson r, p value')
        message('- 1 metadata item saved')
#
for i, gridvarname in enumerate(gridvarnames):
    statsvarname = 'stats_%s' % gridvarname[6:]
    with hdf.File(h5fname, 'r') as h5infile:
        gridvar = np.copy(h5infile[gridvarname])
    message('- read %s %s' % (gridvarname, str(np.shape(gridvar))))
    nyears, nrows, ncols = np.shape(gridvar)
    statscube = np.zeros((7, nrows, ncols))
    statscube[0, :, :] = np.mean(gridvar, axis=0)
    statscube[1, :, :] = np.std(gridvar, axis=0)
    statscube[2, :, :] = np.min(gridvar, axis=0)
    statscube[3, :, :] = np.max(gridvar, axis=0)
    for j in range(nrows):
        for i in range(ncols):
            statscube[4, j, i], statscube[5, j, i], statscube[6, j, i] = \
                regress(nyears, gridvar[:, j, i])
    with hdf.File(h5fname, 'r+') as h5outfile:
        write_to_file(h5outfile, statsvarname, statscube)
message(' ')
#
message('process_NCEI_04b.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_04b.py
