"""
Python script 'process_NCEI_03_preprocess.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Copy *_grids_1.h5 files to *_grids_2.h5

DEPENDENCIES: h5py, numpy

USAGE: '$ python process_NCEI_03_preprocess_sc.py NCEI_WLS_1983 1983 ./grids'

INPUT: '.h5' output files from process_NCEI_02.py
       (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_1.h5')

OUTPUT: Copied '.h5' files (1 / day)
        (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_2.h5')
"""


import os
import sys
import datetime
import glob
import h5py as hdf
import numpy as np


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


message(' ')
message('process_NCEI_03_preprocess.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no input directory indicated, using ./grids')
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
    message('input error: need prefix for weather data .h5 file')
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
message('copying grid files')
for date in dates:
    h5infname = '%s/%d_NCEI_grids_1.h5' % (path, date)
    message('%s' % h5infname)
    h5outfname = '%s/%d_NCEI_grids_2.h5' % (path, date)
    cmdstring = 'cp %s %s' % (h5infname, h5outfname)
    os.system(cmdstring)
    message('--> %s' % h5outfname)
#
message(' ')
message('process_NCEI_03_preprocess.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_03_preprocess.py
