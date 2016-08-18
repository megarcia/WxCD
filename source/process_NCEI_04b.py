"""
Python script 'process_NCEI_04b.py'
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

USAGE: 'python process_NCEI_04b.py 1984 2013 ./analyses'

PURPOSE: Calculating statistics (mean/stdev/trend/p-value over analysis period) for aggregated 
         climatological grids. Numerous variables are addressed.

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O

INPUT: Output file from process_NCEI_04a.py in '.h5' format with aggregated grid datacubes
       (with the naming convention 'analyses/[YYYY]-[YYYY]_derived_clim_grids.h5')

OUTPUT: Same '.h5' file with new calculated statistics datacubes

RUN TIME: ~10 mins per variable
"""

import os, sys, datetime
import h5py as hdf
import numpy as np
from scipy.stats import pearsonr

def message(str):
    print str
    sys.stdout.flush()
    return    

def regress(ny,yvals):
    xvals = np.arange(0,ny)
    coeffs = np.polyfit(xvals, yvals, 1)
    corr, sig = pearsonr(xvals, yvals)
    return coeffs[0], corr, sig

def write_to_file(file,gvar,gdata):
    if gvar in file.keys():
        del file[gvar]
    file.create_dataset(gvar, data=gdata, dtype=np.float32, compression='gzip')
    message('- saving %s %s' % (gvar, str(gdata.shape)))
    return

message(' ')
message('process_NCEI_04b.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no input grids directory path indicated, using ./analyses')
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
file = '%s/%d-%d_derived_clim_grids.h5' % (path, year_begin, year_end)
message('extracting variable information from %s' % file)
h5infile = hdf.File(file,'r') 
varnames = h5infile.keys()
h5infile.close()
gridvarnames = [var for var in varnames if var[:6] == 'grids_']
ngridvars = len(gridvarnames)
message('- found %d collections of climatological grids' % ngridvars)
message(' ')
#
message('calculating statistics grids and writing to %s' % file)
h5outfile = hdf.File(file,'r+') 
if 'stats_order' not in h5outfile['meta'].keys():
    h5outfile.create_dataset('meta/stats_order', data='mean, stdev, min, max, trend, pearson r, p value')
    message('- 1 metadata item saved')
h5outfile.close()
#
for i in range(0,ngridvars):
    gridvarname = gridvarnames[i]
    statsvarname = 'stats_%s' % gridvarname[6:]
    h5infile = hdf.File(file,'r')
    gridvar = np.copy(h5infile[gridvarname])
    h5infile.close()
    message('- read %s %s' % (gridvarname, str(np.shape(gridvar))))
    nyears,nrows,ncols = np.shape(gridvar)
    statscube = np.zeros((7,nrows,ncols))
    statscube[0,:,:] = np.mean(gridvar, axis=0)
    statscube[1,:,:] = np.std(gridvar, axis=0)
    statscube[2,:,:] = np.min(gridvar, axis=0)
    statscube[3,:,:] = np.max(gridvar, axis=0)
    for j in range(0,nrows):
        for i in range(0,ncols):
            statscube[4,j,i], statscube[5,j,i], statscube[6,j,i] = regress(nyears, gridvar[:,j,i])
    h5outfile = hdf.File(file,'r+')
    write_to_file(h5outfile,statsvarname,statscube)
    h5outfile.close()
message(' ')
#
message('process_NCEI_04b.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_04b.py
