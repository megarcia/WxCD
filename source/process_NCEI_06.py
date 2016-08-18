"""
Python script 'process_NCEI_06.py'
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

USAGE: 'python process_NCEI_06.py 1984 2013 ./grids'

PURPOSE: Summarizing key meteorological variables on a daily basis

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O
              The 'Date_Convert' module has no external requirements

INPUT: Output files from process_NCEI_03.py script in '.h5' format (1 / day)
       (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_2.h5')

OUTPUT: New '.csv' files with daily time series (7 variables) <-- not sure if we need these any longer
        New '.h5' file with daily time series and aggregated mean daily values
        (with the naming convention 'analyses/[YYYY]-[YYYY]_clim_values_by_doy.h5')

RUN TIME: ~4 min per data year
"""

import os, sys, datetime, glob
import h5py as hdf
import numpy as np
from Date_Convert import *

def message(str):
    print str
    sys.stdout.flush()
    return    

def makefname(path,yr,doy):
    mmdd = doy_to_date(yr,doy)
    if mmdd < 1000:
        fname = '%s/%d0%d_NCEI_grids_2.h5' % (path,yr,mmdd)
    else:
        fname = '%s/%d%d_NCEI_grids_2.h5' % (path,yr,mmdd)
    return fname

def write_to_file(file,gvar,gdata):
    if gvar in file.keys():
        del file[gvar]
    file.create_dataset(gvar, data=gdata, dtype=np.float32, compression='gzip')
    message('- %s %s' % (gvar, str(gdata.shape)))
    return

message(' ')
message('process_NCEI_06.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no input grids directory path indicated, using ./grids')
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
doy_begin = 1
doy_end = 365
doys = np.arange(doy_begin,doy_end+1).astype(int)
ndoys = len(doys)
#
years = np.arange(year_begin,year_end+1).astype(int)
nyears = len(years)
#
wxlist = glob.glob('%s/*_NCEI_grids_2.h5' % path)
message('found %d weather derivative grid files' % len(wxlist))
message(' ')
#
infile = wxlist[0]
message('extracting grid information from weather derivatives file %s' % infile)
h5file = hdf.File(infile,'r') 
ncols = np.copy(h5file['grid/ncols'])
nrows = np.copy(h5file['grid/nrows'])
h5file.close()
message(' ')
#
prcp_by_year_doy = np.zeros((nyears,ndoys))
tmin_by_year_doy = np.zeros((nyears,ndoys))
tmax_by_year_doy = np.zeros((nyears,ndoys))
tavg_by_year_doy = np.zeros((nyears,ndoys))
cd_by_year_doy = np.zeros((nyears,ndoys))
cdd_by_year_doy = np.zeros((nyears,ndoys))
gdd_by_year_doy = np.zeros((nyears,ndoys))
#
for j in range(0,nyears):
    year = years[j]
    message('processing grids for %d' % year)
    for i in range(0,ndoys):
        doy = doys[i]
        filename = makefname(path, year, doy)
        h5file = hdf.File(filename,'r') 
        grid_prcp = np.copy(h5file['grid_prcp'])
        prcp_by_year_doy[j,i] = np.mean(grid_prcp)
        grid_tmin = np.copy(h5file['grid_tmin'])
        tmin_by_year_doy[j,i] = np.mean(grid_tmin)
        grid_tmax = np.copy(h5file['grid_tmax'])
        tmax_by_year_doy[j,i] = np.mean(grid_tmax)
        grid_tavg = np.copy(h5file['grid_tavg'])
        tavg_by_year_doy[j,i] = np.mean(grid_tavg)
        grid_cd = np.copy(h5file['chill_d'])
        cd_by_year_doy[j,i] = np.mean(grid_cd)
        grid_cdd = np.copy(h5file['chill_dd'])
        cdd_by_year_doy[j,i] = np.mean(grid_cdd)
        grid_gdd = np.copy(h5file['grow_dd'])
        gdd_by_year_doy[j,i] = np.mean(grid_gdd)
        h5file.close()
message(' ')
#
outfile = '%s/../analyses/%d-%d_clim_values_by_doy.h5' % (path, year_begin, year_end)
message('writing metadata and series information to %s' % outfile)
h5outfile = hdf.File(outfile,'w') 
h5outfile.create_dataset('meta/filename', data=outfile)
h5outfile.create_dataset('meta/created', data=datetime.datetime.now().isoformat())
h5outfile.create_dataset('meta/by', data='M. Garcia, UWisconsin-Madison FWE')
h5outfile.create_dataset('meta/last_updated', data=datetime.datetime.now().isoformat())
h5outfile.create_dataset('meta/at', data='climatological and statistics values by DOY')
message('- 5 metadata items saved')
write_to_file(h5outfile,'prcp_by_year_doy',prcp_by_year_doy)
write_to_file(h5outfile,'tmin_by_year_doy',tmin_by_year_doy)
write_to_file(h5outfile,'tmax_by_year_doy',tmax_by_year_doy)
write_to_file(h5outfile,'tavg_by_year_doy',tavg_by_year_doy)
write_to_file(h5outfile,'cd_by_year_doy',cd_by_year_doy)
write_to_file(h5outfile,'cdd_by_year_doy',cd_by_year_doy)
write_to_file(h5outfile,'gdd_by_year_doy',gdd_by_year_doy)
message('- 6 series arrays saved')
h5outfile.close()
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
h5outfile = hdf.File(outfile,'r+') 
write_to_file(h5outfile,'prcp_mean_by_doy',prcp_mean_by_doy)
write_to_file(h5outfile,'prcp_min_by_doy',prcp_min_by_doy)
write_to_file(h5outfile,'prcp_max_by_doy',prcp_max_by_doy)
write_to_file(h5outfile,'tmin_mean_by_doy',tmin_mean_by_doy)
write_to_file(h5outfile,'tmin_min_by_doy',tmin_min_by_doy)
write_to_file(h5outfile,'tmin_max_by_doy',tmin_max_by_doy)
write_to_file(h5outfile,'tmax_mean_by_doy',tmax_mean_by_doy)
write_to_file(h5outfile,'tmax_min_by_doy',tmax_min_by_doy)
write_to_file(h5outfile,'tmax_max_by_doy',tmax_max_by_doy)
write_to_file(h5outfile,'tavg_mean_by_doy',tavg_mean_by_doy)
write_to_file(h5outfile,'tavg_min_by_doy',tavg_min_by_doy)
write_to_file(h5outfile,'tavg_max_by_doy',tavg_max_by_doy)
write_to_file(h5outfile,'cd_mean_by_doy',cd_mean_by_doy)
write_to_file(h5outfile,'cd_min_by_doy',cd_min_by_doy)
write_to_file(h5outfile,'cd_max_by_doy',cd_max_by_doy)
write_to_file(h5outfile,'cdd_mean_by_doy',cdd_mean_by_doy)
write_to_file(h5outfile,'cdd_min_by_doy',cdd_min_by_doy)
write_to_file(h5outfile,'cdd_max_by_doy',cdd_max_by_doy)
write_to_file(h5outfile,'gdd_mean_by_doy',gdd_mean_by_doy)
write_to_file(h5outfile,'gdd_min_by_doy',gdd_min_by_doy)
write_to_file(h5outfile,'gdd_max_by_doy',gdd_max_by_doy)
message('- 21 series saved')
h5outfile.close()
message(' ')
#
message('process_NCEI_06.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_06.py
