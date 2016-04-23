"""
Python script 'process_NCEI_02_preprocess.py'
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

USAGE: 'python process_NCEI_02_preprocess.py NCEI_WLS_19840101-20131231 ./data'

NOTES: [NCEI_WLS_19840101-20131231] is the '_processed.h5' file prefix in your 'data/' directory

PURPOSE: Splits multi-year processed meteorological data file into individual years for parallel use

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O

INPUT: Output file from process_NCEI_01.py script in '.h5' format (one file)
       (with the naming convention 'data/NCEI_WLS_[YYYYMMDD]-[YYYYMMDD]_processed.h5')

OUTPUT: Annual '.h5' files with original meteorological data 

RUN TIME: ~X sec per data year
"""

import os, sys, datetime
import h5py as hdf
import numpy as np

def message(str):
    print str
    sys.stdout.flush()
    return    

message(' ')
message('process_NCEI_02_preprocess.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 3:
    message('input warning: no data directory path indicated, using ./data')
    path = './data'
else:
    path = sys.argv[2]
#
if len(sys.argv) < 2:
    message('input error: need prefix for h5 file containing NCEI weather data')
    sys.exit(1)
else:
    NCEIfname = sys.argv[1]
h5infname = '%s/%s_processed.h5' % (path, NCEIfname)
parts = NCEIfname.split('_')
#
message('reading station and date information from %s' % h5infname)
h5infile = hdf.File(h5infname,'r')
stn_id = np.copy(h5infile['stn_id'])
all_dates = np.copy(h5infile['dates'])
h5infile.close()
all_years = sorted(list(set([int(j//1E4) for j in all_dates])))
message('- identifiers for %d stations found' % len(stn_id))
message('- meteorological data for %d total dates found' % len(all_dates))
message(' ')
#
for year in all_years:
    h5outfname = '%s/%s_%s_%d_processed.h5' % (path, parts[0], parts[1], year)
    h5file = hdf.File(h5outfname,'w')
    h5file.create_dataset('meta/filename', data=h5outfname)
    h5file.create_dataset('meta/created', data=datetime.datetime.now().isoformat())
    h5file.create_dataset('meta/by', data='M. Garcia, UWisconsin-Madison FWE')
    h5file.create_dataset('meta/last_updated', data=datetime.datetime.now().isoformat())
    h5file.create_dataset('meta/at', data='unique stations list, dates list')
    h5file.create_dataset('stn_id', data=stn_id)
    h5file.create_dataset('dates', data=all_dates)
    h5file.close()
    message('initialized processing metadata in %s' % h5outfname)
message(' ')
#
for date in all_dates:
    message('transferring station information and met data for %d' % date)
    h5infile = hdf.File(h5infname,'r')
    datepath = 'metdata/%d' % date
    prcp_stns = np.copy(h5infile['%s/prcp_stns' % datepath])
    prcp_lat = np.copy(h5infile['%s/prcp_lat' % datepath])
    prcp_lon = np.copy(h5infile['%s/prcp_lon' % datepath])
    prcp_vals = np.copy(h5infile['%s/prcp_vals' % datepath])
    tmax_stns = np.copy(h5infile['%s/tmax_stns' % datepath])
    tmax_lat = np.copy(h5infile['%s/tmax_lat' % datepath])
    tmax_lon = np.copy(h5infile['%s/tmax_lon' % datepath])
    tmax_vals = np.copy(h5infile['%s/tmax_vals' % datepath])
    tmin_stns = np.copy(h5infile['%s/tmin_stns' % datepath])
    tmin_lat = np.copy(h5infile['%s/tmin_lat' % datepath])
    tmin_lon = np.copy(h5infile['%s/tmin_lon' % datepath])
    tmin_vals = np.copy(h5infile['%s/tmin_vals' % datepath])
    h5infile.close()
    year = date // 1E4
    h5outfname = '%s/%s_%s_%d_processed.h5' % (path, parts[0], parts[1], year)
    h5outfile = hdf.File(h5outfname,'r+')
    if 'last_updated' in h5outfile['meta'].keys():
        del h5outfile['meta/last_updated']
    h5outfile.create_dataset('meta/last_updated', data=datetime.datetime.now().isoformat())
    if 'at' in h5outfile['meta'].keys():
        del h5outfile['meta/at']
    h5outfile.create_dataset('meta/at', data=date)
    datepath = 'metdata/%d' % date
    if 'metdata' in h5outfile.keys():
        if date in h5outfile['metdata'].keys():
            del h5outfile[datepath]
    h5outfile.create_dataset(datepath + '/prcp_stns', data=prcp_stns)
    h5outfile.create_dataset(datepath + '/prcp_lat', data=prcp_lat)
    h5outfile.create_dataset(datepath + '/prcp_lon', data=prcp_lon)
    h5outfile.create_dataset(datepath + '/prcp_vals', data=prcp_vals)
    h5outfile.create_dataset(datepath + '/tmax_stns', data=tmax_stns)
    h5outfile.create_dataset(datepath + '/tmax_lat', data=tmax_lat)
    h5outfile.create_dataset(datepath + '/tmax_lon', data=tmax_lon)
    h5outfile.create_dataset(datepath + '/tmax_vals', data=tmax_vals)
    h5outfile.create_dataset(datepath + '/tmin_stns', data=tmin_stns)
    h5outfile.create_dataset(datepath + '/tmin_lat', data=tmin_lat)
    h5outfile.create_dataset(datepath + '/tmin_lon', data=tmin_lon)
    h5outfile.create_dataset(datepath + '/tmin_vals', data=tmin_vals)
    h5outfile.close()
message(' ')
#
message('process_NCEI_02_preprocess.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_02_preprocess.py
