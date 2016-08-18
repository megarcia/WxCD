"""
Python script 'process_NCEI_01.py'
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

USAGE: 'python process_NCEI_01.py NCEI_WLS_20000101-20101231 ./data'

PURPOSE: Extraction of daily meteorological station data from cleaned NOAA/NCEI dataset

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O

INPUT: Output files from process_NCEI_00.py script in '.csv' and '.h5' formats (two files)

OUTPUT: Updated input '.h5' file with sorted meteorological data (one existing, zero new files)

RUN TIME: Depends on NCEI data density, ~30 sec per data year
"""

import os, sys, datetime
import h5py as hdf
import numpy as np
import pandas as pd

def message(str):
    print str
    sys.stdout.flush()
    return    

message(' ')
message('process_NCEI_01.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 3:
    message('input warning: no data directory path indicated, using ./data')
    path = './data'
else:
    path = sys.argv[2]
#
if len(sys.argv) < 2:
    message('input error: need prefix for files containing NCEI weather data')
    sys.exit(1)
else:
    NCEIfname = sys.argv[1]
datafile = '%s/%s_cleaned.csv' % (path, NCEIfname)
h5fname = '%s/%s_processed.h5' % (path, NCEIfname)
#
message('reading station and date information from %s' % h5fname)
h5infile = hdf.File(h5fname,'r')
stn_id = np.copy(h5infile['stn_id'])
dates = np.copy(h5infile['dates'])
h5infile.close()
message('- identifiers for %d stations found' % len(stn_id))
message('- meteorological data for %d dates found' % len(dates))
message(' ')
#
message('loading weather observation information from %s' % datafile)
stndata_df = pd.read_csv(datafile, low_memory=False)
ndatarows, ndatacols = np.shape(stndata_df)
message('- read %d total data rows with %d columns' % (ndatarows, ndatacols))
stndata_df = stndata_df.drop(['Unnamed: 0','IDX'], axis=1)
message('- dropped index columns')
#
# sanity (data integrity) check
metvals = ['STATION','STATION_NAME','ELEVATION','LATITUDE','LONGITUDE','DATE','PRCP','PRCP_M_FLAG','PRCP_Q_FLAG','TMAX','TMAX_M_FLAG','TMAX_Q_FLAG','TMIN','TMIN_M_FLAG','TMIN_Q_FLAG']
idxs = list(stndata_df.columns.values)
if idxs == metvals:
    stndata_df = stndata_df.drop(['PRCP_M_FLAG','PRCP_Q_FLAG','TMAX_M_FLAG','TMAX_Q_FLAG','TMIN_M_FLAG','TMIN_Q_FLAG'], axis=1)
    message('- dropped data flag columns')
else:
    message('input error: cleaned NCEI weather data file does not have the expected fields')
    message('   expected %s' % str(metvals))
    message('  but found %s' % str(idxs))
    sys.exit(1)
message(' ')
#
# sort dataset by date and station, and process met values by date
stndata_df = stndata_df.sort_values(by=['DATE','STATION'])
for date in dates:
    message('gathering met data for %d' % date)
    date_df = stndata_df[stndata_df['DATE'] == date]
    nr,nc = np.shape(date_df)
    message('- found %d total rows' % nr)
    #
    prcp_valid_df = date_df[date_df['PRCP'] != -9999]
    prcp_stns = list(prcp_valid_df['STATION'])
    prcp_lat = np.array(prcp_valid_df['LATITUDE'])
    prcp_lon = np.array(prcp_valid_df['LONGITUDE'])
    prcp_vals = np.array(prcp_valid_df['PRCP']) / 100.0  # convert from 0.1mm to cm
    message('-- %d stns with PRCP data (mean %.1f  stdev %.1f  min %.1f  max %.1f)' % (len(prcp_stns), np.mean(prcp_vals), np.std(prcp_vals), np.min(prcp_vals), np.max(prcp_vals)))
    #
    tmax_valid_df = date_df[date_df['TMAX'] != -9999]
    tmax_stns = list(tmax_valid_df['STATION'])
    tmax_lat = np.array(tmax_valid_df['LATITUDE'])
    tmax_lon = np.array(tmax_valid_df['LONGITUDE'])
    tmax_vals = np.array(tmax_valid_df['TMAX']) / 10.0  # convert from 0.1dC to dC
    message('-- %d stns with TMAX data (mean %.1f  stdev %.1f  min %.1f  max %.1f)' % (len(tmax_stns), np.mean(tmax_vals), np.std(tmax_vals), np.min(tmax_vals), np.max(tmax_vals)))
    #
    tmin_valid_df = date_df[date_df['TMIN'] != -9999]
    tmin_stns = list(tmin_valid_df['STATION'])
    tmin_lat = np.array(tmin_valid_df['LATITUDE'])
    tmin_lon = np.array(tmin_valid_df['LONGITUDE'])
    tmin_vals = np.array(tmin_valid_df['TMIN']) / 10.0  # convert from 0.1dC to dC
    message('-- %d stns with TMIN data (mean %.1f  stdev %.1f  min %.1f  max %.1f)' % (len(tmin_stns), np.mean(tmin_vals), np.std(tmin_vals), np.min(tmin_vals), np.max(tmin_vals)))
    #
    message('- saving %d met data to %s' % (date, h5fname))
    h5file = hdf.File(h5fname,'r+')
    if 'last_updated' in h5file['meta'].keys():
        del h5file['meta/last_updated']
    h5file.create_dataset('meta/last_updated', data=datetime.datetime.now().isoformat())
    if 'at' in h5file['meta'].keys():
        del h5file['meta/at']
    h5file.create_dataset('meta/at', data=date)
    datepath = 'metdata/%d' % date
    if 'metdata' in h5file.keys():
        if date in h5file['metdata'].keys():
            del h5file[datepath]
    h5file.create_dataset(datepath + '/prcp_stns', data=prcp_stns)
    h5file.create_dataset(datepath + '/prcp_lat', data=prcp_lat)
    h5file.create_dataset(datepath + '/prcp_lon', data=prcp_lon)
    h5file.create_dataset(datepath + '/prcp_vals', data=prcp_vals)
    h5file.create_dataset(datepath + '/tmax_stns', data=tmax_stns)
    h5file.create_dataset(datepath + '/tmax_lat', data=tmax_lat)
    h5file.create_dataset(datepath + '/tmax_lon', data=tmax_lon)
    h5file.create_dataset(datepath + '/tmax_vals', data=tmax_vals)
    h5file.create_dataset(datepath + '/tmin_stns', data=tmin_stns)
    h5file.create_dataset(datepath + '/tmin_lat', data=tmin_lat)
    h5file.create_dataset(datepath + '/tmin_lon', data=tmin_lon)
    h5file.create_dataset(datepath + '/tmin_vals', data=tmin_vals)
    h5file.close()
    message(' ')
#
message('process_NCEI_01.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_01.py
