"""
Python script 'process_NCEI_02a.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Splits multi-year processed meteorological data file into individual
         years for separate/parallel use

DEPENDENCIES: h5py, numpy

USAGE: 'python process_NCEI_02a.py NCEI_WLS_19830101-20131231 ./data'

NOTES: <NCEI_WLS_19840101-20131231> is the '_processed.h5' file prefix in your
       'data/' directory

INPUT: '.h5' output from process_NCEI_01.py

OUTPUT: annual '.h5' files with original meteorological data
"""


import sys
import datetime
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
message('process_NCEI_02a.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 3:
    message('input warning: no data directory path indicated, using ./data')
    path = './data'
else:
    path = sys.argv[2]
#
if len(sys.argv) < 2:
    message('input error: need prefix for weather station data h5 file')
    sys.exit(1)
else:
    NCEIfname = sys.argv[1]
h5infname = '%s/%s_processed.h5' % (path, NCEIfname)
parts = NCEIfname.split('_')
#
message('reading station and date information from %s' % h5infname)
with hdf.File(h5infname, 'r') as h5infile:
    stn_id = np.copy(h5infile['stn_id'])
    all_dates = np.copy(h5infile['dates'])
all_years = sorted(list(set([int(j // 1E4) for j in all_dates])))
message('- identifiers for %d stations found' % len(stn_id))
message('- meteorological data for %d total dates found' % len(all_dates))
message(' ')
#
for year in all_years:
    h5outfname = '%s/%s_%s_%d_processed.h5' % (path, parts[0], parts[1], year)
    with hdf.File(h5outfname, 'w') as h5file:
        h5file.create_dataset('meta/filename', data=h5outfname)
        h5file.create_dataset('meta/created',
                              data=datetime.datetime.now().isoformat())
        h5file.create_dataset('meta/by',
                              data='M. Garcia, UWisconsin-Madison FWE')
        h5file.create_dataset('meta/last_updated',
                              data=datetime.datetime.now().isoformat())
        h5file.create_dataset('meta/at',
                              data='unique stations list, dates list')
        h5file.create_dataset('stn_id', data=stn_id)
        h5file.create_dataset('dates', data=all_dates)
    message('initialized processing metadata in %s' % h5outfname)
message(' ')
#
for date in all_dates:
    message('transferring station information and met data for %d' % date)
    with hdf.File(h5infname, 'r') as h5infile:
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
    year = date // 1E4
    h5outfname = '%s/%s_%s_%d_processed.h5' % (path, parts[0], parts[1], year)
    with hdf.File(h5outfname, 'r+') as h5outfile:
        if 'last_updated' in h5outfile['meta'].keys():
            del h5outfile['meta/last_updated']
        h5outfile.create_dataset('meta/last_updated',
                                 data=datetime.datetime.now().isoformat())
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
message(' ')
#
message('process_NCEI_02a.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_02a.py
