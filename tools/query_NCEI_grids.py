"""
Python script 'query_NCEI_grids.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
Treat others as you would be treated. Pay it forward. Valar dohaeris.

USAGE: '$ python query_NCEI_grids.py locations_dates.csv'

PURPOSE: Query grids of climatological derivatives at given locations and dates

DEPENDENCIES: Some standard libraries/modules
              The h5py module is required for handling of HDF5 files
              The 'UTM_Geo_Convert' module has its own requirements

INPUT: A '.csv' file containing lat/lon/date of each query
       Output files from process_NCEI_03.py script in '.h5' format
       (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_2.h5')

OUTPUT: New '.csv' file with original location/date input + several new columns

TO DO: add more returned variables
"""


import sys
import datetime
import glob
import h5py as hdf
import numpy as np
import pandas as pd
from UTM_Geo_Convert import geographic_to_utm


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


def pad(x):
    if x < 10:
        char_string = '0%d' % x
    else:
        char_string = '%d' % x
    return char_string


message(' ')
message('query_NCEI_grids.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 2:
    message('need CSV file containing locations and dates for query')
    sys.exit(1)
else:
    csvfile = sys.argv[1]
#
wxlist = glob.glob('grids/*_NCEI_grids_2.h5')
message('found %d weather derivatives files' % len(wxlist))
message(' ')
#
message('reading %s' % csvfile)
arr_df = pd.read_csv(csvfile)
query_lats = arr_df['LAT'].tolist()
query_lons = arr_df['LONG'].tolist()
query_dates = arr_df['DATE'].tolist()
message('- %d query rows found' % len(query_dates))
message(' ')
#
query_D = []
for i in range(0, len(query_dates)):
    j = query_dates[i].split('/')
    if int(j[2]) <= 50:
        yyyy = int(j[2]) + 2000
    else:
        yyyy = int(j[2]) + 1950
    mm = int(j[0])
    dd = int(j[1])
    datestr = '%d%s%s' % (yyyy, str(mm).zfill(2), str(dd).zfill(2))
    query_D.append(datestr)
#
query_E = []
query_N = []
UTMzone = 15
for i in range(0, len(query_lats)):
    UTM_coords = geographic_to_utm(query_lons[i].astype(np.float64),
                                   query_lats[i].astype(np.float64),
                                   UTMzone)
    query_E.append(int(round(UTM_coords[1], 0)))
    query_N.append(int(round(UTM_coords[2], 0)))
#
cd = []
cdd = []
gdd = []
gdd_base0 = []
p_30 = []
p_90 = []
p_180 = []
p_365 = []
for i in range(0, len(query_D)):
    datestr = query_D[i]
    message('lat %.3f lon %.3f date %s  --> UTM zone %d E %d N %d date %s' %
            (query_lats[i], query_lons[i], query_dates[i], UTMzone, query_E[i],
             query_N[i], datestr))
    wx_fname = 'grids/%s_NCEI_grids_2.h5' % datestr
    if wx_fname in wxlist:
        message('- extracting grid info from weather derivatives file %s' %
                wx_fname)
        with hdf.File(wx_fname, 'r') as h5file:
            wx_SEnorthing = np.copy(h5file['grid/min_y'])
            wx_dy = np.copy(h5file['grid/dy'])
            wx_nrows = np.copy(h5file['grid/nrows'])
            wx_NWeasting = np.copy(h5file['grid/min_x'])
            wx_dx = np.copy(h5file['grid/dx'])
            wx_ncols = np.copy(h5file['grid/ncols'])
            #
            query_row = ((query_N[i] - wx_SEnorthing) // wx_dy)
            if query_row < 0:
                message('- location is south of available grid limits')
                continue
            if query_row > wx_nrows:
                message('- location is north of available grid limits')
                continue
            query_col = ((query_E[i] - wx_NWeasting) // wx_dx)
            if query_col < 0:
                message('- location is west of available grid limits')
                continue
            if query_col > wx_ncols:
                message('- location is east of available grid limits')
                continue
            #
            chill_d = np.copy(h5file['chill_d'])
            chill_dd = np.copy(h5file['chill_dd'])
            grow_dd = np.copy(h5file['grow_dd'])
            grow_dd_base0 = np.copy(h5file['grow_dd_base0'])
            prcp_30d = np.copy(h5file['prcp_30d_sum'])
            prcp_90d = np.copy(h5file['prcp_90d_sum'])
            prcp_180d = np.copy(h5file['prcp_180d_sum'])
            prcp_365d = np.copy(h5file['prcp_365d_sum'])
        #
        if int(datestr[4:8] < 701):
            cd.append(chill_d[query_row, query_col])
            cdd.append(chill_dd[query_row, query_col])
        else:
            wx_fname = 'grids/%s%s_NCEI_grids_2.h5' % (datestr[0:4], '0630')
            with hdf.File(wx_fname, 'r') as h5file:
                chill_d_0630 = np.copy(h5file['chill_d'])
                chill_dd_0630 = np.copy(h5file['chill_dd'])
            cd.append(chill_d_0630[query_row, query_col]
                      + chill_d[query_row, query_col])
            cdd.append(chill_dd_0630[query_row, query_col]
                       + chill_dd[query_row, query_col])
        gdd.append(grow_dd[query_row, query_col])
        gdd_base0.append(grow_dd_base0[query_row, query_col])
        p_30.append(prcp_30d[query_row, query_col])
        p_90.append(prcp_90d[query_row, query_col])
        p_180.append(prcp_180d[query_row, query_col])
        p_365.append(prcp_365d[query_row, query_col])
    else:
        cd.append('NA')
        cdd.append('NA')
        gdd.append('NA')
        gdd_base0.append('NA')
        p_30.append('NA')
        p_90.append('NA')
        p_180.append('NA')
        p_365.append('NA')
message(' ')
#
arr_df['UTM15_E'] = query_E
arr_df['UTM15_N'] = query_N
arr_df['CD'] = cd
arr_df['CDD'] = cdd
arr_df['GDD_base0'] = gdd_base0
arr_df['GDD_base5'] = gdd
arr_df['P_30'] = p_30
arr_df['P_90'] = p_90
arr_df['P_180'] = p_180
arr_df['P_365'] = p_365
#
out_csvfile = csvfile[:-4] + '_Clim.csv'
arr_df.to_csv(out_csvfile)
message('wrote %s' % out_csvfile)
message(' ')
#
message('query_NCEI_grids.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end query_NCEI_grids.py
