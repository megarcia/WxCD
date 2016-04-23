"""
Python module 'Teleconnections.py'
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
               Garcia, M., and P.A. Townsend, in review: "Climatological influences on the 
               forest growing season around western Lake Superior, USA." Submitted to J. 
               Geophys. Res. Atmos. on 5 April 2016.
           See also 'README.md', 'CITATION.txt', and 'ACKNOWLEDGEMENTS.txt' for more information.

USAGE: insert 'from Teleconnections import *; line near head of script, then (for example)
           ice_on_doy, ice_off_doy, ice_duration = get_ice_dates(fname, year_begin, (year_end + 1))
           pdo_winter, pdo_spring, pdo_summer, pdo_autumn = get_ao_idx(fname, year_begin, (year_end + 1))
           nino3, nino4, nino34 = get_nino_idx(fname, year_begin, (year_end + 1))

PURPOSE: For importing and processing teleconnection time series from outside data files

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The 'Date_Convert' module has no external requirements

INPUT: Provided by calling script

OUTPUT: Returned to calling script

RUN TIME: negligible
"""

import numpy as np
import pandas as pd
from Date_Convert import *

def get_ice_dates(fname, yr_begin, yr_end):
    print '- reading %s' % fname
    arr_df = pd.read_csv(fname)
    arr_df = arr_df.drop('source', axis=1)
    arr = np.array(arr_df)
    nrows, ncols = np.shape(arr)
    print '-- %d rows found' % nrows
    yr_min = int(arr[0,0])
    yr_max = int(arr[nrows-1,0])
    yrs = np.arange(yr_min,yr_max+1).astype(int)
    ice_on = np.zeros(np.shape(yrs))
    ice_off = np.zeros(np.shape(yrs))
    ice_duration = np.zeros(np.shape(yrs))
    for i in range(1,nrows):
        if arr[i,1] != -999:
            year_on = arr[i,1]
            date_on = arr[i,2] * 100 + arr[i,3]
            doy_on = date_to_doy(year_on, date_on)
            if arr[i,1] == arr[i,0]:
                ice_on[i] = doy_on
            else:
                ice_on[i] = doy_on - 365
            year_off = arr[i,4]
            date_off = arr[i,5] * 100 + arr[i,6]
            ice_off[i] = date_to_doy(year_off, date_off)
        else:
            ice_on[i] = ice_off[i] = date_to_doy(arr[i,0], 601)
        ice_duration[i] = ice_off[i] - ice_on[i]
    if yr_min < yr_begin:
        ice_on = ice_on[(yr_begin - yr_min):]
        ice_off = ice_off[(yr_begin - yr_min):]
        ice_duration = ice_duration[(yr_begin - yr_min):]
    if yr_end <= yr_max:
        ice_on = ice_on[:(len(ice_on) - (yr_max - yr_end))]
        ice_off = ice_off[:(len(ice_off) - (yr_max - yr_end))]
        ice_duration = ice_duration[:(len(ice_duration) - (yr_max - yr_end))]
    return ice_on, ice_off, ice_duration

def get_nino_idx(fname, yr_begin, yr_end):
    print '- reading %s' % fname
    arr_df = pd.read_csv(fname)
    arr = np.array(arr_df)
    nrows, ncols = np.shape(arr)
    print '-- %d rows found' % nrows
    yr_min = int(arr[1,0]) + 1
    yr_max = int(arr[nrows-1,0]) + 1
    yrs = np.arange(yr_min,yr_max).astype(int)
    nino3_sum = np.zeros(np.shape(yrs))
    nino4_sum = np.zeros(np.shape(yrs))
    nino34_sum = np.zeros(np.shape(yrs))
    target_months_a = [12]
    target_months_b = [1, 2, 3]
    for i in range(1,nrows):
        idx = -1
        if int(arr[i,1]) in target_months_a:
            idx = int(arr[i,0]) + 1 - yr_min
        elif int(arr[i,1]) in target_months_b:
            idx = int(arr[i,0]) - yr_min
        if idx >= 0:
            nino3_sum[idx] += arr[i,(ncols-5)]
            nino4_sum[idx] += arr[i,(ncols-3)]
            nino34_sum[idx] += arr[i,(ncols-1)]
    nmonths = len(target_months_a) + len(target_months_b)
    nino3_idx = nino3_sum / nmonths
    nino4_idx = nino4_sum / nmonths
    nino34_idx = nino34_sum / nmonths
    if yr_min < yr_begin:
        nino3_idx = nino3_idx[(yr_begin - yr_min):]
        nino4_idx = nino4_idx[(yr_begin - yr_min):]
        nino34_idx = nino34_idx[(yr_begin - yr_min):]
    if yr_end <= yr_max:
        nino3_idx = nino3_idx[:(len(nino3_idx) - (yr_max - yr_end - 1))]
        nino4_idx = nino4_idx[:(len(nino4_idx) - (yr_max - yr_end - 1))]
        nino34_idx = nino34_idx[:(len(nino34_idx) - (yr_max - yr_end - 1))]
    return nino3_idx, nino4_idx, nino34_idx

def get_ao_idx(fname, yr_begin, yr_end):
    print '- reading %s' % fname
    arr_df = pd.read_csv(fname)
    arr = np.array(arr_df)
    nrows, ncols = np.shape(arr)
    print '-- %d rows found' % nrows
    yr_min = int(arr[1,0]) + 1
    yr_max = int(arr[nrows-1,0]) + 1
    yrs = np.arange(yr_min,yr_max).astype(int)
    ao_winter_sum = np.zeros(np.shape(yrs))
    ao_spring_sum = np.zeros(np.shape(yrs))
    ao_summer_sum = np.zeros(np.shape(yrs))
    ao_autumn_sum = np.zeros(np.shape(yrs))
    winter_months_a = [12]
    winter_months_b = [1, 2]
    spring_months = [3, 4, 5]
    summer_months = [6, 7, 8]
    autumn_months = [9, 10, 11]
    for i in range(1,nrows):
        if int(arr[i,1]) in winter_months_a:
            idx = int(arr[i,0]) - yr_min + 1
            ao_winter_sum[idx] += arr[i,(ncols-1)]
        elif int(arr[i,1]) in winter_months_b:
            idx = int(arr[i,0]) - yr_min
            ao_winter_sum[idx] += arr[i,(ncols-1)]
        elif int(arr[i,1]) in spring_months:
            idx = int(arr[i,0]) - yr_min
            ao_spring_sum[idx] += arr[i,(ncols-1)]
        elif int(arr[i,1]) in summer_months:
            idx = int(arr[i,0]) - yr_min
            ao_summer_sum[idx] += arr[i,(ncols-1)]
        elif int(arr[i,1]) in autumn_months:
            idx = int(arr[i,0]) - yr_min
            ao_autumn_sum[idx] += arr[i,(ncols-1)]
    winter_nmonths = len(winter_months_a) + len(winter_months_b)
    spring_nmonths = len(spring_months)
    summer_nmonths = len(summer_months)
    autumn_nmonths = len(autumn_months)
    ao_winter_idx = ao_winter_sum / winter_nmonths
    ao_spring_idx = ao_spring_sum / spring_nmonths
    ao_summer_idx = ao_summer_sum / summer_nmonths
    ao_autumn_idx = ao_autumn_sum / autumn_nmonths
    if yr_min < yr_begin:
        ao_winter_idx = ao_winter_idx[(yr_begin - yr_min):]
        ao_spring_idx = ao_spring_idx[(yr_begin - yr_min):]
        ao_summer_idx = ao_summer_idx[(yr_begin - yr_min):]
        ao_autumn_idx = ao_autumn_idx[(yr_begin - yr_min):]
    if yr_end <= yr_max:
        ao_winter_idx = ao_winter_idx[:(len(ao_winter_idx) - (yr_max - yr_end - 1))]
        ao_spring_idx = ao_spring_idx[:(len(ao_spring_idx) - (yr_max - yr_end - 1))]
        ao_summer_idx = ao_summer_idx[:(len(ao_summer_idx) - (yr_max - yr_end - 1))]
        ao_autumn_idx = ao_autumn_idx[:(len(ao_autumn_idx) - (yr_max - yr_end - 1))]
    return ao_winter_idx, ao_spring_idx, ao_summer_idx, ao_autumn_idx

# end Teleconnections.py
