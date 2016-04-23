"""
Python script 'process_NCEI_11.py'
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

USAGE: 'python process_NCEI_11.py ecoregion_polygonIDs.txt 1984 2013 ./analyses'

PURPOSE: Summarize climatological statistics for individual ecoregions and perform various
         statistical tests on time series across ecoregions for each variable:
         - Pearson correlation and p-value (similarity between time series)
         - T-test and p-value (difference between time series means)
         - Levene test and p-value (difference between time series variances)
         Calculate similarity (analogous to Euclidean distance in x-y-z space) for each 
         variable and ecoregion pair for use in 'process_NCEI_12.py' clustering procedure.

NOTE: The selected similarity measure uses R and the p-values from all three time series 
      statistical test results. See the reference listed above for more discussion on that 
      formulation. A few other formulations are given here and commented out.

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O

INPUT: Ecoregion time series and stats in '.h5' files from 'process_NCEI_10.py'

OUTPUT: Ecoregion-based time series correlation statistics in '.h5' files

RUN TIME: ~1 min
"""

import os, sys, datetime
import h5py as hdf
import numpy as np
from scipy.stats import pearsonr, ttest_ind, levene

def message(str):
    print str
    sys.stdout.flush()
    return    

message(' ')
message('process_NCEI_11.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
# 
if len(sys.argv) < 5:
    message('input warning: no input analyses directory path indicated, using ./analyses')
    path = './analyses'
else:
    path = sys.argv[4]
# 
if len(sys.argv) < 4:
    year_begin = 1984
    year_end = 2013 
else:
    year_begin = int(sys.argv[2])      
    year_end = int(sys.argv[3])      
#
if len(sys.argv) < 2:
    message('input error: need polygon index/ID txt file')
    sys.exit(1)
else:
	ecopolysfname = '%s/../data/%s' % (path, sys.argv[1])
#
tsfname = '%s/%d-%d_ecoregion_timeseries.h5' % (path,year_begin,year_end)
statsfname = '%s/%d-%d_ecoregion_stats.h5' % (path,year_begin,year_end)
years = np.arange(year_begin,year_end+1).astype(int)
nyears = len(years)
#
message('getting ecoregion designations and variable names from %s' % tsfname)
h5infile = hdf.File(tsfname,'r')
ecoregion_IDs = np.copy(h5infile['ecoregion_IDs'])
necoregions = len(ecoregion_IDs)
varnames = np.copy(h5infile['varnames'])
nvars = len(varnames)
h5infile.close()
#
message('retrieving stats by ecoregion from %s' % statsfname)
stats = np.zeros((necoregions,nvars,7))
ecoregion_area = np.zeros((necoregions))
h5infile = hdf.File(statsfname,'r')
for i in range(0,necoregions):
    message('- %s' % ecoregion_IDs[i])
    stats[i,:,:] = np.copy(h5infile[ecoregion_IDs[i]])
    ecoregion_area[i] = stats[i,0,0]
h5infile.close()
message(' ')
#
message('reorganizing stats by variable and saving to %s' % statsfname)
h5outfile = hdf.File(statsfname,'r+')
for i in range(1,nvars):
    stats_out = np.zeros((necoregions,7))
    message('- %s' % varnames[i])
    stats_out[:,:] = stats[:,i,:]
    if varnames[i] in h5outfile.keys():
        del h5outfile[varnames[i]]
    h5outfile.create_dataset(varnames[i], data=stats_out, dtype=np.float32, compression='gzip')
    csvname = '%s/%d-%d_%s_stats_by_ecoregion.csv' % (path,year_begin,year_end,varnames[i])
    np.savetxt(csvname, stats_out, delimiter=',')
h5outfile.close()
message(' ')
#
message('retrieving time series by ecoregion from %s' % tsfname)
timeseries = np.zeros((necoregions,nvars,nyears))
h5infile = hdf.File(tsfname,'r')
for i in range(0,necoregions):
    message('- %s' % ecoregion_IDs[i])
    timeseries[i,:,:] = np.copy(h5infile[ecoregion_IDs[i]])
h5infile.close()
message(' ')
#
message('reorganizing time series by variable and saving to %s' % tsfname)
h5outfile = hdf.File(tsfname,'r+')
for i in range(1,nvars):
    timeseries_out = np.zeros((necoregions,nyears))
    message('- %s' % varnames[i])
    timeseries_out[:,:] = timeseries[:,i,:]
    if varnames[i] in h5outfile.keys():
        del h5outfile[varnames[i]]
    h5outfile.create_dataset(varnames[i], data=timeseries_out, dtype=np.float32, compression='gzip')
    csvname = '%s/%d-%d_%s_timeseries_by_ecoregion.csv' % (path,year_begin,year_end,varnames[i])
    np.savetxt(csvname, timeseries_out, delimiter=',')
h5outfile.close()
message(' ')
#
np.seterr(invalid='ignore')
message('significant time series correlations and differences (p < 0.05)')
message(' ')
for k in range(1,nvars):
    message('- variable %s' % varnames[k])
    ts_corr_matrix = np.zeros((necoregions,necoregions))
    ts_rsig_matrix = np.zeros((necoregions,necoregions))
    ts_ttest_matrix = np.zeros((necoregions,necoregions))
    ts_tsig_matrix = np.zeros((necoregions,necoregions))
    ts_levene_matrix = np.zeros((necoregions,necoregions))
    ts_lsig_matrix = np.zeros((necoregions,necoregions))
    ts_dist_matrix = np.zeros((necoregions,necoregions))
    for i in range(0,necoregions):
        for j in range(i,necoregions):
            ts_corr_matrix[i,j], ts_rsig_matrix[i,j] = pearsonr(timeseries[i,k,:], timeseries[j,k,:])
            ts_ttest_matrix[i,j], ts_tsig_matrix[i,j] = ttest_ind(timeseries[i,k,:], timeseries[j,k,:], equal_var=False)
            ts_levene_matrix[i,j], ts_lsig_matrix[i,j] = levene(timeseries[i,k,:], timeseries[j,k,:])
            # ts_dist_matrix[i,j] = ts_corr_matrix[i,j] * (1.0 - ts_rsig_matrix[i,j]) * (1.0 - ts_tsig_matrix[i,j]) * (1.0 - ts_lsig_matrix[i,j])
            # ts_dist_matrix[i,j] = np.sqrt((1.0-ts_corr_matrix[i,j])**2 + ts_rsig_matrix[i,j]**2 + ts_tsig_matrix[i,j]**2 + ts_lsig_matrix[i,j]**2)
            # ts_dist_matrix[i,j] = np.sqrt(ts_rsig_matrix[i,j]**2 + (1.0 - ts_tsig_matrix[i,j])**2 + (1.0 - ts_lsig_matrix[i,j])**2)
            ts_dist_matrix[i,j] = np.sqrt((1 - ts_corr_matrix[i,j])**2 + ts_rsig_matrix[i,j]**2 + (1.0 - ts_tsig_matrix[i,j])**2 + (1.0 - ts_lsig_matrix[i,j])**2)
    ts_dist_matrix_nan = np.where(ts_dist_matrix == 0.0, np.nan, ts_dist_matrix)
    dist_mean = np.nanmean(ts_dist_matrix_nan)
    dist_std = np.nanstd(ts_dist_matrix_nan)
    dist_lower_threshold = dist_mean - 1.0 * dist_std
    dist_upper_threshold = dist_mean + 1.0 * dist_std
    message('-- Distance mean = %.3f  std = %.3f  lower threshold = %.3f  upper threshold = %.3f' % (dist_mean,dist_std,dist_lower_threshold,dist_upper_threshold))
    for i in range(0,necoregions):
        for j in range(i,necoregions):
            if (ts_dist_matrix[i,j] < dist_lower_threshold) and (j != i):
                message('-- Similarity: ecoregions %s and %s are similar R = %.3f Rp = %.3f Tp = %.3f Lp = %.3f D = %.3f' % (ecoregion_IDs[i],ecoregion_IDs[j],ts_corr_matrix[i,j],ts_rsig_matrix[i,j],ts_tsig_matrix[i,j],ts_lsig_matrix[i,j],ts_dist_matrix[i,j]))
    for i in range(0,necoregions):
        for j in range(i,necoregions):
            if (ts_dist_matrix[i,j] > dist_upper_threshold) and (j != i):
                message('-- Dissimilarity: ecoregions %s and %s are dissimilar R = %.3f Rp = %.3f Tp = %.3f Lp = %.3f D = %.3f' % (ecoregion_IDs[i],ecoregion_IDs[j],ts_corr_matrix[i,j],ts_rsig_matrix[i,j],ts_tsig_matrix[i,j],ts_lsig_matrix[i,j],ts_dist_matrix[i,j]))
    message('- saving time series cross-testing matrices to %s' % tsfname)
    h5outfile = hdf.File(tsfname,'r+')
    vname = '%s_corr_matrix' % varnames[k]
    if vname in h5outfile.keys():
        del h5outfile[vname]
    h5outfile.create_dataset(vname, data=ts_corr_matrix, dtype=np.float32, compression='gzip')
    message('-- %s' % vname)
    vname = '%s_rsig_matrix' % varnames[k]
    if vname in h5outfile.keys():
        del h5outfile[vname]
    h5outfile.create_dataset(vname, data=ts_rsig_matrix, dtype=np.float32, compression='gzip')
    message('-- %s' % vname)
    vname = '%s_ttest_matrix' % varnames[k]
    if vname in h5outfile.keys():
        del h5outfile[vname]
    h5outfile.create_dataset(vname, data=ts_ttest_matrix, dtype=np.float32, compression='gzip')
    message('-- %s' % vname)
    vname = '%s_tsig_matrix' % varnames[k]
    if vname in h5outfile.keys():
        del h5outfile[vname]
    h5outfile.create_dataset(vname, data=ts_tsig_matrix, dtype=np.float32, compression='gzip')
    message('-- %s' % vname)
    vname = '%s_levene_matrix' % varnames[k]
    if vname in h5outfile.keys():
        del h5outfile[vname]
    h5outfile.create_dataset(vname, data=ts_levene_matrix, dtype=np.float32, compression='gzip')
    message('-- %s' % vname)
    vname = '%s_lsig_matrix' % varnames[k]
    if vname in h5outfile.keys():
        del h5outfile[vname]
    h5outfile.create_dataset(vname, data=ts_lsig_matrix, dtype=np.float32, compression='gzip')
    message('-- %s' % vname)
    vname = '%s_dist_matrix' % varnames[k]
    if vname in h5outfile.keys():
        del h5outfile[vname]
    h5outfile.create_dataset(vname, data=ts_dist_matrix, dtype=np.float32, compression='gzip')
    message('-- %s' % vname)
    h5outfile.close()
    message(' ')
#
message('process_NCEI_11.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_11.py
