"""
Python script 'process_NCEI_14.py'
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

USAGE: 'python process_NCEI_14.py 1984 2013 ./analyses'

PURPOSE: Summarize climatological statistics for ecoregion clusters and perform various
         statistical tests on time series across clusters for each variable:
         - Pearson correlation and p-value (similarity between time series)
         - T-test and p-value (difference between time series means)
         - Levene test and p-value (difference between time series variances)

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O

INPUT: Ecoregion cluster time series and stats in '.h5' files from 'process_NCEI_13.py'

OUTPUT: Aggregated ecoregion cluster-based statistics in '.h5' and '.csv' files

RUN TIME: ~5 sec
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
message('process_NCEI_14.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no input analyses directory path indicated, using ./analyses')
    path = './analyses'
else:
    path = sys.argv[3]
#
if len(sys.argv) < 3:
    year_begin = 1984
    year_end = 2013 
else:
    year_begin = int(sys.argv[1])      
    year_end = int(sys.argv[2])      
#
tsfname = '%s/%d-%d_cluster_timeseries.h5' % (path,year_begin,year_end)
statsfname = '%s/%d-%d_cluster_stats.h5' % (path,year_begin,year_end)
meansdfname = '%s/%d-%d_cluster_meansd.csv' % (path,year_begin,year_end)
trendsfname = '%s/%d-%d_cluster_trends.csv' % (path,year_begin,year_end)
years = np.arange(year_begin,year_end+1).astype(int)
nyears = len(years)
#
message('getting variable names and cluster stats from %s' % statsfname)
h5infile = hdf.File(statsfname,'r')
varnames = np.copy(h5infile['varnames'])
nvars = len(varnames)
nclusters = np.copy(h5infile['nclusters'])
cluster_IDs = []
stats = np.zeros((nclusters,nvars,7))
cluster_area = np.zeros((nclusters))
for i in range(0,nclusters):
    cluster_IDs.append('cluster_%d' % (i+1))
    message('- %s' % cluster_IDs[i])
    stats[i,:,:] = np.copy(h5infile[cluster_IDs[i]])
    cluster_area[i] = stats[i,0,0]
h5infile.close()
message(' ')
#
message('reorganizing stats by variable and saving to %s' % statsfname)
h5outfile = hdf.File(statsfname,'r+')
msf = open(meansdfname,'w')
tf = open(trendsfname,'w')
for i in range(1,nvars):
    stats_out = np.zeros((nclusters,7))
    message('- %s' % varnames[i])
    stats_out[:,:] = stats[:,i,:]
    if varnames[i] in h5outfile.keys():
        del h5outfile[varnames[i]]
    h5outfile.create_dataset(varnames[i], data=stats_out, dtype=np.float32, compression='gzip')
    csvname = '%s/%d-%d_%s_stats_by_cluster.csv' % (path,year_begin,year_end,varnames[i])
    np.savetxt(csvname, stats_out, delimiter=',')
    #
    line = '%s, ' % varnames[i]
    for j in range(0,nclusters):
        line += '%.1f %s %.1f' % (stats_out[j,0], r'$\pm$', stats_out[j,1])
        if j < (nclusters - 1):
            line += ', '
        else:
            line += '\n'
    msf.write(line)
    #
    line = '%s, ' % varnames[i]
    for j in range(0,nclusters):
        if stats_out[j,6] < 0.001:
            sigstr = '***'
        elif stats_out[j,6] < 0.01:
            sigstr = '**'
        elif stats_out[j,6] < 0.05:
            sigstr = '*'
        else:
            sigstr = ''
        line += '%.3f%s' % (stats_out[j,4], sigstr)
        if j < (nclusters - 1):
            line += ', '
        else:
            line += '\n'
    tf.write(line)
h5outfile.close()
msf.close()
tf.close()
message(' ')
#
message('retrieving time series by cluster from %s' % tsfname)
timeseries = np.zeros((nclusters,nvars,nyears))
h5infile = hdf.File(tsfname,'r')
for i in range(0,nclusters):
    message('- %s' % cluster_IDs[i])
    timeseries[i,:,:] = np.copy(h5infile[cluster_IDs[i]])
h5infile.close()
message(' ')
#
message('reorganizing time series by variable and saving to %s' % tsfname)
h5outfile = hdf.File(tsfname,'r+')
for i in range(1,nvars):
    timeseries_out = np.zeros((nclusters,nyears))
    message('- %s' % varnames[i])
    timeseries_out[:,:] = timeseries[:,i,:]
    if varnames[i] in h5outfile.keys():
        del h5outfile[varnames[i]]
    h5outfile.create_dataset(varnames[i], data=timeseries_out, dtype=np.float32, compression='gzip')
    csvname = '%s/%d-%d_%s_timeseries_by_cluster.csv' % (path,year_begin,year_end,varnames[i])
    np.savetxt(csvname, timeseries_out, delimiter=',')
h5outfile.close()
message(' ')
#
np.seterr(invalid='ignore')
message('significant time series correlations and differences (p < 0.05)')
var_dist_matrix = np.zeros((nvars,nclusters,nclusters))
for k in range(1,nvars):
    message('- variable %s' % varnames[k])
    ts_corr_matrix = np.zeros((nclusters,nclusters))
    ts_rsig_matrix = np.zeros((nclusters,nclusters))
    ts_ttest_matrix = np.zeros((nclusters,nclusters))
    ts_tsig_matrix = np.zeros((nclusters,nclusters))
    ts_levene_matrix = np.zeros((nclusters,nclusters))
    ts_lsig_matrix = np.zeros((nclusters,nclusters))
    ts_dist_matrix = np.zeros((nclusters,nclusters))
    for i in range(0,nclusters):
        for j in range(i,nclusters):
            ts_corr_matrix[i,j], ts_rsig_matrix[i,j] = pearsonr(timeseries[i,k,:], timeseries[j,k,:])
            ts_ttest_matrix[i,j], ts_tsig_matrix[i,j] = ttest_ind(timeseries[i,k,:], timeseries[j,k,:], equal_var=False)
            ts_levene_matrix[i,j], ts_lsig_matrix[i,j] = levene(timeseries[i,k,:], timeseries[j,k,:])
            ts_dist_matrix[i,j] = np.sqrt(ts_rsig_matrix[i,j]**2 + (1.0 - ts_tsig_matrix[i,j])**2 + (1.0 - ts_lsig_matrix[i,j])**2)
            if j != i:
                # if ts_rsig_matrix[i,j] < 0.05:
                #     message('-- Correlation: %s and %s are significantly correlated R = %.3f at Rp = %.3f' % (cluster_IDs[i],cluster_IDs[j],ts_corr_matrix[i,j],ts_rsig_matrix[i,j]))
                if ts_corr_matrix[i,j] < 0.0:
                    message('-- Anticorrelation: %s and %s are oppositely correlated R = %.3f at Rp = %.3f' % (cluster_IDs[i],cluster_IDs[j],ts_corr_matrix[i,j],ts_rsig_matrix[i,j]))
                if ts_rsig_matrix[i,j] > 0.05:
                    message('-- Uncorrelation: %s and %s are NOT significantly correlated R = %.3f at Rp = %.3f' % (cluster_IDs[i],cluster_IDs[j],ts_corr_matrix[i,j],ts_rsig_matrix[i,j]))
                if ts_tsig_matrix[i,j] < 0.05:
                    message('-- T-test: %s and %s are significantly different T = %.3f at Tp = %.3f' % (cluster_IDs[i],cluster_IDs[j],ts_ttest_matrix[i,j],ts_tsig_matrix[i,j]))
                if ts_lsig_matrix[i,j] < 0.05:
                    message('-- Levene test: %s and %s have significantly different variances L = %.3f at Lp = %.3f' % (cluster_IDs[i],cluster_IDs[j],ts_levene_matrix[i,j],ts_lsig_matrix[i,j]))
    var_dist_matrix[k,:,:] = ts_dist_matrix[:,:]
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
message('Inter-cluster similarity assessment')
dist_tuples = []
overall_dist_matrix = np.mean(var_dist_matrix, axis=0)
overall_dist_matrix_nan = np.where(overall_dist_matrix == 0.0, np.nan, overall_dist_matrix)
dist_mean = np.nanmean(overall_dist_matrix_nan)
dist_std = np.nanstd(overall_dist_matrix_nan)
dist_lower_threshold = dist_mean - 0.0 * dist_std
dist_upper_threshold = dist_mean + 0.0 * dist_std
message('- Distance overall mean = %.3f std = %.3f lower threshold = %.3f upper threshold = %.3f' % (dist_mean,dist_std,dist_lower_threshold,dist_upper_threshold))
for i in range(0,nclusters):
    for j in range(i,nclusters):
        if (overall_dist_matrix[i,j] < dist_lower_threshold) and (j != i):
            message('- Similarity: %s and %s are overall similar with D = %.3f' % (cluster_IDs[i],cluster_IDs[j],overall_dist_matrix[i,j]))
            dist_tuples.append((overall_dist_matrix[i,j],cluster_IDs[i],cluster_IDs[j]))
for i in range(0,nclusters):
    for j in range(i,nclusters):
        if (overall_dist_matrix[i,j] > dist_upper_threshold) and (j != i):
            message('- Dissimilarity: %s and %s are overall dissimilar with D = %.3f' % (cluster_IDs[i],cluster_IDs[j],overall_dist_matrix[i,j]))
            dist_tuples.append((overall_dist_matrix[i,j],cluster_IDs[i],cluster_IDs[j]))
message(' ')
dist_sorted = sorted(dist_tuples, key=lambda distance: distance[0])
message('Cluster pairs sorted by dissimilarity (distance)')
for pair in dist_sorted:
    message('  %.3f %s %s' % (pair[0],pair[1],pair[2]))
message(' ')
#
message('process_NCEI_14.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_14.py
