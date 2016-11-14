"""
Python script 'process_NCEI_11.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Summarize climatological statistics for individual ecoregions
         Perform various statistical tests on time series across ecoregions
            for each variable:
         - Pearson correlation and p-value (similarity between time series)
         - T-test and p-value (difference between time series means)
         - Levene test and p-value (difference between time series variances)
         Calculate similarity (analogous to Euclidean distance in x-y-z space)
            for each variable + ecoregion pair for use in 'process_NCEI_12.py'
            clustering procedure.

NOTE: The selected similarity measure uses R and the p-values from all three
        time series statistical test results. See the reference listed above
        for more discussion on that formulation. A few of the formulations that
        we tried are also given here and commented out.

DEPENDENCIES: h5py, numpy, scipy

USAGE: '$ python process_NCEI_11.py ecoregion_polygonIDs.txt 1984 2013
          ./analyses'

INPUT: '.h5' files from 'process_NCEI_10.py'

OUTPUT: '.h5' ecoregion-based time series correlation statistics
"""


import sys
import datetime
import h5py as hdf
import numpy as np
from scipy.stats import pearsonr, ttest_ind, levene


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


message(' ')
message('process_NCEI_11.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 5:
    message('input warning: no directory indicated, using ./analyses')
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
tsfname = '%s/%d-%d_ecoregion_timeseries.h5' % (path, year_begin, year_end)
statsfname = '%s/%d-%d_ecoregion_stats.h5' % (path, year_begin, year_end)
#
years = np.arange(year_begin, year_end + 1).astype(int)
nyears = len(years)
#
message('getting ecoregion designations and variable names from %s' % tsfname)
with hdf.File(tsfname, 'r') as h5infile:
    ecoregion_IDs = np.copy(h5infile['ecoregion_IDs'])
    varnames = np.copy(h5infile['varnames'])
necoregions = len(ecoregion_IDs)
nvars = len(varnames)
#
message('retrieving stats by ecoregion from %s' % statsfname)
stats = np.zeros((necoregions, nvars, 7))
ecoregion_area = np.zeros((necoregions))
with hdf.File(statsfname, 'r') as h5infile:
    for i, ecoregion_ID in enumerate(ecoregion_IDs):
        message('- %s' % ecoregion_ID)
        stats[i, :, :] = np.copy(h5infile[ecoregion_ID])
        ecoregion_area[i] = stats[i, 0, 0]
message(' ')
#
message('reorganizing stats by variable and saving to %s' % statsfname)
with hdf.File(statsfname, 'r+') as h5outfile:
    for i in range(1, nvars):
        stats_out = np.zeros((necoregions, 7))
        message('- %s' % varnames[i])
        stats_out[:, :] = stats[:, i, :]
        if varnames[i] in h5outfile.keys():
            del h5outfile[varnames[i]]
        h5outfile.create_dataset(varnames[i], data=stats_out,
                                 dtype=np.float32, compression='gzip')
        csvname = '%s/%d-%d_%s_stats_by_ecoregion.csv' % \
            (path, year_begin, year_end, varnames[i])
        np.savetxt(csvname, stats_out, delimiter=',')
h5outfile.close()
message(' ')
#
message('retrieving time series by ecoregion from %s' % tsfname)
timeseries = np.zeros((necoregions, nvars, nyears))
with hdf.File(tsfname, 'r') as h5infile:
    for i, ecoregion_ID in enumerate(ecoregion_IDs):
        message('- %s' % ecoregion_ID)
        timeseries[i, :, :] = np.copy(h5infile[ecoregion_ID])
message(' ')
#
message('reorganizing time series by variable and saving to %s' % tsfname)
with hdf.File(tsfname, 'r+') as h5outfile:
    for i in range(1, nvars):
        timeseries_out = np.zeros((necoregions, nyears))
        message('- %s' % varnames[i])
        timeseries_out[:, :] = timeseries[:, i, :]
        if varnames[i] in h5outfile.keys():
            del h5outfile[varnames[i]]
        h5outfile.create_dataset(varnames[i], data=timeseries_out,
                                 dtype=np.float32, compression='gzip')
        csvname = '%s/%d-%d_%s_timeseries_by_ecoregion.csv' % \
            (path, year_begin, year_end, varnames[i])
        np.savetxt(csvname, timeseries_out, delimiter=',')
message(' ')
#
np.seterr(invalid='ignore')
message('significant time series correlations and differences (p < 0.05)')
message(' ')
for k in range(1, nvars):
    message('- variable %s' % varnames[k])
    ts_corr_matrix = np.zeros((necoregions, necoregions))
    ts_rsig_matrix = np.zeros((necoregions, necoregions))
    ts_ttest_matrix = np.zeros((necoregions, necoregions))
    ts_tsig_matrix = np.zeros((necoregions, necoregions))
    ts_levene_matrix = np.zeros((necoregions, necoregions))
    ts_lsig_matrix = np.zeros((necoregions, necoregions))
    ts_dist_matrix = np.zeros((necoregions, necoregions))
    for i in range(necoregions):
        for j in range(i, necoregions):
            ts_corr_matrix[i, j], ts_rsig_matrix[i, j] = \
                pearsonr(timeseries[i, k, :], timeseries[j, k, :])
            ts_ttest_matrix[i, j], ts_tsig_matrix[i, j] = \
                ttest_ind(timeseries[i, k, :], timeseries[j, k, :],
                          equal_var=False)
            ts_levene_matrix[i, j], ts_lsig_matrix[i, j] = \
                levene(timeseries[i, k, :], timeseries[j, k, :])
            """
            ts_dist_matrix[i, j] = \
                ts_corr_matrix[i, j] * (1.0 - ts_rsig_matrix[i, j]) * \
                (1.0 - ts_tsig_matrix[i, j]) * (1.0 - ts_lsig_matrix[i, j])
            ts_dist_matrix[i, j] = \
                np.sqrt((1.0 - ts_corr_matrix[i, j])**2 +
                        ts_rsig_matrix[i, j]**2 + ts_tsig_matrix[i, j]**2 +
                        ts_lsig_matrix[i, j]**2)
            ts_dist_matrix[i, j] = \
                np.sqrt(ts_rsig_matrix[i, j]**2 +
                        (1.0 - ts_tsig_matrix[i, j])**2 +
                        (1.0 - ts_lsig_matrix[i, j])**2)
            """
            ts_dist_matrix[i, j] = \
                np.sqrt((1 - ts_corr_matrix[i, j])**2 +
                        ts_rsig_matrix[i, j]**2 +
                        (1.0 - ts_tsig_matrix[i, j])**2 +
                        (1.0 - ts_lsig_matrix[i, j])**2)
    ts_dist_matrix_nan = np.where(ts_dist_matrix == 0.0,
                                  np.nan, ts_dist_matrix)
    dist_mean = np.nanmean(ts_dist_matrix_nan)
    dist_std = np.nanstd(ts_dist_matrix_nan)
    dist_lower_threshold = dist_mean - 1.0 * dist_std
    dist_upper_threshold = dist_mean + 1.0 * dist_std
    message('-- Distance mean = %.3f  std = %.3f  lower threshold = %.3f  \
            upper threshold = %.3f' %
            (dist_mean, dist_std, dist_lower_threshold, dist_upper_threshold))
    for i in range(necoregions):
        for j in range(i, necoregions):
            if (ts_dist_matrix[i, j] < dist_lower_threshold) and (j != i):
                message('-- Similarity: ecoregions %s and %s are similar \
                        R = %.3f Rp = %.3f Tp = %.3f Lp = %.3f D = %.3f' %
                        (ecoregion_IDs[i], ecoregion_IDs[j],
                         ts_corr_matrix[i, j], ts_rsig_matrix[i, j],
                         ts_tsig_matrix[i, j], ts_lsig_matrix[i, j],
                         ts_dist_matrix[i, j]))
    for i in range(necoregions):
        for j in range(i, necoregions):
            if (ts_dist_matrix[i, j] > dist_upper_threshold) and (j != i):
                message('-- Dissimilarity: ecoregions %s and %s are \
                        dissimilar R = %.3f Rp = %.3f Tp = %.3f Lp = %.3f \
                        D = %.3f' %
                        (ecoregion_IDs[i], ecoregion_IDs[j],
                         ts_corr_matrix[i, j], ts_rsig_matrix[i, j],
                         ts_tsig_matrix[i, j], ts_lsig_matrix[i, j],
                         ts_dist_matrix[i, j]))
    message('- saving time series cross-testing matrices to %s' % tsfname)
    with hdf.File(tsfname, 'r+') as h5outfile:
        vname = '%s_corr_matrix' % varnames[k]
        if vname in h5outfile.keys():
            del h5outfile[vname]
        h5outfile.create_dataset(vname, data=ts_corr_matrix,
                                 dtype=np.float32, compression='gzip')
        message('-- %s' % vname)
        vname = '%s_rsig_matrix' % varnames[k]
        if vname in h5outfile.keys():
            del h5outfile[vname]
        h5outfile.create_dataset(vname, data=ts_rsig_matrix,
                                 dtype=np.float32, compression='gzip')
        message('-- %s' % vname)
        vname = '%s_ttest_matrix' % varnames[k]
        if vname in h5outfile.keys():
            del h5outfile[vname]
        h5outfile.create_dataset(vname, data=ts_ttest_matrix,
                                 dtype=np.float32, compression='gzip')
        message('-- %s' % vname)
        vname = '%s_tsig_matrix' % varnames[k]
        if vname in h5outfile.keys():
            del h5outfile[vname]
        h5outfile.create_dataset(vname, data=ts_tsig_matrix,
                                 dtype=np.float32, compression='gzip')
        message('-- %s' % vname)
        vname = '%s_levene_matrix' % varnames[k]
        if vname in h5outfile.keys():
            del h5outfile[vname]
        h5outfile.create_dataset(vname, data=ts_levene_matrix,
                                 dtype=np.float32, compression='gzip')
        message('-- %s' % vname)
        vname = '%s_lsig_matrix' % varnames[k]
        if vname in h5outfile.keys():
            del h5outfile[vname]
        h5outfile.create_dataset(vname, data=ts_lsig_matrix,
                                 dtype=np.float32, compression='gzip')
        message('-- %s' % vname)
        vname = '%s_dist_matrix' % varnames[k]
        if vname in h5outfile.keys():
            del h5outfile[vname]
        h5outfile.create_dataset(vname, data=ts_dist_matrix,
                                 dtype=np.float32, compression='gzip')
        message('-- %s' % vname)
    message(' ')
#
message('process_NCEI_11.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_11.py
