"""
Python script 'process_NCEI_15.py'
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

USAGE: 'python process_NCEI_15.py 1984 2013 ./analyses'

PURPOSE: Cross-correlation of key dates + statistics from variable time series by ecoregion 
         cluster, along with several climatological (teleconnection) indices.

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O
              The 'Stats' and 'Teleconnections' modules have their own requirements

INPUT: Full-grid time series in '.h5' files from 'process_NCEI_10.py'
       Ecoregion cluster time series and stats in '.h5' files from 'process_NCEI_13.py'
       Climatological teleconnection time series from various sources in '.csv' files

OUTPUT: Ecoregion cluster-based time series correlation statistics in '.h5' files

RUN TIME: ~5 sec
"""

import os, sys, datetime
import h5py as hdf
import numpy as np
from Stats import *
from Teleconnections import *

def message(str):
    print str
    sys.stdout.flush()
    return    

def make_correl_code(idx1, idx2, correl):
    if correl > 0:
        code = '%s+%s' % (str(idx1).zfill(3), str(idx2).zfill(3))
    else:
        code = '%s-%s' % (str(idx1).zfill(3), str(idx2).zfill(3))
    return code

message(' ')
message('process_NCEI_15.py started at %s' % datetime.datetime.now().isoformat())
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
clustertsfname = '%s/%d-%d_cluster_timeseries.h5' % (path,year_begin,year_end)
fullgridtsfname = '%s/%d-%d_full_grid_timeseries.h5' % (path,year_begin,year_end)
years = np.arange(year_begin,year_end+1).astype(int)
nyears = len(years)
#
sig_threshold_1 = 0.001
sig_threshold_2 = 0.01
sig_threshold_3 = 0.05
#
varnames = []
varnames.append(('Tavg 90d avg at veq','cluster_series[i+1,:]'))
varnames.append(('Tavg 90d var at veq','cluster_series[i+1,:]'))
varnames.append(('Prcp 90d total at veq','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 0 in 90d prior to veq','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 10 mm in 90d prior to veq','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 25 mm in 90d prior to veq','cluster_series[i+1,:]'))
varnames.append(('Tavg 90d avg at ssol','cluster_series[i+1,:]'))
varnames.append(('Tavg 90d var at ssol','cluster_series[i+1,:]'))
varnames.append(('Prcp 90d total at ssol','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 0 in 90d prior to ssol','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 10 mm in 90d prior to ssol','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 25 mm in 90d prior to ssol','cluster_series[i+1,:]'))
varnames.append(('Tavg 90d avg at aeq','cluster_series[i+1,:]'))
varnames.append(('Tavg 90d var at aeq','cluster_series[i+1,:]'))
varnames.append(('Prcp 90d total at aeq','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 0 in 90d prior to aeq','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 10 mm in 90d prior to aeq','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 25 mm in 90d prior to aeq','cluster_series[i+1,:]'))
varnames.append(('Tavg 90d avg at wsol','cluster_series[i+1,:]'))
varnames.append(('Tavg 90d var at wsol','cluster_series[i+1,:]'))
varnames.append(('Prcp 90d total at wsol','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 0 in 90d prior to wsol','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 10 mm in 90d prior to wsol','cluster_series[i+1,:]'))
varnames.append(('Prcp days > 25 mm in 90d prior to wsol','cluster_series[i+1,:]'))
varnames.append(('CD at beginning of plateau','cluster_series[i+1,:]'))
varnames.append(('Tmin FD at ssol','cluster_series[i+1,:]'))
varnames.append(('Cold season intensity','cluster_series[i+1,:]'))
varnames.append(('DOY at last spring Tmin freeze','cluster_series[i+1,:]'))
varnames.append(('GDD at last spring Tmin freeze','cluster_series[i+1,:]'))
varnames.append(('Prcp 365d total at eoy','cluster_series[i+1,:]'))
varnames.append(('DOY at beginning of plateau','cluster_series[i+1,:]'))
varnames.append(('GDD at beginning of plateau','cluster_series[i+1,:]'))
varnames.append(('GDD btw veq and aeq','cluster_series[i+1,:]'))
varnames.append(('DOY at first autumn Tmin freeze','cluster_series[i+1,:]'))
varnames.append(('DOY at end of plateau','cluster_series[i+1,:]'))
varnames.append(('GDD at end of plateau','cluster_series[i+1,:]'))
varnames.append(('Frost-free calendar days','cluster_series[i+1,:]'))
varnames.append(('Plateau calendar days','cluster_series[i+1,:]'))
varnames.append(('Plateau total GDD','cluster_series[i+1,:]'))
varnames.append(('Growing season intensity','cluster_series[i+1,:]'))
#
varnames.append(('Lake Superior ice-on','ice_on_doy[0:nyears]'))
varnames.append(('Lake Superior ice-off','ice_off_doy[0:nyears]'))
varnames.append(('Lake Superior ice duration','ice_duration[0:nyears]'))
varnames.append(('PDO index for prec DJF','pdo_winter[0:nyears]'))
varnames.append(('PDO index for conc MAM','pdo_spring[0:nyears]'))
varnames.append(('PDO index for conc JJA','pdo_summer[0:nyears]'))
varnames.append(('PDO index for conc SON','pdo_autumn[0:nyears]'))
varnames.append(('PDO index for succ DJF','pdo_winter[1:nyears+1]'))
varnames.append(('Nino3 index for prec DJFM','nino3[0:nyears]'))
varnames.append(('Nino4 index for prec DJFM','nino4[0:nyears]'))
varnames.append(('Nino3.4 index for prec DJFM','nino34[0:nyears]'))
varnames.append(('Nino3.4 index for succ DJFM','nino34[1:nyears+1]'))
varnames.append(('PNA index for prec DJF','pna_winter[0:nyears]'))
varnames.append(('PNA index for conc MAM','pna_spring[0:nyears]'))
varnames.append(('PNA index for conc JJA','pna_summer[0:nyears]'))
varnames.append(('PNA index for conc SON','pna_autumn[0:nyears]'))
varnames.append(('PNA index for succ DJF','pna_winter[1:nyears+1]'))
varnames.append(('AMO index for prec DJF','amo_winter[0:nyears]'))
varnames.append(('AMO index for conc MAM','amo_spring[0:nyears]'))
varnames.append(('AMO index for conc JJA','amo_summer[0:nyears]'))
varnames.append(('AMO index for conc SON','amo_autumn[0:nyears]'))
varnames.append(('AMO index for succ DJF','amo_winter[1:nyears+1]'))
varnames.append(('AO index for prec DJF','ao_winter[0:nyears]'))
varnames.append(('AO index for conc MAM','ao_spring[0:nyears]'))
varnames.append(('AO index for conc JJA','ao_summer[0:nyears]'))
varnames.append(('AO index for conc SON','ao_autumn[0:nyears]'))
varnames.append(('AO index for succ DJF','ao_winter[1:nyears+1]'))
varnames.append(('NAO index for prec DJF','nao_winter[0:nyears]'))
varnames.append(('NAO index for conc MAM','nao_spring[0:nyears]'))
varnames.append(('NAO index for conc JJA','nao_summer[0:nyears]'))
varnames.append(('NAO index for conc SON','nao_autumn[0:nyears]'))
varnames.append(('NAO index for succ DJF','nao_winter[1:nyears+1]'))
nvars = len(varnames)
#
varnames2 = []
varnames2.append(('Tavg 90d avg at veq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Tavg 90d var at veq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp 90d total at veq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 0 in 90d prior to veq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 10 mm in 90d prior to veq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 25 mm in 90d prior to veq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Tavg 90d avg at ssol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Tavg 90d var at ssol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp 90d total at ssol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 0 in 90d prior to ssol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 10 mm in 90d prior to ssol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 25 mm in 90d prior to ssol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Tavg 90d avg at aeq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Tavg 90d var at aeq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp 90d total at aeq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 0 in 90d prior to aeq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 10 mm in 90d prior to aeq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 25 mm in 90d prior to aeq','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Tavg 90d avg at wsol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Tavg 90d var at wsol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp 90d total at wsol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 0 in 90d prior to wsol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 10 mm in 90d prior to wsol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Prcp days > 25 mm in 90d prior to wsol','cluster_series[i+1,0:nyears-1]'))
varnames2.append(('Tavg 90d avg next veq','cluster_series[i-23,1:nyears]'))
varnames2.append(('Tavg 90d var next veq','cluster_series[i-23,1:nyears]'))
varnames2.append(('Prcp 90d total next veq','cluster_series[i-23,1:nyears]'))
varnames2.append(('Prcp days > 0 in 90d prior to next veq','cluster_series[i-23,1:nyears]'))
varnames2.append(('Prcp days > 10 mm in 90d prior to next veq','cluster_series[i-23,1:nyears]'))
varnames2.append(('Prcp days > 25 mm in 90d prior to next veq','cluster_series[i-23,1:nyears]'))
varnames2.append(('CD at beginning of next plateau','cluster_series[i-6,1:nyears]'))
varnames2.append(('Tmin FD at next ssol','cluster_series[i-6,1:nyears]'))
varnames2.append(('Next cold season intensity','cluster_series[i-6,1:nyears]'))
varnames2.append(('Next DOY at last spring Tmin freeze','cluster_series[i-6,1:nyears]'))
varnames2.append(('Next GDD at last spring Tmin freeze','cluster_series[i-6,1:nyears]'))
varnames2.append(('DOY at beginning of next plateau','cluster_series[i-5,1:nyears]'))
varnames2.append(('GDD at beginning of next plateau','cluster_series[i-5,1:nyears]'))
nvars2 = len(varnames2)
#
fname = '%s/../data/NSIDC_MIFL_Superior_Ice.csv' % path
message('retrieving climatology data from %s' % fname)
ice_on_doy, ice_off_doy, ice_duration = get_ice_dates(fname, year_begin, year_end)
mask = np.where(ice_duration != 0, 1, 0)
xvals = np.arange(0,nyears).astype(int)
slope, intercept, r_value, p_value, std_err = linregress(xvals*mask,ice_on_doy*mask)
message('- ice-on date regression: slope = %.2f  r = %.3f  p = %.4f' % (slope,r_value,p_value))
slope, intercept, r_value, p_value, std_err = linregress(xvals*mask,ice_off_doy*mask)
message('- ice-off date regression: slope = %.2f  r = %.3f  p = %.4f' % (slope,r_value,p_value))
slope, intercept, r_value, p_value, std_err = linregress(xvals*mask,ice_duration*mask)
message('- ice duration regression: slope = %.2f  r = %.3f  p = %.4f' % (slope,r_value,p_value))
#
fname = '%s/../data/NOAA_ESRL_PDO_indices.csv' % path
message('retrieving climatology data from %s' % fname)
pdo_winter, pdo_spring, pdo_summer, pdo_autumn = get_ao_idx(fname, year_begin, (year_end + 1))
#
fname = '%s/../data/NCEP_CPC_ENSO_indices.csv' % path
message('retrieving climatology data from %s' % fname)
nino3, nino4, nino34 = get_nino_idx(fname, year_begin, (year_end + 1))
#
fname = '%s/../data/NCEP_CPC_PNA_indices.csv' % path
message('retrieving climatology data from %s' % fname)
pna_winter, pna_spring, pna_summer, pna_autumn = get_ao_idx(fname, year_begin, (year_end + 1))
#
fname = '%s/../data/NOAA_ESRL_AMO_indices.csv' % path
message('retrieving climatology data from %s' % fname)
amo_winter, amo_spring, amo_summer, amo_autumn = get_ao_idx(fname, year_begin, (year_end + 1))
#
fname = '%s/../data/NCEP_CPC_AO_indices.csv' % path
message('retrieving climatology data from %s' % fname)
ao_winter, ao_spring, ao_summer, ao_autumn = get_ao_idx(fname, year_begin, (year_end + 1))
#
fname = '%s/../data/NCEP_CPC_NAO_indices.csv' % path
message('retrieving climatology data from %s' % fname)
nao_winter, nao_spring, nao_summer, nao_autumn = get_ao_idx(fname, year_begin, (year_end + 1))
message(' ')
#
message('getting number of clusters from %s' % clustertsfname)
h5file = hdf.File(clustertsfname,'r') 
nclusters = np.copy(h5file['nclusters'])
h5file.close()
message(' ')
#
for k in range(0,nclusters+1):
    if k < nclusters:
        tsfname = clustertsfname
        cluster_name = 'cluster_%d' % (k+1)
    else:
        tsfname = fullgridtsfname
        cluster_name = 'full_grid'
    message('getting %s variable time series values from %s' % (cluster_name,tsfname))
    h5file = hdf.File(tsfname,'r') 
    cluster_series = np.copy(h5file[cluster_name])
    h5file.close()
    message(' ')
    #
    message('%d indicator variables specified' % nvars)
    var_idx = np.arange(0,nvars).astype(int)
    var_stats = np.zeros((nvars,12))
    var_matrix = np.zeros((nvars,nyears))
    correl_matrix = np.zeros((nvars,nvars))
    sig_matrix = np.zeros((nvars,nvars))
    correl_codes = []
    message('- processing %s individual indicator statistics' % cluster_name)
    for i in range(0,nvars):
        message('-- %s %s' % (str(var_idx[i]).zfill(3), varnames[i][0]))
        var = eval(varnames[i][1])
        var_stats[i,:] = indicator_stats(var,nyears)
        var_matrix[i,:] = var
    message(' ')
    #
    message('- processing %s cross-indicator correlations' % cluster_name)
    for j in range(0,nvars):
        for i in range(j,nvars):
            if i != j:
                corr, sig = pearsonr(var_matrix[j,:], var_matrix[i,:])
                correl_matrix[j,i] = corr
                sig_matrix[j,i] = sig
            else:
                correl_matrix[j,i] = 1.0
    message(' ')
    #
    message('significant %s autocorrelations at p < %.2f for %d-%d' % (cluster_name, sig_threshold_3, year_begin, year_end))
    for i in range(0,nvars):
        if var_stats[i,6] < sig_threshold_3:
            message('  %s %s r2 = %.3f (p = %.4f)' % (str(var_idx[i]).zfill(3), varnames[i][0], var_stats[i,5]**2, var_stats[i,6]))
    message(' ')
    #
    message('significant %s linear regression trends at p < %.2f for %d-%d' % (cluster_name, sig_threshold_3, year_begin, year_end))
    for i in range(0,nvars):
        if var_stats[i,10] < sig_threshold_3:
            message('  %s %s slope = %.4f (r2 = %.3f, p = %.4f)' % (str(var_idx[i]).zfill(3), varnames[i][0], var_stats[i,7], var_stats[i,9]**2, var_stats[i,10]))
    message(' ')
    #
    message('significant %s cross-correlations at p < %.3f for %d-%d' % (cluster_name, sig_threshold_1, year_begin, year_end))
    for j in range(0,nvars-1):
        for i in range(j+1,nvars):
            if sig_matrix[j,i] < sig_threshold_1:
                code = make_correl_code(var_idx[j], var_idx[i], correl_matrix[j,i])
                correl_codes.append(code)
                message('  %s -- %s is correlated with %s (r = %.3f, p = %.4f)' % (code, varnames[j][0], varnames[i][0], correl_matrix[j,i], sig_matrix[j,i]))
    message(' ')
    #
    message('significant %s cross-correlations at p < %.2f for %d-%d' % (cluster_name, sig_threshold_2, year_begin, year_end))
    for j in range(0,nvars-1):
        for i in range(j+1,nvars):
            if (sig_matrix[j,i] >= sig_threshold_1) and (sig_matrix[j,i] < sig_threshold_2):
                code = make_correl_code(var_idx[j], var_idx[i], correl_matrix[j,i])
                correl_codes.append(code)
                message('  %s -- %s is correlated with %s (r = %.3f, p = %.4f)' % (code, varnames[j][0], varnames[i][0], correl_matrix[j,i], sig_matrix[j,i]))
    message(' ')
    #
    message('significant %s cross-correlations at p < %.2f for %d-%d' % (cluster_name, sig_threshold_3, year_begin, year_end))
    for j in range(0,nvars-1):
        for i in range(j+1,nvars):
            if (sig_matrix[j,i] >= sig_threshold_2) and (sig_matrix[j,i] < sig_threshold_3):
                code = make_correl_code(var_idx[j], var_idx[i], correl_matrix[j,i])
                correl_codes.append(code)
                message('  %s -- %s is correlated with %s (r = %.3f, p = %.4f)' % (code, varnames[j][0], varnames[i][0], correl_matrix[j,i], sig_matrix[j,i]))
    message(' ')
    #
    outfile = '%s/%d-%d_%s_correlations.h5' % (path,year_begin,year_end,cluster_name)
    print 'saving collected correlation stats to %s' % outfile
    h5file = hdf.File(outfile,'w')
    h5file.create_dataset('varnames', data=varnames)
    h5file.create_dataset('var_idx', data=var_idx, dtype=np.int32, compression='gzip')
    h5file.create_dataset('var_stats', data=var_stats, dtype=np.float32, compression='gzip')
    h5file.create_dataset('var_matrix', data=var_matrix, dtype=np.float32, compression='gzip')
    h5file.create_dataset('correl_matrix', data=correl_matrix, dtype=np.float32, compression='gzip')
    h5file.create_dataset('sig_matrix', data=sig_matrix, dtype=np.float32, compression='gzip')
    h5file.create_dataset('correl_codes', data=correl_codes)
    h5file.close()
    message(' ')
    #
    message('analyzing over-winter T-P correlations')
    message('%d indicator variables specified' % nvars2)
    var_idx = np.arange(0,nvars2).astype(int)
    var_matrix = np.zeros((nvars2,nyears-1))
    correl_matrix = np.zeros((nvars2,nvars2))
    sig_matrix = np.zeros((nvars2,nvars2))
    correl_codes = []
    for i in range(0,nvars2):
        message('-- %s %s' % (str(var_idx[i]).zfill(3), varnames2[i][0]))
        var = eval(varnames2[i][1])
        var_matrix[i,:] = var
    message(' ')
    message('- processing %s cross-indicator correlations' % cluster_name)
    for j in range(0,nvars2):
        for i in range(j,nvars2):
            if i != j:
                corr, sig = pearsonr(var_matrix[j,:], var_matrix[i,:])
                correl_matrix[j,i] = corr
                sig_matrix[j,i] = sig
            else:
                correl_matrix[j,i] = 1.0
    message(' ')
    #
    message('significant %s cross-correlations at p < %.3f for %d-%d' % (cluster_name, sig_threshold_1, year_begin, year_end))
    for j in range(0,nvars2-1):
        for i in range(j+1,nvars2):
            if sig_matrix[j,i] < sig_threshold_1:
                code = make_correl_code(var_idx[j], var_idx[i], correl_matrix[j,i])
                correl_codes.append(code)
                message('  %s -- %s is correlated with %s (r = %.3f, p = %.4f)' % (code, varnames2[j][0], varnames2[i][0], correl_matrix[j,i], sig_matrix[j,i]))
    message(' ')
    #
    message('significant %s cross-correlations at p < %.2f for %d-%d' % (cluster_name, sig_threshold_2, year_begin, year_end))
    for j in range(0,nvars2-1):
        for i in range(j+1,nvars2):
            if (sig_matrix[j,i] >= sig_threshold_1) and (sig_matrix[j,i] < sig_threshold_2):
                code = make_correl_code(var_idx[j], var_idx[i], correl_matrix[j,i])
                correl_codes.append(code)
                message('  %s -- %s is correlated with %s (r = %.3f, p = %.4f)' % (code, varnames2[j][0], varnames2[i][0], correl_matrix[j,i], sig_matrix[j,i]))
    message(' ')
    #
    message('significant %s cross-correlations at p < %.2f for %d-%d' % (cluster_name, sig_threshold_3, year_begin, year_end))
    for j in range(0,nvars2-1):
        for i in range(j+1,nvars2):
            if (sig_matrix[j,i] >= sig_threshold_2) and (sig_matrix[j,i] < sig_threshold_3):
                code = make_correl_code(var_idx[j], var_idx[i], correl_matrix[j,i])
                correl_codes.append(code)
                message('  %s -- %s is correlated with %s (r = %.3f, p = %.4f)' % (code, varnames2[j][0], varnames2[i][0], correl_matrix[j,i], sig_matrix[j,i]))
    message(' ')
    #
    message('saving collected overwinter correlation stats to %s' % outfile)
    h5file = hdf.File(outfile,'r+')
    h5file.create_dataset('varnames_overwinter', data=varnames2)
    h5file.create_dataset('var_idx_overwinter', data=var_idx, dtype=np.int32, compression='gzip')
    h5file.create_dataset('var_matrix_overwinter', data=var_matrix, dtype=np.float32, compression='gzip')
    h5file.create_dataset('correl_matrix_overwinter', data=correl_matrix, dtype=np.float32, compression='gzip')
    h5file.create_dataset('sig_matrix_overwinter', data=sig_matrix, dtype=np.float32, compression='gzip')
    h5file.create_dataset('correl_codes_overwinter', data=correl_codes)
    h5file.close()
    message(' ')
#
codes_by_cluster = []
all_codes = []
for k in range(0,nclusters):
    cluster_name = 'cluster_%d' % (k+1)
    infile = '%s/%d-%d_%s_correlations.h5' % (path,year_begin,year_end,cluster_name)
    print 'getting correlation codes from %s' % infile
    h5file = hdf.File(infile,'r')
    codes = np.copy(h5file['correl_codes'])
    print '- %d codes' % len(codes)
    codes_by_cluster.append(codes)
    all_codes = sorted(list(set(all_codes) | set(codes)))
message(' ')
print 'total %d codes in common across all clusters' % len(all_codes)
message(' ')
message('indicator codes')
var_idx = np.arange(0,nvars).astype(int)
for i in range(0,nvars):
    message('  %s %s' % (str(var_idx[i]).zfill(3), varnames[i][0]))
message(' ')
for code in all_codes:
    presence = []
    for k in range(0,nclusters):
        if code in codes_by_cluster[k]:
            presence.append('cluster_%d' % (k+1))
    message('%s in %d clusters: %s' % (code,len(presence),str(presence)))
message(' ')
#
message('process_NCEI_15.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_15.py
