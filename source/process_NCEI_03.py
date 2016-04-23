"""
Python script 'process_NCEI_03.py'
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

USAGE: 'python process_NCEI_03.py NCEI_WLS_19840101-20131231 ./grids'

PURPOSE: Temporal accumulation of various gridded daily meteorological data values.
         FD: freezing days for TMIN, TAVG, and TMAX (counted from 1 Jul)
         CD: chilling days using TAVG relative to 5 dC (counted from 1 Jul)
         CDD: chilling degree days using TAVG relative to 5 dC (accumulated from 1 Jul)
         GDD: growing degree days using TAVG relative to 5 dC (accumulated from 1 Jan)
         GDD_base0: growing degree days using TAVG relative to 0 dC (accumulated from 1 Jan)
         Tx_90d: 90-day mean and variance of TMIN, TAVG, and TMAX
         P_Xd: 30-day, 90-day, 180-day, and 365-day precipitation totals
         P_90d_X: 90-day precipitation days at thresholds of 0, 10, and 25 mm

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O
			  The pp module is ParallelPython, employed here for multi-processor application

NOTE 1: The imported module numpy is NOT aliased (e.g. "numpy as np") so that there is no 
        confusion in ParallelPython calls to embedded subroutines.

NOTE 2: This script uses 1 ParallelPython server with a number of CPUs specified below. 
        Make sure the total falls within your number of available processors.

INPUT: Output file from process_NCEI_01.py script in '.h5' format (one file)
       (with the naming convention 'data/NCEI_WLS_[YYYYMMDD]-[YYYYMMDD]_processed.h5')
       Output files from process_NCEI_02.py script in '.h5' format (1 / day)
       (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_1.h5')

OUTPUT: Copied input '.h5' file with new accumulation grids (1 / day)
        (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_2.h5')

RUN TIME: ~85 mins per data year

NOTE: If you have enough processors, you can start this script pretty soon after starting
      the process_NCEI_02.py script and run them in parallel. The 02 script is faster, so 
      03 will not catch up to 02's production of '*_grids_1.h5' files (unless there is an 
      error that causes 02 to stop...). 

TO DO: make a serial version that does not depend on pp (pp availability or user's choice)
"""

import os, sys, datetime
import h5py as hdf
import numpy
import pp

def message(str):
    print str
    sys.stdout.flush()
    return    

def cube_sum(varname,nd,var_cube,var_grid,stns_all,stns):
    var_cube[0:nd-1,:,:] = var_cube[1:nd,:,:]
    var_cube[nd-1,:,:] = var_grid[:,:]
    grid_var_sum = numpy.sum(var_cube, axis=0)
    stns_all_dims = numpy.shape(stns_all)
    if stns_all_dims[0] == nd:
        stns_all[0:nd-1][:] = stns_all[1:nd][:]
        stns_all[nd-1] = stns
    else:
        stns_all.append(stns)
    stns_all_list = []
    stns_all_dims = numpy.shape(stns_all)
    for i in range(0,stns_all_dims[0]):
        stns_all_list = list(set(stns_all_list) | set(stns_all[i]))
    return varname,grid_var_sum,stns_all_list,var_cube,stns_all

def cube_sum_parts(varname,nd,var_cube1,var_cube2,var_grid,stns_all,stns):
    var_cube1_dims = numpy.shape(var_cube1)
    nd1 = var_cube1_dims[0]
    var_cube2_dims = numpy.shape(var_cube2)
    nd2 = var_cube2_dims[0]
    var_cube1[0:nd1-1,:,:] = var_cube1[1:nd1,:,:]
    var_cube1[nd1-1,:,:] = var_cube2[0,:,:]
    var_cube2[0:nd2-1,:,:] = var_cube2[1:nd2,:,:]
    var_cube2[nd2-1,:,:] = var_grid[:,:]
    grid_var_sum = numpy.sum(var_cube1, axis=0) + numpy.sum(var_cube2, axis=0)
    stns_all_dims = numpy.shape(stns_all)
    if stns_all_dims[0] == nd:
        stns_all[0:nd-1][:] = stns_all[1:nd][:]
        stns_all[nd-1] = stns
    else:
        stns_all.append(stns)
    stns_all_list = []
    stns_all_dims = numpy.shape(stns_all)
    for i in range(0,stns_all_dims[0]):
        stns_all_list = list(set(stns_all_list) | set(stns_all[i]))
    return varname,grid_var_sum,stns_all_list,var_cube1,var_cube2,stns_all

def cube_mean_var(varname,nd,var_cube,var_grid,stns_all,stns):
    var_cube[0:nd-1,:,:] = var_cube[1:nd,:,:]
    var_cube[nd-1,:,:] = var_grid[:,:]
    grid_var_mean = numpy.mean(var_cube, axis=0)
    grid_var_var = numpy.var(var_cube, axis=0)
    stns_all_dims = numpy.shape(stns_all)
    if stns_all_dims[0] == nd:
        stns_all[0:nd-1][:] = stns_all[1:nd][:]
        stns_all[nd-1] = stns
    else:
        stns_all.append(stns)
    stns_all_list = []
    stns_all_dims = numpy.shape(stns_all)
    for i in range(0,stns_all_dims[0]):
        stns_all_list = list(set(stns_all_list) | set(stns_all[i]))
    return varname,grid_var_mean,grid_var_var,stns_all_list,var_cube,stns_all

def cube_threshold_count(varname,nd,var_cube,var_grid,thresh):
    thresh_var = numpy.where(var_grid > thresh, 1.0, 0.0)
    var_cube[0:nd-1,:,:] = var_cube[1:nd,:,:]
    var_cube[nd-1,:,:] = thresh_var[:,:]
    grid_var_count = numpy.sum(var_cube, axis=0)
    return varname,grid_var_count,var_cube

def grid_threshold_count(varname,mm,dd,reset_mm,reset_dd,var_grid,grid_var_prev,stns_all=0,stns=0):
    grid_var_new = numpy.where(var_grid < 0.0, 1.0, 0.0)
    if mm == reset_mm and dd == reset_dd:
        grid_var = grid_var_new
        if stns_all != 0:
            stns_all = stns
    else:
        grid_var = grid_var_prev + grid_var_new
        if stns_all != 0:
            stns_all = list(set(stns_all) | set(stns))
    if stns_all == 0:
        return varname,grid_var
    else:
        return varname,grid_var,stns_all

def grid_threshold_accumulate(varname,mm,dd,reset_mm,reset_dd,var_grid,grid_var_prev,stns_all,stns):
    grid_var_new = numpy.where(var_grid > 0.0, var_grid, 0.0)
    if mm == reset_mm and dd == reset_dd:
        grid_var = grid_var_new
        stns_all = stns
    else:
        grid_var = grid_var_prev + grid_var_new
        stns_all = list(set(stns_all) | set(stns))
    return varname,grid_var,stns_all

def write_to_file(file,gvar,gdata,svar,sdata):
    if gvar in file.keys():
        del file[gvar]
    file.create_dataset(gvar, data=gdata, dtype=numpy.float32, compression='gzip')
    stnsdir = 'stns'
    datapath = '%s/%s' % (stnsdir, svar)
    if svar in file[stnsdir].keys():
        del file[datapath]
    file.create_dataset(datapath, data=sdata)
    message('- %s %s with %d stations' % (gvar, str(gdata.shape), len(sdata)))
    return

message(' ')
message('process_NCEI_03.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
ncpus = 8
tbase = 5.0
cd_mm_start = 7
cd_dd_start = 1
cd_start_str = '1 Jul'
gd_mm_start = 1
gd_dd_start = 1
gd_start_str = '1 Jan'
#
if len(sys.argv) < 3:
    message('input warning: no input/output grids directory path indicated, using ./grids')
    path = './grids'
else:
    path = sys.argv[2]
#
if len(sys.argv) < 2:
    message('input error: need prefix for h5 file containing NCEI weather data')
    sys.exit(1)
else:
    NCEIfname = sys.argv[1]
h5infname = '%s/../data/%s_processed.h5' % (path, NCEIfname)
#
message('reading dates information from %s' % h5infname)
h5infile = hdf.File(h5infname,'r')
dates = numpy.copy(h5infile['dates'])
h5infile.close()
message('%d dates to process' % len(dates))
message(' ')
#
h5infname = '%s/%d_NCEI_grids_1.h5' % (path, dates[0])
message('extracting grid information from %s' % h5infname)
h5infile = hdf.File(h5infname,'r')
nrows = numpy.copy(h5infile['grid/nrows'])
ncols = numpy.copy(h5infile['grid/ncols'])
h5infile.close()
message('establishing grids and datacubes for rolling accounting')
tmin_frz_prev = numpy.zeros((nrows,ncols))
tavg_frz_prev = numpy.zeros((nrows,ncols))
tmax_frz_prev = numpy.zeros((nrows,ncols))
chill_d_prev = numpy.zeros((nrows,ncols))
chill_dd_prev = numpy.zeros((nrows,ncols))
grow_dd_prev = numpy.zeros((nrows,ncols))
grow_dd_base0_prev = numpy.zeros((nrows,ncols))
chill_d_stns = []
chill_dd_stns = []
grow_dd_stns = []
grow_dd_base0_stns = []
tmin_90d = numpy.zeros((90,nrows,ncols))
tavg_90d = numpy.zeros((90,nrows,ncols))
tmax_90d = numpy.zeros((90,nrows,ncols))
tmin_90d_stns = []
tavg_90d_stns = []
tmax_90d_stns = [] 
prcp_30d = numpy.zeros((30,nrows,ncols))
prcp_90d = numpy.zeros((90,nrows,ncols))
prcp_90d_nd0 = numpy.zeros((90,nrows,ncols))
prcp_90d_nd10 = numpy.zeros((90,nrows,ncols))
prcp_90d_nd25 = numpy.zeros((90,nrows,ncols))
prcp_180d = numpy.zeros((180,nrows,ncols))
prcp_365d_p1 = numpy.zeros((183,nrows,ncols))
prcp_365d_p2 = numpy.zeros((182,nrows,ncols))
prcp_30d_stns = []
prcp_90d_stns = []
prcp_180d_stns = []
prcp_365d_stns = [] 
#
ppservers0 = ()
job_server0 = pp.Server(ncpus, ppservers=ppservers0)
message('starting ParallelPython (pp) server #0 at %s' % datetime.datetime.now().isoformat())
message('- using %d processors' % job_server0.get_ncpus())
message(' ')
dc_sum = pp.Template(job_server0, cube_sum, depfuncs=(), modules=("numpy",))
dc_sum_parts = pp.Template(job_server0, cube_sum_parts, depfuncs=(), modules=("numpy",))
dc_mean = pp.Template(job_server0, cube_mean_var, depfuncs=(), modules=("numpy",))
dc_thresh_count = pp.Template(job_server0, cube_threshold_count, depfuncs=(), modules=("numpy",))
dg_thresh_count = pp.Template(job_server0, grid_threshold_count, depfuncs=(), modules=("numpy",))
dg_thresh_accum = pp.Template(job_server0, grid_threshold_accumulate, depfuncs=(), modules=("numpy",))
#
for date in dates:
    h5infname = '%s/%d_NCEI_grids_1.h5' % (path, date)
    message('extracting PRCP, TMAX, TMIN, and TAVG grids from %s' % h5infname)
    h5infile = hdf.File(h5infname,'r')
    prcp_stns = numpy.copy(h5infile['stns/prcp_stns'])
    prcp = numpy.copy(h5infile['grid_prcp'])
    tmax_stns = numpy.copy(h5infile['stns/tmax_stns'])
    tmax = numpy.copy(h5infile['grid_tmax'])
    tmin_stns = numpy.copy(h5infile['stns/tmin_stns'])
    tmin = numpy.copy(h5infile['grid_tmin'])
    tavg_stns = numpy.copy(h5infile['stns/tavg_stns'])
    tavg = numpy.copy(h5infile['grid_tavg'])
    h5infile.close()
    #
    year = date // 10000
    month = (date - (year * 10000)) // 100
    day = date - (year * 10000) - (month * 100)
    #
    message('- calculating grid CD, CDD, and GDD (base 5dC and base 0dC)')
    chill_d = tavg - tbase
    chill_dd = tbase - tavg
    grow_dd = tavg - tbase
    grow_dd_base0 = tavg
    message('  mean CD %.1f    mean CDD %.1f    mean GDD %.1f    mean GDD_base0 %.1f' % (numpy.mean(chill_d),numpy.mean(chill_dd),numpy.mean(grow_dd),numpy.mean(grow_dd_base0)))
    #
    jobs_dc_sum = []
    jobs_dc_sum_parts = []
    jobs_dc_mean = []
    jobs_dc_thresh_count = []
    jobs_dg_thresh_count = []
    jobs_dg_thresh_accum = []
    #
    jobs_dg_thresh_count.append(dg_thresh_count.submit('grid_tmin_frz',month,day,cd_mm_start,cd_dd_start,tmin,tmin_frz_prev))
    jobs_dg_thresh_count.append(dg_thresh_count.submit('grid_tavg_frz',month,day,cd_mm_start,cd_dd_start,tavg,tavg_frz_prev))
    jobs_dg_thresh_count.append(dg_thresh_count.submit('grid_tmax_frz',month,day,cd_mm_start,cd_dd_start,tmax,tmax_frz_prev))
    jobs_dg_thresh_count.append(dg_thresh_count.submit('grid_chill_d',month,day,cd_mm_start,cd_dd_start,chill_d,chill_d_prev,chill_d_stns,tavg_stns))
    jobs_dg_thresh_accum.append(dg_thresh_accum.submit('grid_chill_dd',month,day,cd_mm_start,cd_dd_start,chill_dd,chill_dd_prev,chill_dd_stns,tavg_stns))
    jobs_dg_thresh_accum.append(dg_thresh_accum.submit('grid_grow_dd',month,day,gd_mm_start,gd_dd_start,grow_dd,grow_dd_prev,grow_dd_stns,tavg_stns))
    jobs_dg_thresh_accum.append(dg_thresh_accum.submit('grid_grow_dd_base0',month,day,gd_mm_start,gd_dd_start,grow_dd_base0,grow_dd_base0_prev,grow_dd_base0_stns,tavg_stns))
    jobs_dc_mean.append(dc_mean.submit('grid_tmin_90d',90,tmin_90d,tmin,tmin_90d_stns,tmin_stns))
    jobs_dc_mean.append(dc_mean.submit('grid_tavg_90d',90,tavg_90d,tavg,tavg_90d_stns,tavg_stns))
    jobs_dc_mean.append(dc_mean.submit('grid_tmax_90d',90,tmax_90d,tmax,tmax_90d_stns,tmax_stns))
    jobs_dc_sum.append(dc_sum.submit('grid_prcp_30d',30,prcp_30d,prcp,prcp_30d_stns,prcp_stns))
    jobs_dc_sum.append(dc_sum.submit('grid_prcp_90d',90,prcp_90d,prcp,prcp_90d_stns,prcp_stns))
    jobs_dc_sum.append(dc_sum.submit('grid_prcp_180d',180,prcp_180d,prcp,prcp_180d_stns,prcp_stns))
    jobs_dc_sum_parts.append(dc_sum_parts.submit('grid_prcp_365d',365,prcp_365d_p1,prcp_365d_p2,prcp,prcp_365d_stns,prcp_stns))
    jobs_dc_thresh_count.append(dc_thresh_count.submit('grid_prcp_90d_nd0',90,prcp_90d_nd0,prcp,0.0))
    jobs_dc_thresh_count.append(dc_thresh_count.submit('grid_prcp_90d_nd10',90,prcp_90d_nd10,prcp,1.0))
    jobs_dc_thresh_count.append(dc_thresh_count.submit('grid_prcp_90d_nd25',90,prcp_90d_nd25,prcp,2.5))
    #
    for job in jobs_dg_thresh_count:
        returns = job()
        varname = returns[0]
        if varname == 'grid_tmin_frz':
            grid_tmin_frz = returns[1]
            tmin_frz_prev = grid_tmin_frz
            message('- calculated Tmin freezing days (accumulated from %s) mean %.1f' % (cd_start_str, numpy.mean(grid_tmin_frz)))
        elif varname == 'grid_tavg_frz':
            grid_tavg_frz = returns[1]
            tavg_frz_prev = grid_tavg_frz
            message('- calculated Tavg freezing days (accumulated from %s) mean %.1f' % (cd_start_str, numpy.mean(grid_tavg_frz)))
        elif varname == 'grid_tmax_frz':
            grid_tmax_frz = returns[1]
            tmax_frz_prev = grid_tmax_frz
            message('- calculated Tmax freezing days (accumulated from %s) mean %.1f' % (cd_start_str, numpy.mean(grid_tmax_frz)))
        elif varname == 'grid_chill_d':
            grid_chill_d = returns[1]
            chill_d_stns = returns[2]
            chill_d_prev = grid_chill_d
            message('- calculated chilling days (based on Tavg, accumulated from %s) mean %.1f' % (cd_start_str, numpy.mean(grid_chill_d)))
    # 
    for job in jobs_dg_thresh_accum:
        returns = job()
        varname = returns[0]
        if varname == 'grid_chill_dd':
            grid_chill_dd = returns[1]
            chill_dd_stns = returns[2]
            chill_dd_prev = grid_chill_dd
            message('- calculated chilling degree days (based on Tavg, accumulated from %s) mean %.1f' % (cd_start_str, numpy.mean(grid_chill_dd)))
        elif varname == 'grid_grow_dd':
            grid_grow_dd = returns[1]
            grow_dd_stns = returns[2]
            grow_dd_prev = grid_grow_dd
            message('- calculated growing degree days with base 5dC (based on Tavg, accumulated from %s) mean %.1f' % (gd_start_str, numpy.mean(grid_grow_dd)))
        elif varname == 'grid_grow_dd_base0':
            grid_grow_dd_base0 = returns[1]
            grow_dd_base0_stns = returns[2]
            grow_dd_base0_prev = grid_grow_dd_base0
            message('- calculated growing degree days with base 0dC (based on Tavg, accumulated from %s) mean %.1f' % (gd_start_str, numpy.mean(grid_grow_dd_base0)))
    # 
    for job in jobs_dc_mean:
        returns = job()
        varname = returns[0]
        if varname == 'grid_tmin_90d':
            grid_tmin_90d_mean = returns[1]
            grid_tmin_90d_var = returns[2]
            tmin_90d_stns_all = returns[3]
            tmin_90d = returns[4]
            tmin_90d_stns = returns[5]
            message('- calculated updated 90-day running Tmin mean and variance, mean %.1f' % numpy.mean(grid_tmin_90d_mean))
        elif varname == 'grid_tavg_90d':
            grid_tavg_90d_mean = returns[1]
            grid_tavg_90d_var = returns[2]
            tavg_90d_stns_all = returns[3]
            tavg_90d = returns[4]
            tavg_90d_stns = returns[5]
            message('- calculated updated 90-day running Tavg mean and variance, mean %.1f' % numpy.mean(grid_tavg_90d_mean))
        elif varname == 'grid_tmax_90d':
            grid_tmax_90d_mean = returns[1]
            grid_tmax_90d_var = returns[2]
            tmax_90d_stns_all = returns[3]
            tmax_90d = returns[4]
            tmax_90d_stns = returns[5]
            message('- calculated updated 90-day running Tmax mean and variance, mean %.1f' % numpy.mean(grid_tmax_90d_mean))
    #
    for job in jobs_dc_sum:
        returns = job()
        varname = returns[0]
        if varname == 'grid_prcp_30d':
            grid_prcp_30d = returns[1]
            prcp_30d_stns_all = returns[2]
            prcp_30d = returns[3]
            prcp_30d_stns = returns[4]
            message('- calculated updated 30-day running precipitation total, mean %.1f' % numpy.mean(grid_prcp_30d))
        elif varname == 'grid_prcp_90d':
            grid_prcp_90d = returns[1]
            prcp_90d_stns_all = returns[2]
            prcp_90d = returns[3]
            prcp_90d_stns = returns[4]
            message('- calculated updated 90-day running precipitation total, mean %.1f' % numpy.mean(grid_prcp_90d))
        elif varname == 'grid_prcp_180d':
            grid_prcp_180d = returns[1]
            prcp_180d_stns_all = returns[2]
            prcp_180d = returns[3]
            prcp_180d_stns = returns[4]
            message('- calculated updated 180-day running precipitation total, mean %.1f' % numpy.mean(grid_prcp_180d))
    #
    for job in jobs_dc_sum_parts:
        returns = job()
        varname = returns[0]
        if varname == 'grid_prcp_365d':
            grid_prcp_365d = returns[1]
            prcp_365d_stns_all = returns[2]
            prcp_365d_p1 = returns[3]
            prcp_365d_p2 = returns[4]
            prcp_365d_stns = returns[5]
            message('- calculated updated 365-day running precipitation total, mean %.1f' % numpy.mean(grid_prcp_365d))
    # 
    for job in jobs_dc_thresh_count:
        returns = job()
        varname = returns[0]
        if varname == 'grid_prcp_90d_nd0':
            grid_prcp_90d_nd0 = returns[1]
            prcp_90d_nd0 = returns[2]
            message('- calculated updated 90-day running precipitation ndays total for P > 0 mm, mean %.1f' % numpy.mean(grid_prcp_90d_nd0))
        elif varname == 'grid_prcp_90d_nd10':
            grid_prcp_90d_nd10 = returns[1]
            prcp_90d_nd10 = returns[2]
            message('- calculated updated 90-day running precipitation ndays total for P > 10 mm, mean %.1f' % numpy.mean(grid_prcp_90d_nd10))
        elif varname == 'grid_prcp_90d_nd25':
            grid_prcp_90d_nd25 = returns[1]
            prcp_90d_nd25 = returns[2]
            message('- calculated updated 90-day running precipitation ndays total for P > 25 mm, mean %.1f' % numpy.mean(grid_prcp_90d_nd25))
    #
    h5outfname = '%s/%d_NCEI_grids_2.h5' % (path, date)
    cmdstring = 'cp %s %s' % (h5infname, h5outfname)
    os.system(cmdstring)
    message('saving grids to %s' % h5outfname)
    h5outfile = hdf.File(h5outfname,'r+')
    del h5outfile['meta/last_updated']
    h5outfile.create_dataset('meta/last_updated', data=datetime.datetime.now().isoformat())
    del h5outfile['meta/at']
    outstr = '90d temp mean + var, FD, CD, CDD, GDD, GDD_base0, precip running totals, 90d precip day counts'
    h5outfile.create_dataset('meta/at', data=outstr)
    #
    write_to_file(h5outfile,'tmin_frz_days',grid_tmin_frz,'tmin_frz_stns',tmin_stns)
    write_to_file(h5outfile,'tavg_frz_days',grid_tavg_frz,'tavg_frz_stns',tavg_stns)
    write_to_file(h5outfile,'tmax_frz_days',grid_tmax_frz,'tmax_frz_stns',tmax_stns)
    write_to_file(h5outfile,'chill_d',grid_chill_d,'chill_d_stns',chill_d_stns)
    write_to_file(h5outfile,'chill_dd',grid_chill_dd,'chill_dd_stns',chill_dd_stns)
    write_to_file(h5outfile,'grow_dd',grid_grow_dd,'grow_dd_stns',grow_dd_stns)
    write_to_file(h5outfile,'grow_dd_base0',grid_grow_dd_base0,'grow_dd_base0_stns',grow_dd_base0_stns)
    write_to_file(h5outfile,'tmin_90d_avg',grid_tmin_90d_mean,'tmin_90d_stns',tmin_90d_stns_all)
    write_to_file(h5outfile,'tmin_90d_var',grid_tmin_90d_var,'tmin_90d_stns',tmin_90d_stns_all)
    write_to_file(h5outfile,'tavg_90d_avg',grid_tavg_90d_mean,'tavg_90d_stns',tavg_90d_stns_all)
    write_to_file(h5outfile,'tavg_90d_var',grid_tavg_90d_var,'tavg_90d_stns',tavg_90d_stns_all)
    write_to_file(h5outfile,'tmax_90d_avg',grid_tmax_90d_mean,'tmax_90d_stns',tmax_90d_stns_all)
    write_to_file(h5outfile,'tmax_90d_var',grid_tmax_90d_var,'tmax_90d_stns',tmax_90d_stns_all)
    write_to_file(h5outfile,'prcp_30d_sum',grid_prcp_30d,'prcp_30d_stns',prcp_30d_stns_all)
    write_to_file(h5outfile,'prcp_90d_sum',grid_prcp_90d,'prcp_90d_stns',prcp_90d_stns_all)
    write_to_file(h5outfile,'prcp_180d_sum',grid_prcp_180d,'prcp_180d_stns',prcp_180d_stns_all)
    write_to_file(h5outfile,'prcp_365d_sum',grid_prcp_365d,'prcp_365d_stns',prcp_365d_stns_all)
    write_to_file(h5outfile,'prcp_90d_nd0_sum',grid_prcp_90d_nd0,'prcp_90d_nd0_stns',prcp_90d_stns_all)
    write_to_file(h5outfile,'prcp_90d_nd10_sum',grid_prcp_90d_nd10,'prcp_90d_nd10_stns',prcp_90d_stns_all)
    write_to_file(h5outfile,'prcp_90d_nd25_sum',grid_prcp_90d_nd25,'prcp_90d_nd25_stns',prcp_90d_stns_all)
    h5outfile.close()
    message(' ')
#
job_server0.destroy()
message('server #0 parallel processing completed at %s' % datetime.datetime.now().isoformat())
message(' ')
#
message('process_NCEI_03.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_03.py
