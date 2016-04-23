"""
Python script 'process_NCEI_03_serial_checkpointed.py'
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

USAGE: 'python process_NCEI_03_serial_checkpointed.py NCEI_WLS_19830101-20151031 1983 ./grids'

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
              The pickle module is used here for year-begin and -end I/O
              The h5py module is required for handling of HDF5 files, used here for I/O

INPUT: Output file from process_NCEI_01.py script in '.h5' format (one file)
       (with the naming convention 'data/NCEI_WLS_[YYYYMMDD]-[YYYYMMDD]_processed.h5')
       Output files from process_NCEI_02.py script in '.h5' format (1 / day)
       (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_1.h5')

OUTPUT: Copied input '.h5' file with new accumulation grids (1 / day)
        (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_2.h5')
        Year-end '.h5' and '.pickle' files with rolling accounted variables

RUN TIME: ~X mins per data year

NOTE: You can start this script pretty soon after starting the process_NCEI_02.py script 
      and run them in parallel. The 02 script is faster, so 03 will not catch up to 02's 
      production of '*_grids_1.h5' files (unless there is an error that causes 02 to stop...) 
"""

import os, sys, datetime, glob, pickle
import h5py as hdf
import numpy as np

def message(str):
    print str
    sys.stdout.flush()
    return    

def cube_sum(nd,var_cube,var_grid,stns_all,stns):
    var_cube[0:nd-1,:,:] = var_cube[1:nd,:,:]
    var_cube[nd-1,:,:] = var_grid[:,:]
    grid_var_sum = np.sum(var_cube, axis=0)
    stns_all_dims = np.shape(stns_all)
    if stns_all_dims[0] == nd:
        stns_all[0:nd-1][:] = stns_all[1:nd][:]
        stns_all[nd-1] = stns
    else:
        stns_all.append(stns)
    stns_all_list = []
    stns_all_dims = np.shape(stns_all)
    for i in range(0,stns_all_dims[0]):
        stns_all_list = list(set(stns_all_list) | set(stns_all[i]))
    return grid_var_sum,stns_all_list,var_cube,stns_all

def cube_sum_parts(nd,var_cube1,var_cube2,var_grid,stns_all,stns):
    var_cube1_dims = np.shape(var_cube1)
    nd1 = var_cube1_dims[0]
    var_cube2_dims = np.shape(var_cube2)
    nd2 = var_cube2_dims[0]
    var_cube1[0:nd1-1,:,:] = var_cube1[1:nd1,:,:]
    var_cube1[nd1-1,:,:] = var_cube2[0,:,:]
    var_cube2[0:nd2-1,:,:] = var_cube2[1:nd2,:,:]
    var_cube2[nd2-1,:,:] = var_grid[:,:]
    grid_var_sum = np.sum(var_cube1, axis=0) + np.sum(var_cube2, axis=0)
    stns_all_dims = np.shape(stns_all)
    if stns_all_dims[0] == nd:
        stns_all[0:nd-1][:] = stns_all[1:nd][:]
        stns_all[nd-1] = stns
    else:
        stns_all.append(stns)
    stns_all_list = []
    stns_all_dims = np.shape(stns_all)
    for i in range(0,stns_all_dims[0]):
        stns_all_list = list(set(stns_all_list) | set(stns_all[i]))
    return grid_var_sum,stns_all_list,var_cube1,var_cube2,stns_all

def cube_mean_var(nd,var_cube,var_grid,stns_all,stns):
    var_cube[0:nd-1,:,:] = var_cube[1:nd,:,:]
    var_cube[nd-1,:,:] = var_grid[:,:]
    grid_var_mean = np.mean(var_cube, axis=0)
    grid_var_var = np.var(var_cube, axis=0)
    stns_all_dims = np.shape(stns_all)
    if stns_all_dims[0] == nd:
        stns_all[0:nd-1][:] = stns_all[1:nd][:]
        stns_all[nd-1] = stns
    else:
        stns_all.append(stns)
    stns_all_list = []
    stns_all_dims = np.shape(stns_all)
    for i in range(0,stns_all_dims[0]):
        stns_all_list = list(set(stns_all_list) | set(stns_all[i]))
    return grid_var_mean,grid_var_var,stns_all_list,var_cube,stns_all

def cube_threshold_count(nd,var_cube,var_grid,thresh):
    thresh_var = np.where(var_grid > thresh, 1.0, 0.0)
    var_cube[0:nd-1,:,:] = var_cube[1:nd,:,:]
    var_cube[nd-1,:,:] = thresh_var[:,:]
    grid_var_count = np.sum(var_cube, axis=0)
    return grid_var_count,var_cube

def grid_threshold_count(mm,dd,reset_mm,reset_dd,var_grid,grid_var_prev,stns_all=0,stns=0):
    grid_var_new = np.where(var_grid < 0.0, 1.0, 0.0)
    if mm == reset_mm and dd == reset_dd:
        grid_var = grid_var_new
        if stns_all != 0:
            stns_all = stns
    else:
        grid_var = grid_var_prev + grid_var_new
        if stns_all != 0:
            stns_all = list(set(stns_all) | set(stns))
    if stns_all == 0:
        return grid_var
    else:
        return grid_var,stns_all

def grid_threshold_accumulate(mm,dd,reset_mm,reset_dd,var_grid,grid_var_prev,stns_all,stns):
    grid_var_new = np.where(var_grid > 0.0, var_grid, 0.0)
    if mm == reset_mm and dd == reset_dd:
        grid_var = grid_var_new
        stns_all = stns
    else:
        grid_var = grid_var_prev + grid_var_new
        stns_all = list(set(stns_all) | set(stns))
    return grid_var,stns_all

def write_to_file(file,gvar,gdata,svar,sdata):
    if gvar in file.keys():
        del file[gvar]
    file.create_dataset(gvar, data=gdata, dtype=np.float32, compression='gzip')
    stnsdir = 'stns'
    datapath = '%s/%s' % (stnsdir, svar)
    if svar in file[stnsdir].keys():
        del file[datapath]
    file.create_dataset(datapath, data=sdata)
    message('- %s %s with %d stations' % (gvar, str(gdata.shape), len(sdata)))
    return

def write_stn_lists(path,year,var,contents):
    fname = '%s/%d_%s.pickle' (path, year, var)
    f = open(fname, 'wb')
    pickle.dump(contents, f)
    f.close()
    return

def get_stn_lists(path,year,var):
    fname = '%s/%d_%s.pickle' (path, year, var)
    f = open(fname, 'rb')
    contents = pickle.load(f)
    f.close()
    return contents

message(' ')
message('process_NCEI_03_serial_checkpointed.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
tbase = 5.0
cd_mm_start = 7
cd_dd_start = 1
cd_start_str = '1 Jul'
gd_mm_start = 1
gd_dd_start = 1
gd_start_str = '1 Jan'
#
if len(sys.argv) < 4:
    message('input warning: no input grids directory path indicated, using ./grids')
    path = './grids'
else:
    path = sys.argv[3]
#
if len(sys.argv) < 3:
    message('input error: need year to process')
    sys.exit(1)
else:
    this_year = int(sys.argv[2])
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
all_dates = np.copy(h5infile['dates'])
h5infile.close()
message('- information for %d total dates found' % len(all_dates))
dates = sorted([j for j in all_dates if int(j//1E4) == this_year])
message('- processing %d dates in %d' % (len(dates), this_year))
message(' ')
#
prev_year = this_year - 1
vars_files = sorted(glob.glob('%s/*_year_end_vars.h5' % path))
use_vars_file = False
if len(vars_files) > 0:
    for vars_file in vars_files:
        if str(prev_year) in vars_file:
            use_vars_file = True
            varfname = vars_file
            break
if use_vars_file:  # if rolling accounting variable files exist to be carried over from previous year
    message('extracting rolling accounting variable grids and datacubes from %s' % varfname)
    h5infile = hdf.File(varfname,'r')
    nrows = np.copy(h5infile['nrows'])
    ncols = np.copy(h5infile['ncols'])
    tmin_frz_prev = np.copy(h5infile['tmin_frz_prev'])
    tavg_frz_prev = np.copy(h5infile['tavg_frz_prev'])
    tmax_frz_prev = np.copy(h5infile['tmax_frz_prev'])
    chill_d_prev = np.copy(h5infile['chill_d_prev'])
    chill_dd_prev = np.copy(h5infile['chill_dd_prev'])
    grow_dd_prev = np.copy(h5infile['grow_dd_prev'])
    grow_dd_base0_prev = np.copy(h5infile['grow_dd_base0_prev'])
    tmin_90d = np.copy(h5infile['tmin_90d'])
    tavg_90d = np.copy(h5infile['tavg_90d'])
    tmax_90d = np.copy(h5infile['tmax_90d'])
    prcp_30d = np.copy(h5infile['prcp_30d'])
    prcp_90d = np.copy(h5infile['prcp_90d'])
    prcp_90d_nd0 = np.copy(h5infile['prcp_90d_nd0'])
    prcp_90d_nd10 = np.copy(h5infile['prcp_90d_nd10'])
    prcp_90d_nd25 = np.copy(h5infile['prcp_90d_nd25'])
    prcp_180d = np.copy(h5infile['prcp_180d'])
    prcp_365d_p1 = np.copy(h5infile['prcp_365d_p1'])
    prcp_365d_p2 = np.copy(h5infile['prcp_365d_p2'])
    h5infile.close()
    message('extracting rolling accounting variable station lists')
    chill_d_stns = get_stn_lists(path, prev_year, 'chill_d_stns')
    chill_dd_stns = get_stn_lists(path, prev_year, 'chill_dd_stns')
    grow_dd_stns = get_stn_lists(path, prev_year, 'grow_dd_stns')
    grow_dd_base0_stns = get_stn_lists(path, prev_year, 'grow_dd_base0_stns')
    tmin_90d_stns = get_stn_lists(path, prev_year, 'tmin_90d_stns')
    tavg_90d_stns = get_stn_lists(path, prev_year, 'tavg_90d_stns')
    tmax_90d_stns = get_stn_lists(path, prev_year, 'tmax_90d_stns')
    prcp_30d_stns = get_stn_lists(path, prev_year, 'prcp_30d_stns')
    prcp_90d_stns = get_stn_lists(path, prev_year, 'prcp_90d_stns')
    prcp_180d_stns = get_stn_lists(path, prev_year, 'prcp_180d_stns')
    prcp_365d_stns = get_stn_lists(path, prev_year, 'prcp_365d_stns') 
else:  # otherwise, initialize those variable spaces
    h5infname = '%s/%d_NCEI_grids_1.h5' % (path, dates[0])
    message('extracting grid information from %s' % h5infname)
    h5infile = hdf.File(h5infname,'r')
    nrows = np.copy(h5infile['grid/nrows'])
    ncols = np.copy(h5infile['grid/ncols'])
    h5infile.close()
    message('establishing rolling accounting variable grids, datacubes, and station lists')
    tmin_frz_prev = np.zeros((nrows,ncols))
    tavg_frz_prev = np.zeros((nrows,ncols))
    tmax_frz_prev = np.zeros((nrows,ncols))
    chill_d_prev = np.zeros((nrows,ncols))
    chill_dd_prev = np.zeros((nrows,ncols))
    grow_dd_prev = np.zeros((nrows,ncols))
    grow_dd_base0_prev = np.zeros((nrows,ncols))
    tmin_90d = np.zeros((90,nrows,ncols))
    tavg_90d = np.zeros((90,nrows,ncols))
    tmax_90d = np.zeros((90,nrows,ncols))
    prcp_30d = np.zeros((30,nrows,ncols))
    prcp_90d = np.zeros((90,nrows,ncols))
    prcp_90d_nd0 = np.zeros((90,nrows,ncols))
    prcp_90d_nd10 = np.zeros((90,nrows,ncols))
    prcp_90d_nd25 = np.zeros((90,nrows,ncols))
    prcp_180d = np.zeros((180,nrows,ncols))
    prcp_365d_p1 = np.zeros((183,nrows,ncols))
    prcp_365d_p2 = np.zeros((182,nrows,ncols))
    chill_d_stns = []
    chill_dd_stns = []
    grow_dd_stns = []
    grow_dd_base0_stns = []
    tmin_90d_stns = []
    tavg_90d_stns = []
    tmax_90d_stns = [] 
    prcp_30d_stns = []
    prcp_90d_stns = []
    prcp_180d_stns = []
    prcp_365d_stns = [] 
message(' ')
#
for date in dates:
    h5infname = '%s/%d_NCEI_grids_1.h5' % (path, date)
    message('extracting PRCP, TMAX, TMIN, and TAVG grids from %s' % h5infname)
    h5infile = hdf.File(h5infname,'r')
    prcp_stns = np.copy(h5infile['stns/prcp_stns'])
    prcp = np.copy(h5infile['grid_prcp'])
    tmax_stns = np.copy(h5infile['stns/tmax_stns'])
    tmax = np.copy(h5infile['grid_tmax'])
    tmin_stns = np.copy(h5infile['stns/tmin_stns'])
    tmin = np.copy(h5infile['grid_tmin'])
    tavg_stns = np.copy(h5infile['stns/tavg_stns'])
    tavg = np.copy(h5infile['grid_tavg'])
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
    message('  mean CD %.1f    mean CDD %.1f    mean GDD %.1f    mean GDD_base0 %.1f' % (np.mean(chill_d),np.mean(chill_dd),np.mean(grow_dd),np.mean(grow_dd_base0)))
    #
    grid_tmin_frz = grid_threshold_count(month,day,cd_mm_start,cd_dd_start,tmin,tmin_frz_prev)
    tmin_frz_prev = grid_tmin_frz
    message('- calculated Tmin freezing days (accumulated from %s) mean %.1f' % (cd_start_str, np.mean(grid_tmin_frz)))
    #
    grid_tavg_frz = grid_threshold_count(month,day,cd_mm_start,cd_dd_start,tavg,tavg_frz_prev)
    tavg_frz_prev = grid_tavg_frz
    message('- calculated Tavg freezing days (accumulated from %s) mean %.1f' % (cd_start_str, np.mean(grid_tavg_frz)))
    #
    grid_tmax_frz = grid_threshold_count(month,day,cd_mm_start,cd_dd_start,tmax,tmax_frz_prev)
    tmax_frz_prev = grid_tmax_frz
    message('- calculated Tmax freezing days (accumulated from %s) mean %.1f' % (cd_start_str, np.mean(grid_tmax_frz)))
    #
    grid_chill_d, chill_d_stns = grid_threshold_count(month,day,cd_mm_start,cd_dd_start,chill_d,chill_d_prev,chill_d_stns,tavg_stns)
    chill_d_prev = grid_chill_d
    message('- calculated chilling days (based on Tavg, accumulated from %s) mean %.1f' % (cd_start_str, np.mean(grid_chill_d)))
    #
    grid_chill_dd, chill_dd_stns = grid_threshold_accumulate(month,day,cd_mm_start,cd_dd_start,chill_dd,chill_dd_prev,chill_dd_stns,tavg_stns)
    chill_dd_prev = grid_chill_dd
    message('- calculated chilling degree days (based on Tavg, accumulated from %s) mean %.1f' % (cd_start_str, np.mean(grid_chill_dd)))
    #
    grid_grow_dd, grow_dd_stns = grid_threshold_accumulate(month,day,gd_mm_start,gd_dd_start,grow_dd,grow_dd_prev,grow_dd_stns,tavg_stns)
    grow_dd_prev = grid_grow_dd
    message('- calculated growing degree days with base 5dC (based on Tavg, accumulated from %s) mean %.1f' % (gd_start_str, np.mean(grid_grow_dd)))
    #
    grid_grow_dd_base0, grow_dd_base0_stns = grid_threshold_accumulate(month,day,gd_mm_start,gd_dd_start,grow_dd_base0,grow_dd_base0_prev,grow_dd_base0_stns,tavg_stns)
    grow_dd_base0_prev = grid_grow_dd_base0
    message('- calculated growing degree days with base 0dC (based on Tavg, accumulated from %s) mean %.1f' % (gd_start_str, np.mean(grid_grow_dd_base0)))
    #
    grid_tmin_90d_mean, grid_tmin_90d_var, tmin_90d_stns_all, tmin_90d, tmin_90d_stns = cube_mean_var(90,tmin_90d,tmin,tmin_90d_stns,tmin_stns)
    message('- calculated updated 90-day running Tmin mean and variance, mean %.1f' % np.mean(grid_tmin_90d_mean))
    #
    grid_tavg_90d_mean, grid_tavg_90d_var, tavg_90d_stns_all, tavg_90d, tavg_90d_stns = cube_mean_var(90,tavg_90d,tavg,tavg_90d_stns,tavg_stns)
    message('- calculated updated 90-day running Tavg mean and variance, mean %.1f' % np.mean(grid_tavg_90d_mean))
    #
    grid_tmax_90d_mean, grid_tmax_90d_var, tmax_90d_stns_all, tmax_90d, tmax_90d_stns = cube_mean_var(90,tmax_90d,tmax,tmax_90d_stns,tmax_stns)
    message('- calculated updated 90-day running Tmax mean and variance, mean %.1f' % np.mean(grid_tmax_90d_mean))
    #
    grid_prcp_30d, prcp_30d_stns_all, prcp_30d, prcp_30d_stns = cube_sum(30,prcp_30d,prcp,prcp_30d_stns,prcp_stns)
    message('- calculated updated 30-day running precipitation total, mean %.1f' % np.mean(grid_prcp_30d))
    #
    grid_prcp_90d, prcp_90d_stns_all, prcp_90d, prcp_90d_stns = cube_sum(90,prcp_90d,prcp,prcp_90d_stns,prcp_stns)
    message('- calculated updated 90-day running precipitation total, mean %.1f' % np.mean(grid_prcp_90d))
    #
    grid_prcp_180d, prcp_180d_stns_all, prcp_180d, prcp_180d_stns = cube_sum(180,prcp_180d,prcp,prcp_180d_stns,prcp_stns)
    message('- calculated updated 180-day running precipitation total, mean %.1f' % np.mean(grid_prcp_180d))
    #
    grid_prcp_365d, prcp_365d_stns_all, prcp_365d_p1, prcp_365d_p2, prcp_365d_stns = cube_sum_parts(365,prcp_365d_p1,prcp_365d_p2,prcp,prcp_365d_stns,prcp_stns)
    message('- calculated updated 365-day running precipitation total, mean %.1f' % np.mean(grid_prcp_365d))
    #
    grid_prcp_90d_nd0, prcp_90d_nd0 = cube_threshold_count(90,prcp_90d_nd0,prcp,0.0)
    message('- calculated updated 90-day running precipitation ndays total for P > 0 mm, mean %.1f' % np.mean(grid_prcp_90d_nd0))
    #
    grid_prcp_90d_nd10, prcp_90d_nd10 = cube_threshold_count(90,prcp_90d_nd10,prcp,1.0)
    message('- calculated updated 90-day running precipitation ndays total for P > 10 mm, mean %.1f' % np.mean(grid_prcp_90d_nd10))
    #
    grid_prcp_90d_nd25, prcp_90d_nd25 = cube_threshold_count(90,prcp_90d_nd25,prcp,2.5)
    message('- calculated updated 90-day running precipitation ndays total for P > 25 mm, mean %.1f' % np.mean(grid_prcp_90d_nd25))
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
# save rolling accounting variables for next year's run
varfname = '%s/%d_year_end_vars.h5' % (path, this_year)
message('saving rolling accounting variable grids and datacubes to %s' % varfname)
h5outfile = hdf.File(varfname,'w')
h5outfile.create_dataset('nrows', data=nrows)
h5outfile.create_dataset('ncols', data=ncols)
h5outfile.create_dataset('tmin_frz_prev', data=tmin_frz_prev, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('tavg_frz_prev', data=tavg_frz_prev, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('tmax_frz_prev', data=tmax_frz_prev, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('chill_d_prev', data=chill_d_prev, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('chill_dd_prev', data=chill_dd_prev, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('grow_dd_prev', data=grow_dd_prev, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('grow_dd_base0_prev', data=grow_dd_base0_prev, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('tmin_90d', data=tmin_90d, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('tavg_90d', data=tavg_90d, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('tmax_90d', data=tmax_90d, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('prcp_30d', data=prcp_30d, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('prcp_90d', data=prcp_90d, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('prcp_90d_nd0', data=prcp_90d_nd0, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('prcp_90d_nd10', data=prcp_90d_nd10, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('prcp_90d_nd25', data=prcp_90d_nd25, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('prcp_180d', data=prcp_180d, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('prcp_365d_p1', data=prcp_365d_p1, dtype=np.float32, compression='gzip')
h5outfile.create_dataset('prcp_365d_p2', data=prcp_365d_p2, dtype=np.float32, compression='gzip')
h5outfile.close()
message('saving rolling accounting variable station lists')
write_stn_lists(path, this_year, 'chill_d_stns', chill_d_stns)
write_stn_lists(path, this_year, 'chill_dd_stns', chill_dd_stns)
write_stn_lists(path, this_year, 'grow_dd_stns', grow_dd_stns)
write_stn_lists(path, this_year, 'grow_dd_base0_stns', grow_dd_base0_stns)
write_stn_lists(path, this_year, 'tmin_90d_stns', tmin_90d_stns)
write_stn_lists(path, this_year, 'tavg_90d_stns', tavg_90d_stns)
write_stn_lists(path, this_year, 'tmax_90d_stns', tmax_90d_stns)
write_stn_lists(path, this_year, 'prcp_30d_stns', prcp_30d_stns)
write_stn_lists(path, this_year, 'prcp_90d_stns', prcp_90d_stns)
write_stn_lists(path, this_year, 'prcp_180d_stns', prcp_180d_stns)
write_stn_lists(path, this_year, 'prcp_365d_stns', prcp_365d_stns) 
#
message('process_NCEI_03_serial_checkpointed.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_03_serial_checkpointed.py
