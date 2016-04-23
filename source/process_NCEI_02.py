"""
Python script 'process_NCEI_02.py'
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

USAGE: 'python process_NCEI_02.py NLCD_2011_WLS_UTM15N NCEI_WLS_19840101-20131231 ./grids 480 RBF 1'

NOTES: [NCEI_WLS_19840101-20131231] is the '_processed.h5' file prefix in your 'data/' directory
       [480] is the (default) output grid resolution in meters
       [RBF] is the (default) interpolation method (see below for choices)
       [1] is the (default) flag to plot the resulting daily grids 

PURPOSE: Gridded interpolation of daily PRCP/TMAX/TMIN station data from cleaned NOAA/NCEI 
         dataset via user's method of choice. 
         RBF: radial basis functions, a SciPy built-in method
         CSP: cubic splines via griddata, a SciPy built-in method
         BSP: bivariate cubic B-splines, a SciPy built-in method
         IDW: inverse-distance-squared for temperature and -cubed for precipitation

NOTES: RBF is likely the most economical and useful, as it interpolates over the whole domain 
       at once using only 3 processors (one for each input variable) with an excellent daily 
       processing time. CSP is just about as fast but may not cover the whole desired domain,
       forming a convex hull around the available stations for each grid, can get weird in
       areas with low station density, and produces some weird patterns between stations. 
       BSP does cover the whole desired domain but also gets weird in areas with low station 
       density, especially near the edges of the domain. IDW is likely the least economical, 
       as it tiles the desired domain and uses 9 processors in order to achieve the same 
       daily processing time as the other methods, and maps often show isolated bulls-eye 
       patterns that are difficult to smooth out.

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for I/O
			  The pp module is ParallelPython, employed here for multi-processor application
              The 'UTM_Geo_Convert', 'Interpolation', 'Plots' modules have their own requirements

NOTE: This script uses 1 ParallelPython server with a number of CPUs that depends on the 
      selected interpolation method (9 for IDW, 3 for all others). Make sure the total falls 
      within your number of available processors 

INPUT: Output file from process_NCEI_01.py script in '.h5' format (one file)
       (with the naming convention 'data/NCEI_WLS_[YYYYMMDD]-[YYYYMMDD]_processed.h5')
       A header file corresponding to your study region's clipped NLCD 30m map boundaries

NOTE: If you don't have an NLCD map of your study region, you can spoof the '.hdr' file
      by modifying the provided example. Note the variables that are necessary for this
      to work: UTMzone, nrows, ncols, min_y, max_y, min_x, max_x, dy, dx. You only need to 
      change NROWS, NCOLS, ULXMAP, ULYMAP, XDIM, and YDIM for this (UTMzone is in the file 
      name, min_y and max_x are calculated). Other information like BYTEORDER, NBANDS, etc. 
      are not important, but still need to be present with values that will be read by the 
      external routine. 

OUTPUT: Daily '.h5' file with original meteorological data and gridded fields (1 / day)
        (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_1.h5')
        Corresponding plotted fields in 'images/*.png' (4 / day) 

RUN TIME: ~55 mins per data year using RBF interpolation method
"""

import os, sys, datetime
import h5py as hdf
import numpy as np
import pp
import Interpolation
from UTM_Geo_Convert import *
from Read_Header_Files import *
from Plots import p_map_plot, t_map_plot

def message(str):
    print str
    sys.stdout.flush()
    return    

def get_stn_coords(lats,lons,UTMz):
    easts = []
    norths = []
    for i in range(0,len(lats)):
        UTM_coords = geographic_to_utm(lons[i].astype(np.float64), lats[i].astype(np.float64), UTMz)
        easts.append(int(round(UTM_coords[1],0)))
        norths.append(int(round(UTM_coords[2],0)))
    return easts, norths

def check_duplicates(stns,east,north,vals):
    nstns = len(stns)
    # first remove known bad stations
    bad_stns = ['GHCND:USC00214884','GHCND:CA006041109']
    bad_idxs = []
    for i in range(0,nstns):
        if stns[i] in bad_stns:
            bad_idxs.append(i)
    nbad = len(bad_idxs)
    if nbad > 0:
        bad_idxs = sorted(bad_idxs, reverse=True)
        for idx in bad_idxs:
            stns = np.delete(stns,idx)
            east = np.delete(east,idx)
            north = np.delete(north,idx)
            vals = np.delete(vals,idx)
    # now check for duplicate locations
    nstns = len(stns)
    ndups = 0
    e_n_pairs = []
    for i in range(0,nstns):
        e_n_pairs.append((east[i],north[i]))
    unique_pairs = list(set(e_n_pairs))
    if nstns != len(unique_pairs):
        dup_loc_idxs = []
        for i in range(0,len(unique_pairs)):
            n_entries = 0
            for j in range(0,len(e_n_pairs)):
                if e_n_pairs[j] == unique_pairs[i]:
                    n_entries += 1
                    if n_entries > 1:
                        dup_loc_idxs.append(j)
        ndups = len(dup_loc_idxs)
        if ndups > 0:
            dup_loc_idxs = sorted(dup_loc_idxs, reverse=True)
            for idx in dup_loc_idxs:
                stns = np.delete(stns,idx)
                east = np.delete(east,idx)
                north = np.delete(north,idx)
                vals = np.delete(vals,idx)
    return stns, east, north, vals, ndups, nbad

def process_tile_idw(tid,t_grid_x,t_grid_y,t_shape,stns_e,stns_n,stns_vals,exp):
    t_var_flat = Interpolation.simple_idw(stns_e, stns_n, stns_vals, t_grid_x, t_grid_y, exp)
    t_var = t_var_flat.reshape(t_shape)
    return tid, t_var

def process_grid_csp(grid_x,grid_y,grid_shape,stns_e,stns_n,stns_vals):
    grid_var = Interpolation.scipy_griddata(stns_e, stns_n, stns_vals, grid_x, grid_y)
    grid_var = grid_var.T
    return grid_var

def process_grid_bsp(grid_x,grid_y,grid_shape,stns_e,stns_n,stns_vals):
    grid_var = Interpolation.scipy_bspline(stns_e, stns_n, stns_vals, grid_x, grid_y)
    return grid_var

def process_grid_rbf(grid_x,grid_y,grid_shape,stns_e,stns_n,stns_vals):
    grid_var = Interpolation.scipy_rbf(stns_e, stns_n, stns_vals, grid_x, grid_y)
    return grid_var

def process_tiles(jobs,ts,tshape,grid_var):
    for job in jobs:
        returns = job()
        tnum = returns[0]     # tile number
        t_var = returns[1]    # tile grid
        t_minc = ts[tnum][4]  # tile min full-grid column number
        t_maxc = ts[tnum][5]  # tile max full-grid column number
        t_minr = ts[tnum][6]  # tile min full-grid row number
        t_maxr = ts[tnum][7]  # tile max full-grid row number
        t_nr = tshape[0]      # tile rows
        t_nc = tshape[1]      # tile columns        
        grid_var[t_minr:t_maxr, t_minc:t_maxc+1] = t_var[0:t_nr, 0:t_nc] 
    return grid_var

def process_grid(jobs):
    for job in jobs:
        grid_var = job()
    return grid_var

message(' ')
message('process_NCEI_02.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 7:
    message('no plot flag indicated, setting plots = True')
    plots = 1
else:
    plots = int(sys.argv[6])
#
if len(sys.argv) < 6:
    message('no interpolation method indicated, using RBF (most economical choice)')
    interp_method = 'RBF'
elif sys.argv[4] in ['RBF','CSP','BSP','IDW']:
    interp_method = sys.argv[5]
else:
    message('interpolation method %s is not available; use RBF/CSP/BSP/IDW' % sys.argv[4])
    sys.exit(1)
#
if len(sys.argv) < 5:
    message('no output grid resolution indicated, setting output at 500m grid spacing')
    dx_out = dy_out = 500
else:
    dx_out = dy_out = int(sys.argv[4])
#
if len(sys.argv) < 4:
    message('input warning: no output grids directory path indicated, using ./grids')
    path = './grids'
else:
    path = sys.argv[3]
#
if len(sys.argv) < 3:
    message('input error: need prefix for h5 file containing NCEI weather data')
    sys.exit(1)
else:
    NCEIfname = sys.argv[2]
h5infname = '%s/../data/%s_processed.h5' % (path, NCEIfname)
#
if len(sys.argv) < 2:
    message('input error: need prefix for hdr file containing NLCD map boundaries')
    sys.exit(1)
else:
    NLCDhname = sys.argv[1]
NLCDfname = '%s/../data/%s.hdr' % (path, NLCDhname)
#
message('extracting header information from %s' % NLCDfname)
UTMzone, nrows, ncols, min_y, max_y, min_x, max_x, dy, dx = get_bil_hdr_info(NLCDfname)
if UTMzone == 0:
    message('error: UTM zone not found or improperly specified in file name')
    sys.exit(1)
message('- landcover grid information')
message('-- UTM zone %d' % UTMzone)
message('-- dx = dy = %d' % dx)
message('-- Y from %d to %d (%d rows)' % (min_y, max_y, nrows))
message('-- X from %d to %d (%d cols)' % (min_x, max_x, ncols))
dy = dy_out
min_y = min_y - dy
y_dist = max_y - min_y
nrows = (y_dist // dy) + 1
max_y = min_y + (nrows * dy)
dx = dx_out
min_x = min_x - dx
x_dist = max_x - min_x
ncols = (x_dist // dx) + 1
max_x = min_x + (ncols * dx)
message('- adjusted grid extent')
message('-- dx = dy = %d' % dx)
message('-- Y from %d to %d (%d rows)' % (min_y, max_y, nrows))
message('-- X from %d to %d (%d cols)' % (min_x, max_x, ncols))
message(' ')
#
if interp_method == 'RBF':
    # interpolate using the whole target grid
    message('generating target grid for radial basis function (RBF) interpolation')
    grid_x = np.arange(min_x, max_x, dx)
    grid_y = np.arange(min_y, max_y, dy)
    grid_x, grid_y = np.meshgrid(grid_x, grid_y)
    grid_shape = grid_x.shape
elif interp_method == 'CSP':
    # interpolate using the whole target grid
    message('generating target grid for 2D cubic spline (CSP) interpolation')
    grid_x, grid_y = np.mgrid[min_x:max_x:dx, min_y:max_y:dy]
    grid_shape = grid_x.shape
elif interp_method == 'BSP':
    # interpolate using the whole target grid
    message('generating target grid for bivariate cubic spline (BSP) interpolation')
    grid_x, grid_y = np.mgrid[min_x:max_x:dx, min_y:max_y:dy]
    grid_shape = grid_x.shape
    grid_x, grid_y = grid_x[:,0], grid_y[0,:]
elif interp_method == 'IDW':
    # interpolate on individual tiles
    ntiles = 9  # must be a perfect square (for now, at least)
    message('generating target grid for IDW interpolation in %d tile(s)' % ntiles)
    ntiledivs = int(np.sqrt(ntiles))
    x_tile_npts = x_dist // dx // ntiledivs
    y_tile_npts = y_dist // dy // ntiledivs
    tiles = []
    for i in range(0,ntiledivs):
        tile_min_x = min_x + (i * x_tile_npts * dx)
        tile_min_col = (tile_min_x - min_x) // dx
        tile_max_x = tile_min_x + (x_tile_npts * dx)
        tile_max_col = tile_min_col + x_tile_npts
        for j in range(0,ntiledivs):
            tile_min_y = min_y + (j * y_tile_npts * dy)
            tile_min_row = (tile_min_y - min_y) // dy
            tile_max_y = tile_min_y + (y_tile_npts * dy)
            tile_max_row = tile_min_row + y_tile_npts + 1
            tile_grid_x = np.arange(tile_min_x, tile_max_x+dx, dx)
            tile_grid_y = np.arange(tile_min_y, tile_max_y+dy, dy)
            tile_grid_x, tile_grid_y = np.meshgrid(tile_grid_x, tile_grid_y)
            tile_shape = tile_grid_x.shape
            tile_grid_x, tile_grid_y = tile_grid_x.flatten(), tile_grid_y.flatten()
            tiles.append([tile_min_x, tile_max_x, tile_min_y, tile_max_y, 
                          tile_min_col, tile_max_col, tile_min_row, tile_max_row,
                          tile_grid_x, tile_grid_y, tile_shape])
            message('- tile %d' % len(tiles))
            message('-- Y from %d to %d (rows %d to %d)' % (tile_min_y, tile_max_y, tile_min_row, tile_max_row))
            message('-- X from %d to %d (cols %d to %d)' % (tile_min_x, tile_max_x, tile_min_col, tile_max_col))
            message('-- tile grid dimensions %s' % str(tile_shape))
message(' ')
#
if interp_method in ['RBF','CSP','BSP']:
    ncpus = 3
elif interp_method == 'IDW':
    ncpus = ntiles
ppservers0 = ()
job_server0 = pp.Server(ncpus, ppservers=ppservers0)
message('starting ParallelPython (pp) server #0 at %s' % datetime.datetime.now().isoformat())
message('- using %d processors' % job_server0.get_ncpus())
if interp_method == 'RBF':
    prcp_interp = tmax_interp = tmin_interp = pp.Template(job_server0, process_grid_rbf, depfuncs=(), modules=("Interpolation",))
elif interp_method == 'CSP':
    prcp_interp = tmax_interp = tmin_interp = pp.Template(job_server0, process_grid_csp, depfuncs=(), modules=("Interpolation",))
elif interp_method == 'BSP':
    prcp_interp = tmax_interp = tmin_interp = pp.Template(job_server0, process_grid_bsp, depfuncs=(), modules=("Interpolation",))
elif interp_method == 'IDW':
    prcp_interp = tmax_interp = tmin_interp = pp.Template(job_server0, process_tile_idw, depfuncs=(), modules=("Interpolation",))
message(' ')
#
message('reading station and date information from %s' % h5infname)
h5infile = hdf.File(h5infname,'r')
stn_id = np.copy(h5infile['stn_id'])
dates = np.copy(h5infile['dates'])
h5infile.close()
message('- identifiers for %d stations found' % len(stn_id))
message('- meteorological data for %d dates found' % len(dates))
message(' ')
#
for date in dates:
    message('reading station information and met data for %d' % date)
    h5infile = hdf.File(h5infname,'r')
    datepath = 'metdata/%d' % date
    #
    prcp_stns = np.copy(h5infile[datepath + '/prcp_stns'])
    prcp_lat = np.copy(h5infile[datepath + '/prcp_lat'])
    prcp_lon = np.copy(h5infile[datepath + '/prcp_lon'])
    prcp_east, prcp_north = get_stn_coords(prcp_lat, prcp_lon, UTMzone)
    prcp_vals = np.copy(h5infile[datepath + '/prcp_vals'])
    prcp_stns,prcp_east,prcp_north,prcp_vals,nduplicates,nbadstns = check_duplicates(prcp_stns,prcp_east,prcp_north,prcp_vals)
    message('- PRCP data for %d unique stations found' % len(prcp_stns))
    if nbadstns > 0:
        message('-- %d known bad stations were removed' % nbadstns)
    if nduplicates > 0:
        message('-- %d duplicate locations were removed' % nduplicates)
    #
    tmax_stns = np.copy(h5infile[datepath + '/tmax_stns'])
    tmax_lat = np.copy(h5infile[datepath + '/tmax_lat'])
    tmax_lon = np.copy(h5infile[datepath + '/tmax_lon'])
    tmax_east, tmax_north = get_stn_coords(tmax_lat, tmax_lon, UTMzone)
    tmax_vals = np.copy(h5infile[datepath + '/tmax_vals'])
    tmax_stns,tmax_east,tmax_north,tmax_vals,nduplicates,nbadstns = check_duplicates(tmax_stns,tmax_east,tmax_north,tmax_vals)
    message('- TMAX data for %d unique stations found' % len(tmax_stns))
    if nbadstns > 0:
        message('-- %d known bad stations were removed' % nbadstns)
    if nduplicates > 0:
        message('-- %d duplicate locations were removed' % nduplicates)
    #
    tmin_stns = np.copy(h5infile[datepath + '/tmin_stns'])
    tmin_lat = np.copy(h5infile[datepath + '/tmin_lat'])
    tmin_lon = np.copy(h5infile[datepath + '/tmin_lon'])
    tmin_east, tmin_north = get_stn_coords(tmin_lat, tmin_lon, UTMzone)
    tmin_vals = np.copy(h5infile[datepath + '/tmin_vals'])
    tmin_stns,tmin_east,tmin_north,tmin_vals,nduplicates,nbadstns = check_duplicates(tmin_stns,tmin_east,tmin_north,tmin_vals)
    message('- TMIN data for %d unique stations found' % len(tmin_stns))
    if nbadstns > 0:
        message('-- %d known bad stations were removed' % nbadstns)
    if nduplicates > 0:
        message('-- %d duplicate locations were removed' % nduplicates)
    h5infile.close()
    #
    message('interpolating meteorological values for %d' % date)
    jobs_prcp_interp = []
    jobs_tmax_interp = []
    jobs_tmin_interp = []
    if interp_method in ['RBF','CSP','BSP']:
        jobs_prcp_interp.append(prcp_interp.submit(grid_x,grid_y,grid_shape,prcp_east,prcp_north,prcp_vals))
        jobs_tmax_interp.append(tmax_interp.submit(grid_x,grid_y,grid_shape,tmax_east,tmax_north,tmax_vals))
        jobs_tmin_interp.append(tmin_interp.submit(grid_x,grid_y,grid_shape,tmin_east,tmin_north,tmin_vals))
        #
        grid_prcp = process_grid(jobs_prcp_interp)
        grid_prcp = np.where(grid_prcp < 0.0, 0.0, grid_prcp) # ensure prcp is a positive definite field
        message('- retrieved interpolated PRCP grid')
        #
        grid_tmax = process_grid(jobs_tmax_interp)
        message('- retrieved interpolated TMAX grid')
        #
        grid_tmin = process_grid(jobs_tmin_interp)
        message('- retrieved interpolated TMIN grid')
    elif interp_method == 'IDW':
        for t in range(0,ntiles):
            tile_grid_x = tiles[t][8]
            tile_grid_y = tiles[t][9]
            tile_shape = tiles[t][10]
            jobs_prcp_interp.append(prcp_interp.submit(t,tile_grid_x,tile_grid_y,tile_shape,prcp_east,prcp_north,prcp_vals,3))
            jobs_tmax_interp.append(tmax_interp.submit(t,tile_grid_x,tile_grid_y,tile_shape,tmax_east,tmax_north,tmax_vals,2))
            jobs_tmin_interp.append(tmin_interp.submit(t,tile_grid_x,tile_grid_y,tile_shape,tmin_east,tmin_north,tmin_vals,2))
        message('- %d tiles queued for interpolation' % ntiles)
        #
        grid_prcp = np.zeros((nrows,ncols))
        grid_prcp = process_tiles(jobs_prcp_interp, tiles, tile_shape, grid_prcp)
        grid_prcp = np.where(grid_prcp < 0.0, 0.0, grid_prcp) # ensure prcp is a positive definite field
        message('- retrieved interpolated PRCP grid tiles')
        #
        grid_tmax = np.zeros((nrows,ncols))
        grid_tmax = process_tiles(jobs_tmax_interp, tiles, tile_shape, grid_tmax)
        message('- retrieved interpolated TMAX grid tiles')
        #
        grid_tmin = np.zeros((nrows,ncols))
        grid_tmin = process_tiles(jobs_tmin_interp, tiles, tile_shape, grid_tmin)
        message('- retrieved interpolated TMIN grid tiles')
    #
    grid_tavg = (grid_tmin + grid_tmax) / 2.0
    tavg_stns = list(set(tmax_stns) | set(tmin_stns))
    message('- calculated TAVG grid and station set')
    #
    h5outfname = '%s/%d_NCEI_grids_1.h5' % (path, date)
    message('writing grids to %s' % h5outfname)
    sys.stdout.flush()
    h5outfile = hdf.File(h5outfname,'w')
    h5outfile.create_dataset('meta/filename', data=h5outfname)
    h5outfile.create_dataset('meta/created', data=datetime.datetime.now().isoformat())
    h5outfile.create_dataset('meta/by', data='M. Garcia, UWisconsin-Madison FWE')
    h5outfile.create_dataset('meta/last_updated', data=datetime.datetime.now().isoformat())
    h5outfile.create_dataset('meta/at', data='prcp/tmax/tmin/tavg grids')
    message('- saved processing metadata items')
    #
    h5outfile.create_dataset('grid/UTMzone', data=UTMzone)
    h5outfile.create_dataset('grid/min_x', data=min_x)
    h5outfile.create_dataset('grid/max_x', data=max_x)
    h5outfile.create_dataset('grid/x_dist', data=x_dist)
    h5outfile.create_dataset('grid/dx', data=dx)
    h5outfile.create_dataset('grid/ncols', data=ncols)
    h5outfile.create_dataset('grid/min_y', data=min_y)
    h5outfile.create_dataset('grid/max_y', data=max_y)
    h5outfile.create_dataset('grid/y_dist', data=y_dist)
    h5outfile.create_dataset('grid/dy', data=dy)
    h5outfile.create_dataset('grid/nrows', data=nrows)
    h5outfile.create_dataset('grid/interp_method', data=interp_method)
    message('- saved grid definition metadata items')
    #
    h5outfile.create_dataset('stns/stn_id', data=stn_id)
    h5outfile.create_dataset('stns/prcp_stns', data=prcp_stns)
    h5outfile.create_dataset('stns/prcp_lat', data=prcp_lat)
    h5outfile.create_dataset('stns/prcp_lon', data=prcp_lon)
    h5outfile.create_dataset('stns/prcp_easting', data=prcp_east)
    h5outfile.create_dataset('stns/prcp_northing', data=prcp_north)
    h5outfile.create_dataset('stns/prcp_value', data=prcp_vals)
    h5outfile.create_dataset('stns/tmax_stns', data=tmax_stns)
    h5outfile.create_dataset('stns/tmax_lat', data=tmax_lat)
    h5outfile.create_dataset('stns/tmax_lon', data=tmax_lon)
    h5outfile.create_dataset('stns/tmax_easting', data=tmax_east)
    h5outfile.create_dataset('stns/tmax_northing', data=tmax_north)
    h5outfile.create_dataset('stns/tmax_value', data=tmax_vals)
    h5outfile.create_dataset('stns/tmin_stns', data=tmin_stns)
    h5outfile.create_dataset('stns/tmin_lat', data=tmin_lat)
    h5outfile.create_dataset('stns/tmin_lon', data=tmin_lon)
    h5outfile.create_dataset('stns/tmin_easting', data=tmin_east)
    h5outfile.create_dataset('stns/tmin_northing', data=tmin_north)
    h5outfile.create_dataset('stns/tmin_value', data=tmin_vals)
    h5outfile.create_dataset('stns/tavg_stns', data=tavg_stns)
    message('- saved input station data items')
    #
    h5outfile.create_dataset('grid_prcp', data=grid_prcp, dtype=np.float32, compression='gzip')
    message('- saved PRCP grid %s with %d stations' % (str(grid_prcp.shape), len(prcp_stns)))
    h5outfile.create_dataset('grid_tmax', data=grid_tmax, dtype=np.float32, compression='gzip')
    message('- saved TMAX grid %s with %d stations' % (str(grid_tmax.shape), len(tmax_stns)))
    h5outfile.create_dataset('grid_tmin', data=grid_tmin, dtype=np.float32, compression='gzip')
    message('- saved TMIN grid %s with %d stations' % (str(grid_tmin.shape), len(tmin_stns)))
    h5outfile.create_dataset('grid_tavg', data=grid_tavg, dtype=np.float32, compression='gzip')
    message('- saved TAVG grid %s with %d stations' % (str(grid_tavg.shape), len(tavg_stns)))
    h5outfile.close()
    #
    if plots:
        message('plotting grids')
        titlestr = '%d Precip (cm) via %s' % (date, interp_method)
        filename = '%s/../images/%d_%s_prcp.png' % (path, date, interp_method.lower())
        p_map_plot(prcp_east, prcp_north, prcp_vals, min_x, max_x, min_y, max_y, grid_prcp, UTMzone, titlestr, filename)
        titlestr = '%d Tmax (%sC) via %s' % (date, r'$^\circ$', interp_method)
        filename = '%s/../images/%d_%s_tmax.png' % (path, date, interp_method.lower())
        t_map_plot(tmax_east, tmax_north, tmax_vals, min_x, max_x, min_y, max_y, grid_tmax, UTMzone, titlestr, filename)
        titlestr = '%d Tmin (%sC) via %s' % (date, r'$^\circ$', interp_method)
        filename = '%s/../images/%d_%s_tmin.png' % (path, date, interp_method.lower())
        t_map_plot(tmin_east, tmin_north, tmin_vals, min_x, max_x, min_y, max_y, grid_tmin, UTMzone, titlestr, filename)
        titlestr = '%d Tavg (%sC) via %s' % (date, r'$^\circ$', interp_method)
        filename = '%s/../images/%d_%s_tavg.png' % (path, date, interp_method.lower())
        t_map_plot(tmin_east, tmin_north, tmin_vals, min_x, max_x, min_y, max_y, grid_tavg, UTMzone, titlestr, filename, stations=0)
    message(' ')
#
job_server0.destroy()
message('server #0 parallel processing completed at %s' % datetime.datetime.now().isoformat())
message(' ')
#
message('process_NCEI_02.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_02.py
