"""
Python script 'process_NCEI_04a.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Aggregating grids for specific dates and time periods over the
         analysis period, e.g. CD values at VEQ (vernal equinox) and other
         seasonal boundaries, finding CD between VEQ and SSOL, finding
         beginning and end of CD plateau, calculating length (both days and
         GDD) of CD plateau. Numerous variables are addressed.

DEPENDENCIES: h5py, numpy
              The 'Date_Convert' module has no external requirements

USAGE: '$ python process_NCEI_04a.py 1984 2013 ./grids'

INPUT: '.h5' output files from process_NCEI_03.py
       (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_2.h5')

OUTPUT: '.h5' file with aggregated grid datacubes
        (with the naming convention
         'analyses/[YYYY]-[YYYY]_derived_clim_grids.h5')
"""


import sys
import datetime
import glob
import h5py as hdf
import numpy as np
from Date_Convert import date_to_doy


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


def getgrids(fname, nr, nc, invars, n=0):
    if n == 0:
        nvars = len(invars)
    else:
        nvars = n
    grids = np.zeros((nvars, nr, nc))
    with hdf.File(fname, 'r') as h5file:
        for k, var in enumerate(invars):
            grids[k, :, :] = np.copy(h5file[var])
    return grids


def getpgrid(fname, nr, nc):
    grid = np.zeros((nr, nc))
    with hdf.File(fname, 'r') as h5file:
        grid[:, :] = np.copy(h5file['prcp_365d_sum'])
    return grid


def write_to_file(h5file, gvar, gdata):
    if gvar in h5file.keys():
        del h5file[gvar]
    h5file.create_dataset(gvar, data=gdata, dtype=np.float32,
                          compression='gzip')
    message('- %s %s' % (gvar, str(gdata.shape)))
    return


message(' ')
message('process_NCEI_04a.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no input directory path indicated, using ./grids')
    path = './grids'
else:
    path = sys.argv[3]
#
if len(sys.argv) < 3:
    message('no dates specified, analyzing 1984-2013 period')
    year_begin = 1984
    year_end = 2013
else:
    year_begin = int(sys.argv[1])
    year_end = int(sys.argv[2])
years = np.arange(year_begin, year_end + 1).astype(int)
nyears = len(years)
date_vequinox = 320
date_ssolstice = 621
date_aequinox = 922
date_wsolstice = 1221
date_eoy = 1231
#
wxlist = glob.glob('%s/*_NCEI_grids_2.h5' % path)
message('found %d weather derivative grid files' % len(wxlist))
message(' ')
#
infile = wxlist[0]
message('extracting grid information from weather derivatives file %s' %
        infile)
with hdf.File(infile, 'r') as h5infile:
    UTM_zone = np.copy(h5infile['grid/UTMzone'])
    NWeasting = np.copy(h5infile['grid/min_x'])
    NWnorthing = np.copy(h5infile['grid/max_y'])
    SEeasting = np.copy(h5infile['grid/max_x'])
    SEnorthing = np.copy(h5infile['grid/min_y'])
    dx = np.copy(h5infile['grid/dx'])
    dy = np.copy(h5infile['grid/dy'])
    ncols = np.copy(h5infile['grid/ncols'])
    nrows = np.copy(h5infile['grid/nrows'])
message('- UTM zone %d' % UTM_zone)
message('- Y from %d to %d (%d rows)' % (SEnorthing, NWnorthing, nrows))
message('- X from %d to %d (%d cols)' % (NWeasting, SEeasting, ncols))
message('- dx = dy = %d' % dx)
message(' ')
UTM_bounds = [NWeasting, NWnorthing, SEeasting, SEnorthing]
#
outfile = '%s/../analyses/%d-%d_derived_clim_grids.h5' % \
    (path, year_begin, year_end)
message('writing metadata and grid information to %s' % outfile)
with hdf.File(outfile, 'w') as h5outfile:
    h5outfile.create_dataset('meta/filename', data=outfile)
    h5outfile.create_dataset('meta/created',
                             data=datetime.datetime.now().isoformat())
    h5outfile.create_dataset('meta/by',
                             data='M. Garcia, UWisconsin-Madison FWE')
    h5outfile.create_dataset('meta/last_updated',
                             data=datetime.datetime.now().isoformat())
    h5outfile.create_dataset('meta/at',
                             data='derived climatological annual and \
                                  statistics grids')
    message('- 5 metadata items saved')
    h5outfile.create_dataset('grid/UTM_zone', data=UTM_zone)
    h5outfile.create_dataset('grid/min_x', data=NWeasting)
    h5outfile.create_dataset('grid/max_x', data=SEeasting)
    h5outfile.create_dataset('grid/dx', data=dx)
    h5outfile.create_dataset('grid/ncols', data=ncols)
    h5outfile.create_dataset('grid/min_y', data=SEnorthing)
    h5outfile.create_dataset('grid/max_y', data=NWnorthing)
    h5outfile.create_dataset('grid/dy', data=dy)
    h5outfile.create_dataset('grid/nrows', data=nrows)
    message('- 9 grid parameters saved')
    h5outfile.create_dataset('year_begin', data=year_begin)
    h5outfile.create_dataset('year_end', data=year_end)
    h5outfile.create_dataset('years', data=years)
    message('- 3 series parameters saved')
h5outfile.close()
message(' ')
#
inputgridvars = ['chill_d', 'chill_dd', 'grow_dd', 'grow_dd_base0',
                 'tmin_frz_days', 'tavg_frz_days', 'tmax_frz_days',
                 'tmin_90d_avg', 'tmin_90d_var', 'tavg_90d_avg',
                 'tavg_90d_var', 'tmax_90d_avg', 'tmax_90d_var',
                 'prcp_90d_sum', 'prcp_90d_nd0_sum', 'prcp_90d_nd10_sum',
                 'prcp_90d_nd25_sum']
ninputgridvars = len(inputgridvars)
#
# establish variable grids at EQ and SOL dates
inputgrids_at_veq = np.zeros((ninputgridvars, nyears, nrows, ncols))
inputgrids_at_ssol = np.zeros((ninputgridvars, nyears, nrows, ncols))
inputgrids_at_aeq = np.zeros((ninputgridvars, nyears, nrows, ncols))
inputgrids_at_wsol = np.zeros((ninputgridvars, nyears, nrows, ncols))
#
grids_intensity_winter = np.zeros((nyears, nrows, ncols))
grids_cd_veq_to_ssol = np.zeros((nyears, nrows, ncols))
grids_cdd_veq_to_ssol = np.zeros((nyears, nrows, ncols))
grids_gdd_veq_to_ssol = np.zeros((nyears, nrows, ncols))
grids_gdd_base0_veq_to_ssol = np.zeros((nyears, nrows, ncols))
grids_tmin_frz_veq_to_ssol = np.zeros((nyears, nrows, ncols))
grids_tavg_frz_veq_to_ssol = np.zeros((nyears, nrows, ncols))
grids_tmax_frz_veq_to_ssol = np.zeros((nyears, nrows, ncols))
#
grids_cd_ssol_to_aeq = np.zeros((nyears, nrows, ncols))
grids_cdd_ssol_to_aeq = np.zeros((nyears, nrows, ncols))
grids_gdd_veq_to_aeq = np.zeros((nyears, nrows, ncols))
grids_gdd_base0_veq_to_aeq = np.zeros((nyears, nrows, ncols))
grids_gdd_ssol_to_aeq = np.zeros((nyears, nrows, ncols))
grids_gdd_base0_ssol_to_aeq = np.zeros((nyears, nrows, ncols))
grids_tmin_frz_ssol_to_aeq = np.zeros((nyears, nrows, ncols))
grids_tavg_frz_ssol_to_aeq = np.zeros((nyears, nrows, ncols))
grids_tmax_frz_ssol_to_aeq = np.zeros((nyears, nrows, ncols))
#
grids_cd_aeq_to_wsol = np.zeros((nyears, nrows, ncols))
grids_cdd_aeq_to_wsol = np.zeros((nyears, nrows, ncols))
grids_gdd_aeq_to_wsol = np.zeros((nyears, nrows, ncols))
grids_gdd_base0_aeq_to_wsol = np.zeros((nyears, nrows, ncols))
grids_tmin_frz_aeq_to_wsol = np.zeros((nyears, nrows, ncols))
grids_tavg_frz_aeq_to_wsol = np.zeros((nyears, nrows, ncols))
grids_tmax_frz_aeq_to_wsol = np.zeros((nyears, nrows, ncols))
#
grids_prcp_365d_at_eoy = np.zeros((nyears, nrows, ncols))
#
for j, year in enumerate(years):
    message('processing Tavg, FD, CD, CDD, GDD and precip grids for %d' % year)
    #
    # winter grids ending at spring equinox (VEQ)
    doy_vequinox = date_to_doy(year, date_vequinox)
    filename = '%s/%d%d_NCEI_grids_2.h5' % \
        (path, year, str(doy_vequinox).zfill(4))
    message('- vernal equinox (%s)' % filename)
    grids_at_veq = getgrids(filename, nrows, ncols, inputgridvars)
    for i, var in enumerate(inputgridvars):
        inputgrids_at_veq[i, j, :, :] = grids_at_veq[i, :, :]
        message('-- %d %s mean %.1f' %
                (year, var, np.mean(inputgrids_at_veq[i, j, :, :])))
    #
    # spring grids ending at summer solstice (SSOL)
    doy_ssolstice = date_to_doy(year, date_ssolstice)
    filename = '%s/%d%d_NCEI_grids_2.h5' % \
        (path, year, str(doy_ssolstice).zfill(4))
    message('- summer solstice (%s)' % filename)
    grids_at_ssol = getgrids(filename, nrows, ncols, inputgridvars)
    for i, var in enumerate(inputgridvars):
        inputgrids_at_ssol[i, j, :, :] = grids_at_ssol[i, :, :]
        message('-- %d %s mean %.1f' %
                (year, var, np.mean(inputgrids_at_ssol[i, j, :, :])))
    #
    # summer grids ending at autumnal equinox (AEQ)
    doy_aequinox = date_to_doy(year, date_aequinox)
    filename = '%s/%d%d_NCEI_grids_2.h5' % \
        (path, year, str(doy_aequinox).zfill(4))
    message('- autumnal equinox (%s)' % filename)
    grids_at_aeq = getgrids(filename, nrows, ncols, inputgridvars)
    for i, var in enumerate(inputgridvars):
        inputgrids_at_aeq[i, j, :, :] = grids_at_aeq[i, :, :]
        message('-- %d %s mean %.1f' %
                (year, var, np.mean(inputgrids_at_aeq[i, j, :, :])))
    #
    # autumn grids ending at winter solstice (WSOL)
    doy_wsolstice = date_to_doy(year, date_wsolstice)
    filename = '%s/%d%d_NCEI_grids_2.h5' % \
        (path, year, str(doy_wsolstice).zfill(4))
    message('- winter solstice (%s)' % filename)
    grids_at_wsol = getgrids(filename, nrows, ncols, inputgridvars)
    for i, var in enumerate(inputgridvars):
        inputgrids_at_wsol[i, j, :, :] = grids_at_wsol[i, :, :]
        message('-- %d %s mean %.1f' %
                (year, var, np.mean(inputgrids_at_wsol[i, j, :, :])))
    #
    message('- spring seasonal accumulations')
    grids_intensity_winter[j, :, :] = \
        inputgrids_at_ssol[1, j, :, :] / inputgrids_at_ssol[0, j, :, :]
    message('-- %d intensity_winter mean %.1f' %
            (year, np.mean(grids_intensity_winter[j, :, :])))
    grids_cd_veq_to_ssol[j, :, :] = \
        inputgrids_at_ssol[0, j, :, :] - inputgrids_at_veq[0, j, :, :]
    message('-- %d cd_veq_to_ssol mean %.1f' %
            (year, np.mean(grids_cd_veq_to_ssol[j, :, :])))
    grids_cdd_veq_to_ssol[j, :, :] = \
        inputgrids_at_ssol[1, j, :, :] - inputgrids_at_veq[1, j, :, :]
    message('-- %d cdd_veq_to_ssol mean %.1f' %
            (year, np.mean(grids_cdd_veq_to_ssol[j, :, :])))
    grids_gdd_veq_to_ssol[j, :, :] = \
        inputgrids_at_ssol[2, j, :, :] - inputgrids_at_veq[2, j, :, :]
    message('-- %d gdd_veq_to_ssol mean %.1f' %
            (year, np.mean(grids_gdd_veq_to_ssol[j, :, :])))
    grids_tmin_frz_veq_to_ssol[j, :, :] = \
        inputgrids_at_ssol[3, j, :, :] - inputgrids_at_veq[3, j, :, :]
    message('-- %d tmin_frz_veq_to_ssol mean %.1f' %
            (year, np.mean(grids_tmin_frz_veq_to_ssol[j, :, :])))
    grids_tavg_frz_veq_to_ssol[j, :, :] = \
        inputgrids_at_ssol[4, j, :, :] - inputgrids_at_veq[4, j, :, :]
    message('-- %d tavg_frz_veq_to_ssol mean %.1f' %
            (year, np.mean(grids_tavg_frz_veq_to_ssol[j, :, :])))
    grids_tmax_frz_veq_to_ssol[j, :, :] = \
        inputgrids_at_ssol[5, j, :, :] - inputgrids_at_veq[5, j, :, :]
    message('-- %d tmax_frz_veq_to_ssol mean %.1f' %
            (year, np.mean(grids_tmax_frz_veq_to_ssol[j, :, :])))
    #
    message('- summer seasonal accumulations')
    grids_cd_ssol_to_aeq[j, :, :] = \
        inputgrids_at_aeq[0, j, :, :] - inputgrids_at_ssol[0, j, :, :]
    message('-- %d cd_ssol_to_aeq mean %.1f' %
            (year, np.mean(grids_cd_ssol_to_aeq[j, :, :])))
    grids_cdd_ssol_to_aeq[j, :, :] = \
        inputgrids_at_aeq[1, j, :, :] - inputgrids_at_ssol[1, j, :, :]
    message('-- %d cdd_ssol_to_aeq mean %.1f' %
            (year, np.mean(grids_cdd_ssol_to_aeq[j, :, :])))
    grids_gdd_veq_to_aeq[j, :, :] = \
        inputgrids_at_aeq[2, j, :, :] - inputgrids_at_veq[2, j, :, :]
    message('-- %d gdd_veq_to_aeq mean %.1f' %
            (year, np.mean(grids_gdd_veq_to_aeq[j, :, :])))
    grids_gdd_ssol_to_aeq[j, :, :] = \
        inputgrids_at_aeq[2, j, :, :] - inputgrids_at_ssol[2, j, :, :]
    message('-- %d gdd_ssol_to_aeq mean %.1f' %
            (year, np.mean(grids_gdd_ssol_to_aeq[j, :, :])))
    grids_tmin_frz_ssol_to_aeq[j, :, :] = \
        inputgrids_at_aeq[3, j, :, :] - inputgrids_at_ssol[3, j, :, :]
    message('-- %d tmin_frz_ssol_to_aeq mean %.1f' %
            (year, np.mean(grids_tmin_frz_ssol_to_aeq[j, :, :])))
    grids_tavg_frz_ssol_to_aeq[j, :, :] = \
        inputgrids_at_aeq[4, j, :, :] - inputgrids_at_ssol[4, j, :, :]
    message('-- %d tavg_frz_ssol_to_aeq mean %.1f' %
            (year, np.mean(grids_tavg_frz_ssol_to_aeq[j, :, :])))
    grids_tmax_frz_ssol_to_aeq[j, :, :] = \
        inputgrids_at_aeq[5, j, :, :] - inputgrids_at_ssol[5, j, :, :]
    message('-- %d tmax_frz_ssol_to_aeq mean %.1f' %
            (year, np.mean(grids_tmax_frz_ssol_to_aeq[j, :, :])))
    #
    message('- autumn seasonal accumulations')
    grids_cd_aeq_to_wsol[j, :, :] = \
        inputgrids_at_wsol[0, j, :, :] - inputgrids_at_aeq[0, j, :, :]
    message('-- %d cd_aeq_to_wsol mean %.1f' %
            (year, np.mean(grids_cd_aeq_to_wsol[j, :, :])))
    grids_cdd_aeq_to_wsol[j, :, :] = \
        inputgrids_at_wsol[1, j, :, :] - inputgrids_at_aeq[1, j, :, :]
    message('-- %d cdd_aeq_to_wsol mean %.1f' %
            (year, np.mean(grids_cdd_aeq_to_wsol[j, :, :])))
    grids_gdd_aeq_to_wsol[j, :, :] = \
        inputgrids_at_wsol[2, j, :, :] - inputgrids_at_aeq[2, j, :, :]
    message('-- %d gdd_aeq_to_wsol mean %.1f' %
            (year, np.mean(grids_gdd_aeq_to_wsol[j, :, :])))
    grids_tmin_frz_aeq_to_wsol[j, :, :] = \
        inputgrids_at_wsol[3, j, :, :] - inputgrids_at_aeq[3, j, :, :]
    message('-- %d tmin_frz_aeq_to_wsol mean %.1f' %
            (year, np.mean(grids_tmin_frz_aeq_to_wsol[j, :, :])))
    grids_tavg_frz_aeq_to_wsol[j, :, :] = \
        inputgrids_at_wsol[4, j, :, :] - inputgrids_at_aeq[4, j, :, :]
    message('-- %d tavg_frz_aeq_to_wsol mean %.1f' %
            (year, np.mean(grids_tavg_frz_aeq_to_wsol[j, :, :])))
    grids_tmax_frz_aeq_to_wsol[j, :, :] = \
        inputgrids_at_wsol[5, j, :, :] - inputgrids_at_aeq[5, j, :, :]
    message('-- %d tmax_frz_aeq_to_wsol mean %.1f' %
            (year, np.mean(grids_tmax_frz_aeq_to_wsol[j, :, :])))
    #
    # full-year precip ending with calendar year (EOY)
    doy_eoy = date_to_doy(year, date_eoy)
    filename = '%s/%d%d_NCEI_grids_2.h5' % (path, year, str(doy_eoy).zfill(4))
    message('- end of year (%s)' % filename)
    grids_prcp_365d_at_eoy[j, :, :] = getpgrid(filename, nrows, ncols)
    message('-- %d prcp_365d_at_eoy mean %.1f' %
            (year, np.mean(grids_prcp_365d_at_eoy[j, :, :])))
    message(' ')
#
print 'writing derived climatological grids to %s' % outfile
with hdf.File(outfile, 'r+') as h5outfile:
    for i in range(0, ninputgridvars):
        varname = 'grids_%s_at_veq' % inputgridvars[i]
        vargrids = inputgrids_at_veq[i, :, :, :]
        write_to_file(h5outfile, varname, vargrids)
    message('-- %d %d-year climatological grids at VEQ saved' %
            (ninputgridvars, nyears))
    #
    for i in range(0, ninputgridvars):
        varname = 'grids_%s_at_ssol' % inputgridvars[i]
        vargrids = inputgrids_at_ssol[i, :, :, :]
        write_to_file(h5outfile, varname, vargrids)
    message('-- %d %d-year climatological grids at SSOL saved' %
            (ninputgridvars, nyears))
    #
    for i in range(0, ninputgridvars):
        varname = 'grids_%s_at_aeq' % inputgridvars[i]
        vargrids = inputgrids_at_aeq[i, :, :, :]
        write_to_file(h5outfile, varname, vargrids)
    message('-- %d %d-year climatological grids at AEQ saved' %
            (ninputgridvars, nyears))
    #
    for i in range(0, ninputgridvars):
        varname = 'grids_%s_at_wsol' % inputgridvars[i]
        vargrids = inputgrids_at_wsol[i, :, :, :]
        write_to_file(h5outfile, varname, vargrids)
    message('-- %d %d-year climatological grids at WSOL saved' %
            (ninputgridvars, nyears))
    #
    write_to_file(h5outfile, 'grids_intensity_winter', grids_intensity_winter)
    message('-- 1 %d-year climatological winter intensity grids saved' %
            nyears)
    #
    write_to_file(h5outfile, 'grids_cd_veq_to_ssol', grids_cd_veq_to_ssol)
    write_to_file(h5outfile, 'grids_cdd_veq_to_ssol', grids_cdd_veq_to_ssol)
    write_to_file(h5outfile, 'grids_gdd_veq_to_ssol', grids_gdd_veq_to_ssol)
    write_to_file(h5outfile, 'grids_tmin_frz_veq_to_ssol',
                  grids_tmin_frz_veq_to_ssol)
    write_to_file(h5outfile, 'grids_tavg_frz_veq_to_ssol',
                  grids_tavg_frz_veq_to_ssol)
    write_to_file(h5outfile, 'grids_tmax_frz_veq_to_ssol',
                  grids_tmax_frz_veq_to_ssol)
    message('-- 6 %d-year climatological difference grids across spring \
            saved' % nyears)
    #
    write_to_file(h5outfile, 'grids_cd_ssol_to_aeq', grids_cd_ssol_to_aeq)
    write_to_file(h5outfile, 'grids_cdd_ssol_to_aeq', grids_cdd_ssol_to_aeq)
    write_to_file(h5outfile, 'grids_gdd_veq_to_aeq', grids_gdd_veq_to_aeq)
    write_to_file(h5outfile, 'grids_gdd_ssol_to_aeq', grids_gdd_ssol_to_aeq)
    write_to_file(h5outfile, 'grids_tmin_frz_ssol_to_aeq',
                  grids_tmin_frz_ssol_to_aeq)
    write_to_file(h5outfile, 'grids_tavg_frz_ssol_to_aeq',
                  grids_tavg_frz_ssol_to_aeq)
    write_to_file(h5outfile, 'grids_tmax_frz_ssol_to_aeq',
                  grids_tmax_frz_ssol_to_aeq)
    message('-- 7 %d-year climatological difference grids across summer \
            saved' % nyears)
    #
    write_to_file(h5outfile, 'grids_cd_aeq_to_wsol', grids_cd_aeq_to_wsol)
    write_to_file(h5outfile, 'grids_cdd_aeq_to_wsol', grids_cdd_aeq_to_wsol)
    write_to_file(h5outfile, 'grids_gdd_aeq_to_wsol', grids_gdd_aeq_to_wsol)
    write_to_file(h5outfile, 'grids_tmin_frz_aeq_to_wsol',
                  grids_tmin_frz_aeq_to_wsol)
    write_to_file(h5outfile, 'grids_tavg_frz_aeq_to_wsol',
                  grids_tavg_frz_aeq_to_wsol)
    write_to_file(h5outfile, 'grids_tmax_frz_aeq_to_wsol',
                  grids_tmax_frz_aeq_to_wsol)
    message('-- 6 %d-year climatological difference grids across autumn \
            saved' % nyears)
    #
    write_to_file(h5outfile, 'grids_prcp_365d_at_eoy', grids_prcp_365d_at_eoy)
    message('-- 1 %d-year climatological precipitation total grids saved' %
            nyears)
message(' ')
#
# establish grids to analyze frost-free growing season
grids_doy_last_spring_tmin_frz = np.zeros((nyears, nrows, ncols))
grids_cd_last_spring_tmin_frz = np.zeros((nyears, nrows, ncols))
grids_cdd_last_spring_tmin_frz = np.zeros((nyears, nrows, ncols))
grids_gdd_last_spring_tmin_frz = np.zeros((nyears, nrows, ncols))
grids_doy_first_autumn_tmin_frz = np.zeros((nyears, nrows, ncols))
grids_cd_first_autumn_tmin_frz = np.zeros((nyears, nrows, ncols))
grids_cdd_first_autumn_tmin_frz = np.zeros((nyears, nrows, ncols))
grids_gdd_first_autumn_tmin_frz = np.zeros((nyears, nrows, ncols))
grids_frost_free_season_days = np.zeros((nyears, nrows, ncols))
grids_frost_free_season_gdd = np.zeros((nyears, nrows, ncols))
#
# establish grids to analyze CD plateau growing season
grids_doy_plateau_begin = np.zeros((nyears, nrows, ncols))
grids_gdd_plateau_begin = np.zeros((nyears, nrows, ncols))
grids_doy_plateau_end = np.zeros((nyears, nrows, ncols))
grids_gdd_plateau_end = np.zeros((nyears, nrows, ncols))
grids_days_plateau_length = np.zeros((nyears, nrows, ncols))
grids_gdd_plateau_length = np.zeros((nyears, nrows, ncols))
grids_gdd_veq_to_plateau = np.zeros((nyears, nrows, ncols))
grids_gdd_aeq_to_plateau_end = np.zeros((nyears, nrows, ncols))
grids_gdd_plateau_to_wsol = np.zeros((nyears, nrows, ncols))
grids_intensity_plateau = np.zeros((nyears, nrows, ncols))
#
for j in range(0, nyears):
    year = years[j]
    doy_vequinox = date_to_doy(year, date_vequinox)
    doy_ssolstice = date_to_doy(year, date_ssolstice)
    doy_aequinox = date_to_doy(year, date_aequinox)
    doy_eoy = date_to_doy(year, date_eoy)
    doy_wsolstice = date_to_doy(year, date_wsolstice)
    message('processing %d for seasonal indicators' % year)
    #
    message('- last spring freeze and plateau begin dates')
    doy_last_spring_tmin_frz = np.zeros((nrows, ncols))
    cd_last_spring_tmin_frz = np.zeros((nrows, ncols))
    cdd_last_spring_tmin_frz = np.zeros((nrows, ncols))
    gdd_last_spring_tmin_frz = np.zeros((nrows, ncols))
    doy_plateau_begin = np.zeros((nrows, ncols))
    gdd_plateau_begin = np.zeros((nrows, ncols))
    last_grid_tmin_frz = np.zeros((nrows, ncols))
    last_grid_cd = np.zeros((nrows, ncols))
    # search through a specific 90-day window
    for doy in range(doy_vequinox, doy_ssolstice):
        filename = '%s/%d%d_NCEI_grids_2.h5' % (path, year, str(doy).zfill(4))
        grids_at_doy = getgrids(filename, nrows, ncols, inputgridvars, 5)
        grid_cd = grids_at_doy[0, :, :]
        grid_cdd = grids_at_doy[1, :, :]
        grid_gdd = grids_at_doy[2, :, :]
        grid_tmin_frz = grids_at_doy[4, :, :]
        doy_last_spring_tmin_frz = \
            np.where(grid_tmin_frz != last_grid_tmin_frz, float(doy),
                     doy_last_spring_tmin_frz)
        cd_last_spring_tmin_frz = \
            np.where(doy_last_spring_tmin_frz == doy, grid_cd,
                     cd_last_spring_tmin_frz)
        cdd_last_spring_tmin_frz = \
            np.where(doy_last_spring_tmin_frz == doy, grid_cdd,
                     cdd_last_spring_tmin_frz)
        gdd_last_spring_tmin_frz = \
            np.where(doy_last_spring_tmin_frz == doy, grid_gdd,
                     gdd_last_spring_tmin_frz)
        doy_plateau_begin = \
            np.where(grid_cd != last_grid_cd, float(doy), doy_plateau_begin)
        gdd_plateau_begin = \
            np.where(doy_plateau_begin == doy, grid_gdd, gdd_plateau_begin)
        last_grid_tmin_frz = grid_tmin_frz
        last_grid_cd = grid_cd
    grids_doy_last_spring_tmin_frz[j, :, :] = doy_last_spring_tmin_frz[:, :]
    message('-- %d doy_last_spring_tmin_frz mean %.1f' %
            (year, np.mean(grids_doy_last_spring_tmin_frz[j, :, :])))
    grids_cd_last_spring_tmin_frz[j, :, :] = cd_last_spring_tmin_frz[:, :]
    message('-- %d cd_last_spring_tmin_frz mean %.1f' %
            (year, np.mean(grids_cd_last_spring_tmin_frz[j, :, :])))
    grids_cdd_last_spring_tmin_frz[j, :, :] = cdd_last_spring_tmin_frz[:, :]
    message('-- %d cdd_last_spring_tmin_frz mean %.1f' %
            (year, np.mean(grids_cdd_last_spring_tmin_frz[j, :, :])))
    grids_gdd_last_spring_tmin_frz[j, :, :] = gdd_last_spring_tmin_frz[:, :]
    message('-- %d gdd_last_spring_tmin_frz mean %.1f' %
            (year, np.mean(grids_gdd_last_spring_tmin_frz[j, :, :])))
    grids_doy_plateau_begin[j, :, :] = doy_plateau_begin[:, :]
    message('-- %d doy_plateau_begin mean %.1f' %
            (year, np.mean(grids_doy_plateau_begin[j, :, :])))
    grids_gdd_plateau_begin[j, :, :] = gdd_plateau_begin[:, :]
    message('-- %d gdd_plateau_begin mean %.1f' %
            (year, np.mean(grids_gdd_plateau_begin[j, :, :])))
    #
    message('- first autumn freeze and plateau end dates')
    doy_first_autumn_tmin_frz = np.zeros((nrows, ncols))
    cd_first_autumn_tmin_frz = np.zeros((nrows, ncols))
    cdd_first_autumn_tmin_frz = np.zeros((nrows, ncols))
    gdd_first_autumn_tmin_frz = np.zeros((nrows, ncols))
    doy_plateau_end = np.zeros((nrows, ncols))
    gdd_plateau_end = np.zeros((nrows, ncols))
    last_grid_cd = np.zeros((nrows, ncols))
    last_grid_tmin_frz = np.zeros((nrows, ncols))
    # search through a specific 90-day window
    for doy in range(doy_wsolstice - 30, doy_aequinox - 30, -1):
        filename = '%s/%d%d_NCEI_grids_2.h5' % (path, year, str(doy).zfill(4))
        grids_at_doy = getgrids(filename, nrows, ncols, inputgridvars, 5)
        grid_cd = grids_at_doy[0, :, :]
        grid_cdd = grids_at_doy[1, :, :]
        grid_gdd = grids_at_doy[2, :, :]
        grid_tmin_frz = grids_at_doy[4, :, :]
        doy_first_autumn_tmin_frz = \
            np.where(grid_tmin_frz != last_grid_tmin_frz, float(doy + 1),
                     doy_first_autumn_tmin_frz)
        cd_first_autumn_tmin_frz = \
            np.where(doy_first_autumn_tmin_frz == doy + 1, grid_cd,
                     cd_first_autumn_tmin_frz)
        cdd_first_autumn_tmin_frz = \
            np.where(doy_first_autumn_tmin_frz == doy + 1, grid_cdd,
                     cdd_first_autumn_tmin_frz)
        gdd_first_autumn_tmin_frz = \
            np.where(doy_first_autumn_tmin_frz == doy + 1, grid_gdd,
                     gdd_first_autumn_tmin_frz)
        doy_plateau_end = \
            np.where(grid_cd != last_grid_cd, float(doy), doy_plateau_end)
        gdd_plateau_end = \
            np.where(doy_plateau_end == doy, grid_gdd, gdd_plateau_end)
        last_grid_cd = grid_cd
        last_grid_tmin_frz = grid_tmin_frz
    grids_doy_first_autumn_tmin_frz[j, :, :] = doy_first_autumn_tmin_frz[:, :]
    message('-- %d doy_first_autumn_tmin_frz mean %.1f' %
            (year, np.mean(grids_doy_first_autumn_tmin_frz[j, :, :])))
    grids_cd_first_autumn_tmin_frz[j, :, :] = cd_first_autumn_tmin_frz[:, :]
    message('-- %d cd_first_autumn_tmin_frz mean %.1f' %
            (year, np.mean(grids_cd_first_autumn_tmin_frz[j, :, :])))
    grids_cdd_first_autumn_tmin_frz[j, :, :] = cdd_first_autumn_tmin_frz[:, :]
    message('-- %d cdd_first_autumn_tmin_frz mean %.1f' %
            (year, np.mean(grids_cdd_first_autumn_tmin_frz[j, :, :])))
    grids_gdd_first_autumn_tmin_frz[j, :, :] = gdd_first_autumn_tmin_frz[:, :]
    message('-- %d gdd_first_autumn_tmin_frz mean %.1f' %
            (year, np.mean(grids_gdd_first_autumn_tmin_frz[j, :, :])))
    grids_doy_plateau_end[j, :, :] = doy_plateau_end[:, :]
    message('-- %d doy_plateau_end mean %.1f' %
            (year, np.mean(grids_doy_plateau_end[j, :, :])))
    grids_gdd_plateau_end[j, :, :] = gdd_plateau_end[:, :]
    message('-- %d gdd_plateau_end mean %.1f' %
            (year, np.mean(grids_gdd_plateau_end[j, :, :])))
    #
    message('- frost-free and plateau-based growing seasons')
    grids_frost_free_season_days[j, :, :] = \
        grids_doy_first_autumn_tmin_frz[j, :, :] \
        - grids_doy_last_spring_tmin_frz[j, :, :]
    message('-- %d frost_free_season_days mean %.1f' %
            (year, np.mean(grids_frost_free_season_days[j, :, :])))
    grids_frost_free_season_gdd[j, :, :] = \
        grids_gdd_first_autumn_tmin_frz[j, :, :] \
        - grids_gdd_last_spring_tmin_frz[j, :, :]
    message('-- %d frost_free_season_gdd mean %.1f' %
            (year, np.mean(grids_frost_free_season_gdd[j, :, :])))
    grids_days_plateau_length[j, :, :] = \
        grids_doy_plateau_end[j, :, :] - grids_doy_plateau_begin[j, :, :]
    message('-- %d days_plateau_length mean %.1f' %
            (year, np.mean(grids_days_plateau_length[j, :, :])))
    grids_gdd_plateau_length[j, :, :] = \
        grids_gdd_plateau_end[j, :, :] - grids_gdd_plateau_begin[j, :, :]
    message('-- %d gdd_plateau_length mean %.1f' %
            (year, np.mean(grids_gdd_plateau_length[j, :, :])))
    grids_gdd_veq_to_plateau[j, :, :] = \
        grids_gdd_plateau_begin[j, :, :] - inputgrids_at_veq[2, j, :, :]
    message('-- %d gdd_veq_to_plateau mean %.1f' %
            (year, np.mean(grids_gdd_veq_to_plateau[j, :, :])))
    grids_gdd_aeq_to_plateau_end[j, :, :] = \
        grids_gdd_plateau_end[j, :, :] - inputgrids_at_aeq[2, j, :, :]
    message('-- %d gdd_aeq_to_plateau_end mean %.1f' %
            (year, np.mean(grids_gdd_aeq_to_plateau_end[j, :, :])))
    grids_gdd_plateau_to_wsol[j, :, :] = \
        inputgrids_at_wsol[2, j, :, :] - grids_gdd_plateau_end[j, :, :]
    message('-- %d gdd_plateau_to_wsol mean %.1f' %
            (year, np.mean(grids_gdd_plateau_to_wsol[j, :, :])))
    grids_intensity_plateau[j, :, :] = \
        grids_gdd_plateau_length[j, :, :] / grids_days_plateau_length[j, :, :]
    message('-- %d intensity_plateau mean %.1f' %
            (year, np.mean(grids_intensity_plateau[j, :, :])))
    message(' ')
#
message('writing growing season indicator grids to %s' % outfile)
with hdf.File(outfile, 'r+') as h5outfile:
    write_to_file(h5outfile, 'grids_doy_last_spring_tmin_frz',
                  grids_doy_last_spring_tmin_frz)
    write_to_file(h5outfile, 'grids_cd_last_spring_tmin_frz',
                  grids_cd_last_spring_tmin_frz)
    write_to_file(h5outfile, 'grids_cdd_last_spring_tmin_frz',
                  grids_cdd_last_spring_tmin_frz)
    write_to_file(h5outfile, 'grids_gdd_last_spring_tmin_frz',
                  grids_gdd_last_spring_tmin_frz)
    write_to_file(h5outfile, 'grids_doy_first_autumn_tmin_frz',
                  grids_doy_first_autumn_tmin_frz)
    write_to_file(h5outfile, 'grids_cd_first_autumn_tmin_frz',
                  grids_cd_first_autumn_tmin_frz)
    write_to_file(h5outfile, 'grids_cdd_first_autumn_tmin_frz',
                  grids_cdd_first_autumn_tmin_frz)
    write_to_file(h5outfile, 'grids_gdd_first_autumn_tmin_frz',
                  grids_gdd_first_autumn_tmin_frz)
    write_to_file(h5outfile, 'grids_frost_free_season_days',
                  grids_frost_free_season_days)
    write_to_file(h5outfile, 'grids_frost_free_season_gdd',
                  grids_frost_free_season_gdd)
    message('-- 10 collections of %d-year frost-based climatological grids \
            saved' % nyears)
    #
    write_to_file(h5outfile, 'grids_doy_plateau_begin',
                  grids_doy_plateau_begin)
    write_to_file(h5outfile, 'grids_gdd_plateau_begin',
                  grids_gdd_plateau_begin)
    write_to_file(h5outfile, 'grids_doy_plateau_end',
                  grids_doy_plateau_end)
    write_to_file(h5outfile, 'grids_gdd_plateau_end',
                  grids_gdd_plateau_end)
    write_to_file(h5outfile, 'grids_days_plateau_length',
                  grids_days_plateau_length)
    write_to_file(h5outfile, 'grids_gdd_plateau_length',
                  grids_gdd_plateau_length)
    write_to_file(h5outfile, 'grids_gdd_veq_to_plateau',
                  grids_gdd_veq_to_plateau)
    write_to_file(h5outfile, 'grids_gdd_aeq_to_plateau_end',
                  grids_gdd_aeq_to_plateau_end)
    write_to_file(h5outfile, 'grids_gdd_plateau_to_wsol',
                  grids_gdd_plateau_to_wsol)
    write_to_file(h5outfile, 'grids_intensity_plateau',
                  grids_intensity_plateau)
    message('-- 10 collections of %d-year plateau-based climatological grids \
            saved' % nyears)
message(' ')
#
message('process_NCEI_04a.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_04a.py
