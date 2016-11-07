"""
Python script 'process_NCEI_07.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Plotting meteorological variables on daily/seasonal/annual/decadal
         basis

NOTE: These are special plots (not maps) for publication figures

DEPENDENCIES: h5py, numpy
              'Date_Convert' module has no external requirements
              'Plots' module has its own requirements

USAGE: '$ python process_NCEI_07.py 1984 2013 ./analyses'


INPUT: '.csv' output files from process_NCEI_05.py
       (with the naming convention
        'analyses/[YYYY]-[YYYY]_[XXXX]_season[al]_stats.csv')
       '.h5' output file from process_NCEI_06.py
       (with the naming convention
        'analyses/[YYYY]-[YYYY]_clim_values_by_doy.h5')

OUTPUT: '.png' in specified output directory
"""


import sys
import datetime
import h5py as hdf
import numpy as np
from Plots import p_ts_daily_plot, t_ts_daily_plot, ttt_ts_daily_plot, \
    shifted_cd_plot, gdd_plot, p_ts_annual_plot, t_ts_bydecade_plot, \
    p_ts_bydecade_plot, p_ts_daily_var_plot


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


message(' ')
message('process_NCEI_07.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input warning: no input directory indicated, using ./analyses')
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
years = np.arange(year_begin, year_end + 1).astype(int)
nyears = len(years)
#
doy_begin = 1
doy_end = 365
doys = np.arange(doy_begin, doy_end + 1).astype(int)
ndoys = len(doys)
#
h5infname = '%s/%d-%d_clim_values_by_doy.h5' % (path, year_begin, year_end)
message('getting series information from %s' % h5infname)
with hdf.File(h5infname, 'r') as h5infile:
    tmin_by_year_doy = np.copy(h5infile['tmin_by_year_doy'])
    tavg_by_year_doy = np.copy(h5infile['tavg_by_year_doy'])
    tmax_by_year_doy = np.copy(h5infile['tmax_by_year_doy'])
    prcp_by_year_doy = np.copy(h5infile['prcp_by_year_doy'])
    cd_by_year_doy = np.copy(h5infile['cd_by_year_doy'])
    cdd_by_year_doy = np.copy(h5infile['cdd_by_year_doy'])
    gdd_by_year_doy = np.copy(h5infile['gdd_by_year_doy'])
    prcp_mean_by_doy = np.copy(h5infile['prcp_mean_by_doy'])
    prcp_min_by_doy = np.copy(h5infile['prcp_min_by_doy'])
    prcp_max_by_doy = np.copy(h5infile['prcp_max_by_doy'])
    tmin_mean_by_doy = np.copy(h5infile['tmin_mean_by_doy'])
    tmin_min_by_doy = np.copy(h5infile['tmin_min_by_doy'])
    tmin_max_by_doy = np.copy(h5infile['tmin_max_by_doy'])
    tmax_mean_by_doy = np.copy(h5infile['tmax_mean_by_doy'])
    tmax_min_by_doy = np.copy(h5infile['tmax_min_by_doy'])
    tmax_max_by_doy = np.copy(h5infile['tmax_max_by_doy'])
    tavg_mean_by_doy = np.copy(h5infile['tavg_mean_by_doy'])
    tavg_min_by_doy = np.copy(h5infile['tavg_min_by_doy'])
    tavg_max_by_doy = np.copy(h5infile['tavg_max_by_doy'])
    cd_mean_by_doy = np.copy(h5infile['cd_mean_by_doy'])
    cd_min_by_doy = np.copy(h5infile['cd_min_by_doy'])
    cd_max_by_doy = np.copy(h5infile['cd_max_by_doy'])
    gdd_mean_by_doy = np.copy(h5infile['gdd_mean_by_doy'])
    gdd_min_by_doy = np.copy(h5infile['gdd_min_by_doy'])
    gdd_max_by_doy = np.copy(h5infile['gdd_max_by_doy'])
message(' ')
#
# summary plots of 4 variables with mean, min, max traces
titlestr = '%d-%d study area average Precipitation (cm)' % \
    (year_begin, year_end)
filename = '%s/figures/%d-%d_Prcp.png' % (path, year_begin, year_end)
p_ts_daily_plot(prcp_mean_by_doy, prcp_min_by_doy, prcp_max_by_doy,
                titlestr, filename)
titlestr = '%d-%d study area average Tmin (%sC)' % \
    (year_begin, year_end, r'$^\circ$')
filename = '%s/figures/%d-%d_Tmin.png' % (path, year_begin, year_end)
t_ts_daily_plot(tmin_mean_by_doy, tmin_min_by_doy, tmin_max_by_doy,
                titlestr, filename, fill=0)
titlestr = '%d-%d study area average Tmax (%sC)' % \
    (year_begin, year_end, r'$^\circ$')
filename = '%s/figures/%d-%d_Tmax.png' % (path, year_begin, year_end)
t_ts_daily_plot(tmax_mean_by_doy, tmax_min_by_doy, tmax_max_by_doy,
                titlestr, filename, fill=0)
titlestr = '%d-%d study area average Tavg (%sC)' % \
    (year_begin, year_end, r'$^\circ$')
filename = '%s/figures/%d-%d_Tavg.png' % (path, year_begin, year_end)
t_ts_daily_plot(tavg_mean_by_doy, tavg_min_by_doy, tavg_max_by_doy,
                titlestr, filename, tbase=5.0)
#
# get some useful statistics from process_NCEI_05.py output
fname = '%s/%d-%d_overall_stats.csv' % (path, year_begin, year_end)
allstats = np.genfromtxt(fname, delimiter=',')
#
# summary plot of tavg_mean, tmin_mean, and tmax_mean traces
# (shaded, with tbase threshold marked)
# TO DO: replace hard-coded allstats row numbers
mean_plateau_begin = allstats[41, 0]
mean_plateau_end = allstats[43, 0]
mean_last_frost = allstats[35, 0]
mean_first_frost = allstats[38, 0]
titlestr = '%d-%d study area average Tmin, Tavg, Tmax (%sC)' % \
    (year_begin, year_end, r'$^\circ$')
filename = '%s/figures/%d-%d_Tmin_Tavg_Tmax.png' % (path, year_begin, year_end)
ttt_ts_daily_plot(tavg_mean_by_doy, tmin_mean_by_doy, tmax_mean_by_doy,
                  mean_plateau_begin, mean_plateau_end, mean_last_frost,
                  mean_first_frost, titlestr, filename, tbase=5.0)
#
# summary plot of CD mean, min, max traces shifted with 1 Jul at beginning
# TO DO: replace hard-coded allstats row numbers
mean_plateau_cd = allstats[36, 0]
titlestr = '%d-%d study area average Chilling Days' % (year_begin, year_end)
filename = '%s/figures/%d-%d_CD.png' % (path, year_begin, year_end)
shifted_cd_plot(cd_mean_by_doy, cd_min_by_doy, cd_max_by_doy,
                mean_plateau_begin, mean_plateau_end, mean_plateau_cd,
                titlestr, filename)
#
# summary plot of GDD mean, min, max traces (no shift required)
titlestr = '%d-%d study area average Growing Degree Days' % \
    (year_begin, year_end)
filename = '%s/figures/%d-%d_GDD.png' % (path, year_begin, year_end)
gdd_plot(gdd_mean_by_doy, gdd_min_by_doy, gdd_max_by_doy, mean_plateau_begin,
         mean_plateau_end, titlestr, filename)
#
# get some useful precip statistics from process_NCEI_05.py output
fname = '%s/%d-%d_prcp_seasonal_stats.csv' % (path, year_begin, year_end)
prcp_seasonal = np.genfromtxt(fname, delimiter=',')
#
# summary plot of seasonal prcp
titlestr = '%d-%d study area average Precipitation (cm)' % \
    (year_begin, year_end)
filename = '%s/figures/%d-%d_Prcp_seasons.png' % (path, year_begin, year_end)
p_ts_annual_plot(prcp_seasonal[:, 0], prcp_seasonal[:, 8],
                 prcp_seasonal[:, 16], prcp_seasonal[:, 24],
                 prcp_seasonal[:, 32], prcp_seasonal[:, 33],
                 year_begin, year_end, titlestr, filename)
#
# decadal summaries
# TO DO: replace hard-coded 10-year segments with arbitrary segments
tmin_decade_1 = np.mean(tmin_by_year_doy[0:9, :], axis=0)
tmin_decade_2 = np.mean(tmin_by_year_doy[10:19, :], axis=0)
tmin_decade_3 = np.mean(tmin_by_year_doy[20:29, :], axis=0)
tavg_decade_1 = np.mean(tavg_by_year_doy[0:9, :], axis=0)
tavg_decade_2 = np.mean(tavg_by_year_doy[10:19, :], axis=0)
tavg_decade_3 = np.mean(tavg_by_year_doy[20:29, :], axis=0)
tmax_decade_1 = np.mean(tmax_by_year_doy[0:9, :], axis=0)
tmax_decade_2 = np.mean(tmax_by_year_doy[10:19, :], axis=0)
tmax_decade_3 = np.mean(tmax_by_year_doy[20:29, :], axis=0)
prcp_decade_1 = np.mean(prcp_by_year_doy[0:9, :], axis=0)
prcp_decade_2 = np.mean(prcp_by_year_doy[10:19, :], axis=0)
prcp_decade_3 = np.mean(prcp_by_year_doy[20:29, :], axis=0)
prcp_all_mean = np.mean(prcp_by_year_doy[:, :], axis=0)
prcp_all_std = np.std(prcp_by_year_doy[:, :], axis=0)
prcp_all_max = np.max(prcp_by_year_doy[:, :], axis=0)
#
titlestr = '%d-%d study area average Tmin (%sC) by decade' % \
    (year_begin, year_end, r'$^\circ$')
filename = '%s/figures/%d-%d_decadal_Tmin.png' % (path, year_begin, year_end)
t_ts_bydecade_plot(tmin_decade_1, tmin_decade_2, tmin_decade_3,
                   titlestr, filename)
titlestr = '%d-%d study area average Tavg (%sC) by decade' % \
    (year_begin, year_end, r'$^\circ$')
filename = '%s/figures/%d-%d_decadal_Tavg.png' % (path, year_begin, year_end)
t_ts_bydecade_plot(tavg_decade_1, tavg_decade_2, tavg_decade_3,
                   titlestr, filename)
titlestr = '%d-%d study area average Tmax (%sC) by decade' % \
    (year_begin, year_end, r'$^\circ$')
filename = '%s/figures/%d-%d_decadal_Tmax.png' % (path, year_begin, year_end)
t_ts_bydecade_plot(tmax_decade_1, tmax_decade_2, tmax_decade_3,
                   titlestr, filename)
titlestr = '%d-%d study area average Prcp (cm) by decade' % \
    (year_begin, year_end)
filename = '%s/figures/%d-%d_Prcp.png' % (path, year_begin, year_end)
p_ts_bydecade_plot(prcp_decade_1, prcp_decade_2, prcp_decade_3,
                   titlestr, filename)
titlestr = '%d-%d study area average Prcp (cm)' % (year_begin, year_end)
filename = '%s/figures/%d-%d_Prcp2.png' % (path, year_begin, year_end)
p_ts_daily_var_plot(prcp_all_mean, prcp_all_std, prcp_all_max,
                    titlestr, filename)
message(' ')
#
message('process_NCEI_07.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_07.py
