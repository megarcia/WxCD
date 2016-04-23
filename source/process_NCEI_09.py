"""
Python script 'process_NCEI_09.py'
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

USAGE: 'python process_NCEI_09.py 1984 2013 ./analyses 1 1'

PURPOSE: Generate maps of numerous variables on annual and summary bases

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for input
              The 'Plots' module has its own requirements

INPUT: Aggregated grid datacubes and statistics grids in '.h5' file from 'process_NCEI_04.py'
       (with the naming convention 'analyses/[YYYY]-[YYYY]_derived_clim_grids.h5')
       Ecoregion maps and grid info in 'clipped_ecoregions.h5' from 'process_NCEI_08.py'

OUTPUT: Lots and lots of maps in '.png' format in your 'analyses/X_maps' folders

RUN TIME: ~75 sec per data year
"""

import os, sys, datetime
import h5py as hdf
import numpy as np
from Plots import masked_map_plot_geo

def message(str):
    print str
    sys.stdout.flush()
    return    

message(' ')
message('process_NCEI_09.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 6:
    annual_plots = 0
    summary_plots = 1
else:
    annual_plots = int(sys.argv[4])
    summary_plots = int(sys.argv[5])
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
years = np.arange(year_begin,year_end+1).astype(int)
nyears = len(years)
#
# get working area size/shape/location from clipped ecoregion maps file
infile = '%s/../data/clipped_ecoregions.h5' % path
message('extracting grid information and land mask from ecoregion maps file %s' % infile)
h5infile = hdf.File(infile,'r') 
UTM_zone = np.copy(h5infile['grid/UTM_zone'])
clip_bounds = np.copy(h5infile['grid/clip_bounds'])
UTM_bounds = clip_bounds[:4]
landmask = np.copy(h5infile['landmask_reduced'])
h5infile.close()
message(' ')
#
infile = '%s/%d-%d_derived_clim_grids.h5' % (path, year_begin, year_end)
if annual_plots:
    message('getting gridded annual fields from %s' % infile)
    h5infile = hdf.File(infile,'r') 
    grids_tavg_90d_avg_at_veq = np.copy(h5infile['grids_tavg_90d_avg_at_veq'])
    grids_tavg_90d_var_at_veq = np.copy(h5infile['grids_tavg_90d_var_at_veq'])
    grids_prcp_90d_sum_at_veq = np.copy(h5infile['grids_prcp_90d_sum_at_veq'])
    grids_prcp_90d_nd0_at_veq = np.copy(h5infile['grids_prcp_90d_nd0_sum_at_veq'])
    grids_prcp_90d_nd10_at_veq = np.copy(h5infile['grids_prcp_90d_nd10_sum_at_veq'])
    grids_prcp_90d_nd25_at_veq = np.copy(h5infile['grids_prcp_90d_nd25_sum_at_veq'])
    grids_tavg_90d_avg_at_ssol = np.copy(h5infile['grids_tavg_90d_avg_at_ssol'])
    grids_tavg_90d_var_at_ssol = np.copy(h5infile['grids_tavg_90d_var_at_ssol'])
    grids_prcp_90d_sum_at_ssol = np.copy(h5infile['grids_prcp_90d_sum_at_ssol'])
    grids_prcp_90d_nd0_at_ssol = np.copy(h5infile['grids_prcp_90d_nd0_sum_at_ssol'])
    grids_prcp_90d_nd10_at_ssol = np.copy(h5infile['grids_prcp_90d_nd10_sum_at_ssol'])
    grids_prcp_90d_nd25_at_ssol = np.copy(h5infile['grids_prcp_90d_nd25_sum_at_ssol'])
    grids_tavg_90d_avg_at_aeq = np.copy(h5infile['grids_tavg_90d_avg_at_aeq'])
    grids_tavg_90d_var_at_aeq = np.copy(h5infile['grids_tavg_90d_var_at_aeq'])
    grids_prcp_90d_sum_at_aeq = np.copy(h5infile['grids_prcp_90d_sum_at_aeq'])
    grids_prcp_90d_nd0_at_aeq = np.copy(h5infile['grids_prcp_90d_nd0_sum_at_aeq'])
    grids_prcp_90d_nd10_at_aeq = np.copy(h5infile['grids_prcp_90d_nd10_sum_at_aeq'])
    grids_prcp_90d_nd25_at_aeq = np.copy(h5infile['grids_prcp_90d_nd25_sum_at_aeq'])
    grids_tavg_90d_avg_at_wsol = np.copy(h5infile['grids_tavg_90d_avg_at_wsol'])
    grids_tavg_90d_var_at_wsol = np.copy(h5infile['grids_tavg_90d_var_at_wsol'])
    grids_prcp_90d_sum_at_wsol = np.copy(h5infile['grids_prcp_90d_sum_at_wsol'])
    grids_prcp_90d_nd0_at_wsol = np.copy(h5infile['grids_prcp_90d_nd0_sum_at_wsol'])
    grids_prcp_90d_nd10_at_wsol = np.copy(h5infile['grids_prcp_90d_nd10_sum_at_wsol'])
    grids_prcp_90d_nd25_at_wsol = np.copy(h5infile['grids_prcp_90d_nd25_sum_at_wsol'])
    grids_chill_d_at_ssol = np.copy(h5infile['grids_chill_d_at_ssol'])
    grids_tmin_frz_days_at_ssol = np.copy(h5infile['grids_tmin_frz_days_at_ssol'])
    grids_intensity_winter = np.copy(h5infile['grids_intensity_winter'])
    grids_gdd_veq_to_aeq = np.copy(h5infile['grids_gdd_veq_to_aeq'])
    grids_prcp_365d_at_eoy = np.copy(h5infile['grids_prcp_365d_at_eoy'])
    grids_doy_last_spring_tmin_frz = np.copy(h5infile['grids_doy_last_spring_tmin_frz'])
    grids_gdd_last_spring_tmin_frz = np.copy(h5infile['grids_gdd_last_spring_tmin_frz'])
    grids_doy_plateau_begin = np.copy(h5infile['grids_doy_plateau_begin'])
    grids_gdd_plateau_begin = np.copy(h5infile['grids_gdd_plateau_begin'])
    grids_doy_first_autumn_tmin_frz = np.copy(h5infile['grids_doy_first_autumn_tmin_frz'])
    grids_doy_plateau_end = np.copy(h5infile['grids_doy_plateau_end'])
    grids_gdd_plateau_end = np.copy(h5infile['grids_gdd_plateau_end'])
    grids_frost_free_season_days = np.copy(h5infile['grids_frost_free_season_days'])
    grids_days_plateau_length = np.copy(h5infile['grids_days_plateau_length'])
    grids_gdd_plateau_length = np.copy(h5infile['grids_gdd_plateau_length'])
    grids_intensity_plateau = np.copy(h5infile['grids_intensity_plateau'])
    h5infile.close()
    message(' ')
    #
    for j in range(0,nyears):
        year = years[j]
        message('plotting selected annual grids for %d' % year)
        #
        titlestr = '%d winter Tavg mean' % year
        fname = '%s/annual_maps/%d_tavg_90d_avg_at_veq_map.png' % (path, year)
        masked_map_plot_geo(grids_tavg_90d_avg_at_veq[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d winter Tavg variance' % year
        fname = '%s/annual_maps/%d_tavg_90d_var_at_veq_map.png' % (path, year)
        masked_map_plot_geo(grids_tavg_90d_var_at_veq[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d winter Precip (cm)' % year
        fname = '%s/annual_maps/%d_prcp_90d_sum_at_veq_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_sum_at_veq[j,:,:],landmask,UTM_zone,UTM_bounds,'YlGnBu','tight',titlestr,fname)
        titlestr = '%d winter Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd0_at_veq_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd0_at_veq[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d winter moderate Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd10_at_veq_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd10_at_veq[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d winter heavy Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd25_at_veq_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd25_at_veq[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        #
        titlestr = '%d spring Tavg mean' % year
        fname = '%s/annual_maps/%d_tavg_90d_avg_at_ssol_map.png' % (path, year)
        masked_map_plot_geo(grids_tavg_90d_avg_at_ssol[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d spring Tavg variance' % year
        fname = '%s/annual_maps/%d_tavg_90d_var_at_ssol_map.png' % (path, year)
        masked_map_plot_geo(grids_tavg_90d_var_at_ssol[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d spring Precip (cm)' % year
        fname = '%s/annual_maps/%d_prcp_90d_sum_at_ssol_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_sum_at_ssol[j,:,:],landmask,UTM_zone,UTM_bounds,'YlGnBu','tight',titlestr,fname)
        titlestr = '%d spring Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd0_at_ssol_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd0_at_ssol[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d spring moderate Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd10_at_ssol_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd10_at_ssol[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d spring heavy Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd25_at_ssol_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd25_at_ssol[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        #
        titlestr = '%d summer Tavg mean' % year
        fname = '%s/annual_maps/%d_tavg_90d_avg_at_aeq_map.png' % (path, year)
        masked_map_plot_geo(grids_tavg_90d_avg_at_aeq[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d summer Tavg variance' % year
        fname = '%s/annual_maps/%d_tavg_90d_var_at_aeq_map.png' % (path, year)
        masked_map_plot_geo(grids_tavg_90d_var_at_aeq[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d summer Precip (cm)' % year
        fname = '%s/annual_maps/%d_prcp_90d_sum_at_aeq_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_sum_at_aeq[j,:,:],landmask,UTM_zone,UTM_bounds,'YlGnBu','tight',titlestr,fname)
        titlestr = '%d summer Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd0_at_aeq_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd0_at_aeq[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d summer moderate Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd10_at_aeq_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd10_at_aeq[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d summer heavy Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd25_at_aeq_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd25_at_aeq[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        #
        titlestr = '%d autumn Tavg mean' % year
        fname = '%s/annual_maps/%d_tavg_90d_avg_at_wsol_map.png' % (path, year)
        masked_map_plot_geo(grids_tavg_90d_avg_at_wsol[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d autumn Tavg variance' % year
        fname = '%s/annual_maps/%d_tavg_90d_var_at_wsol_map.png' % (path, year)
        masked_map_plot_geo(grids_tavg_90d_var_at_wsol[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d autumn Precip (cm)' % year
        fname = '%s/annual_maps/%d_prcp_90d_sum_at_wsol_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_sum_at_wsol[j,:,:],landmask,UTM_zone,UTM_bounds,'YlGnBu','tight',titlestr,fname)
        titlestr = '%d autumn Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd0_at_wsol_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd0_at_wsol[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d autumn moderate Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd10_at_wsol_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd10_at_wsol[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d autumn heavy Precip days' % year
        fname = '%s/annual_maps/%d_prcp_90d_nd25_at_wsol_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_90d_nd25_at_wsol[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        #
        titlestr = '%d plateau CD' % year
        fname = '%s/annual_maps/%d_plateau_CD_map.png' % (path, year)
        masked_map_plot_geo(grids_chill_d_at_ssol[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d winter freezing nights' % year
        fname = '%s/annual_maps/%d_winter_tmin_frz_map.png' % (path, year)
        masked_map_plot_geo(grids_tmin_frz_days_at_ssol[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d winter intensity (CDD/CD)' % year
        fname = '%s/annual_maps/%d_winter_intensity_map.png' % (path, year)
        masked_map_plot_geo(grids_intensity_winter[j,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname)
        titlestr = '%d last spring Tmin freeze DOY' % year
        fname = '%s/annual_maps/%d_DOY_last_spring_tmin_freeze_map.png' % (path, year)
        masked_map_plot_geo(grids_doy_last_spring_tmin_frz[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','tight',titlestr,fname)
        titlestr = '%d last spring Tmin freeze GDD' % year
        fname = '%s/annual_maps/%d_gdd_last_spring_tmin_freeze_map.png' % (path, year)
        masked_map_plot_geo(grids_gdd_last_spring_tmin_frz[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','tight',titlestr,fname)
        #
        titlestr = '%d precipitation (cm)' % year
        fname = '%s/annual_maps/%d_prcp_365d_at_eoy_map.png' % (path, year)
        masked_map_plot_geo(grids_prcp_365d_at_eoy[j,:,:],landmask,UTM_zone,UTM_bounds,'YlGnBu','tight',titlestr,fname)
        #
        titlestr = '%d plateau beginning DOY' % year
        fname = '%s/annual_maps/%d_DOY_plateau_begin_map.png' % (path, year)
        masked_map_plot_geo(grids_doy_plateau_begin[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','tight',titlestr,fname)
        titlestr = '%d plateau beginning GDD' % year
        fname = '%s/annual_maps/%d_gdd_plateau_begin_map.png' % (path, year)
        masked_map_plot_geo(grids_gdd_plateau_begin[j,:,:],landmask,UTM_zone,UTM_bounds,'Reds','tight',titlestr,fname)
        titlestr = '%d last spring Tmin freeze DOY' % year
        fname = '%s/annual_maps/%d_DOY_last_spring_tmin_freeze_map.png' % (path, year)
        masked_map_plot_geo(grids_doy_last_spring_tmin_frz[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','tight',titlestr,fname)
        titlestr = '%d last spring Tmin freeze GDD' % year
        fname = '%s/annual_maps/%d_gdd_last_spring_tmin_freeze_map.png' % (path, year)
        masked_map_plot_geo(grids_gdd_last_spring_tmin_frz[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d spring + summer GDD' % year
        fname = '%s/annual_maps/%d_veq_to_aeq_GDD_map.png' % (path, year)
        masked_map_plot_geo(grids_gdd_veq_to_aeq[j,:,:],landmask,UTM_zone,UTM_bounds,'Reds','tight',titlestr,fname)
        titlestr = '%d first autumn Tmin freeze DOY' % year
        fname = '%s/annual_maps/%d_DOY_first_autumn_tmin_freeze_map.png' % (path, year)
        masked_map_plot_geo(grids_doy_first_autumn_tmin_frz[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d plateau ending DOY' % year
        fname = '%s/annual_maps/%d_DOY_plateau_end_map.png' % (path, year)
        masked_map_plot_geo(grids_doy_plateau_end[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d plateau ending GDD' % year
        fname = '%s/annual_maps/%d_gdd_plateau_end_map.png' % (path, year)
        masked_map_plot_geo(grids_gdd_plateau_end[j,:,:],landmask,UTM_zone,UTM_bounds,'Reds','tight',titlestr,fname)
        titlestr = '%d frost-free season length (days)' % year
        fname = '%s/annual_maps/%d_frost_free_season_days_map.png' % (path, year)
        masked_map_plot_geo(grids_frost_free_season_days[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d plateau length (days)' % year
        fname = '%s/annual_maps/%d_days_plateau_length_map.png' % (path, year)
        masked_map_plot_geo(grids_days_plateau_length[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        titlestr = '%d plateau GDD' % year
        fname = '%s/annual_maps/%d_gdd_plateau_length_map.png' % (path, year)
        masked_map_plot_geo(grids_gdd_plateau_length[j,:,:],landmask,UTM_zone,UTM_bounds,'Reds','tight',titlestr,fname)
        titlestr = '%d growing season intensity (GDD/day)' % year
        fname = '%s/annual_maps/%d_plateau_intensity_map.png' % (path, year)
        masked_map_plot_geo(grids_intensity_plateau[j,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname)
        message(' ')
#
if summary_plots:
    message('getting gridded summary statistics fields from %s' % infile)
    h5infile = hdf.File(infile,'r') 
    stats_tavg_90d_avg_at_veq = np.copy(h5infile['stats_tavg_90d_avg_at_veq'])
    stats_tavg_90d_var_at_veq = np.copy(h5infile['stats_tavg_90d_var_at_veq'])
    stats_prcp_90d_sum_at_veq = np.copy(h5infile['stats_prcp_90d_sum_at_veq'])
    stats_prcp_90d_nd0_at_veq = np.copy(h5infile['stats_prcp_90d_nd0_sum_at_veq'])
    stats_prcp_90d_nd10_at_veq = np.copy(h5infile['stats_prcp_90d_nd10_sum_at_veq'])
    stats_prcp_90d_nd25_at_veq = np.copy(h5infile['stats_prcp_90d_nd25_sum_at_veq'])
    stats_tavg_90d_avg_at_ssol = np.copy(h5infile['stats_tavg_90d_avg_at_ssol'])
    stats_tavg_90d_var_at_ssol = np.copy(h5infile['stats_tavg_90d_var_at_ssol'])
    stats_prcp_90d_sum_at_ssol = np.copy(h5infile['stats_prcp_90d_sum_at_ssol'])
    stats_prcp_90d_nd0_at_ssol = np.copy(h5infile['stats_prcp_90d_nd0_sum_at_ssol'])
    stats_prcp_90d_nd10_at_ssol = np.copy(h5infile['stats_prcp_90d_nd10_sum_at_ssol'])
    stats_prcp_90d_nd25_at_ssol = np.copy(h5infile['stats_prcp_90d_nd25_sum_at_ssol'])
    stats_tavg_90d_avg_at_aeq = np.copy(h5infile['stats_tavg_90d_avg_at_aeq'])
    stats_tavg_90d_var_at_aeq = np.copy(h5infile['stats_tavg_90d_var_at_aeq'])
    stats_prcp_90d_sum_at_aeq = np.copy(h5infile['stats_prcp_90d_sum_at_aeq'])
    stats_prcp_90d_nd0_at_aeq = np.copy(h5infile['stats_prcp_90d_nd0_sum_at_aeq'])
    stats_prcp_90d_nd10_at_aeq = np.copy(h5infile['stats_prcp_90d_nd10_sum_at_aeq'])
    stats_prcp_90d_nd25_at_aeq = np.copy(h5infile['stats_prcp_90d_nd25_sum_at_aeq'])
    stats_tavg_90d_avg_at_wsol = np.copy(h5infile['stats_tavg_90d_avg_at_wsol'])
    stats_tavg_90d_var_at_wsol = np.copy(h5infile['stats_tavg_90d_var_at_wsol'])
    stats_prcp_90d_sum_at_wsol = np.copy(h5infile['stats_prcp_90d_sum_at_wsol'])
    stats_prcp_90d_nd0_at_wsol = np.copy(h5infile['stats_prcp_90d_nd0_sum_at_wsol'])
    stats_prcp_90d_nd10_at_wsol = np.copy(h5infile['stats_prcp_90d_nd10_sum_at_wsol'])
    stats_prcp_90d_nd25_at_wsol = np.copy(h5infile['stats_prcp_90d_nd25_sum_at_wsol'])
    stats_chill_d_at_veq = np.copy(h5infile['stats_chill_d_at_veq'])
    stats_chill_d_at_ssol = np.copy(h5infile['stats_chill_d_at_ssol'])
    stats_intensity_winter = np.copy(h5infile['stats_intensity_winter'])
    stats_prcp_365d_at_eoy = np.copy(h5infile['stats_prcp_365d_at_eoy'])
    stats_doy_last_spring_tmin_frz = np.copy(h5infile['stats_doy_last_spring_tmin_frz'])
    stats_gdd_last_spring_tmin_frz = np.copy(h5infile['stats_gdd_last_spring_tmin_frz'])
    stats_doy_first_autumn_tmin_frz = np.copy(h5infile['stats_doy_first_autumn_tmin_frz'])
    stats_frost_free_season_days = np.copy(h5infile['stats_frost_free_season_days'])
    stats_doy_plateau_begin = np.copy(h5infile['stats_doy_plateau_begin'])
    stats_gdd_plateau_begin = np.copy(h5infile['stats_gdd_plateau_begin'])
    stats_doy_plateau_end = np.copy(h5infile['stats_doy_plateau_end'])
    stats_gdd_plateau_end = np.copy(h5infile['stats_gdd_plateau_end'])
    stats_days_plateau_length = np.copy(h5infile['stats_days_plateau_length'])
    stats_gdd_plateau_length = np.copy(h5infile['stats_gdd_plateau_length'])
    stats_intensity_plateau = np.copy(h5infile['stats_intensity_plateau'])
    h5infile.close()
    message(' ')
    #
    message('plotting selected %d-year summary statistics grids' % nyears)
    titlestr = '%d-%d winter Tavg mean' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_avg_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_avg_at_veq[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d winter Tavg mean trend' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_avg_trend_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_avg_at_veq[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d winter Tavg variance' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_var_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_var_at_veq[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d winter Precip (cm)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_avg_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_sum_at_veq[0,:,:],landmask,UTM_zone,UTM_bounds,'YlGnBu','tight',titlestr,fname,1,1)
    titlestr = '%d-%d winter Precip trend' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_avg_trend_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_sum_at_veq[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d winter Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd0_avg_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd0_at_veq[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d winter moderate Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd10_avg_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd10_at_veq[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d winter heavy Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd25_avg_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd25_at_veq[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    #
    titlestr = '%d-%d spring Tavg mean' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_avg_at_ssol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_avg_at_ssol[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d spring Tavg mean trend' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_avg_trend_at_ssol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_avg_at_ssol[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d spring Tavg variance' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_var_at_ssol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_var_at_ssol[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d spring Precip (cm)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_avg_at_ssol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_sum_at_ssol[0,:,:],landmask,UTM_zone,UTM_bounds,'YlGnBu','tight',titlestr,fname,1,1)
    titlestr = '%d-%d spring Precip trend' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_avg_trend_at_ssol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_sum_at_ssol[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d spring Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd0_avg_at_ssol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd0_at_ssol[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d spring moderate Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd10_avg_at_ssol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd10_at_ssol[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d spring heavy Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd25_avg_at_ssol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd25_at_ssol[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    #
    titlestr = '%d-%d summer Tavg mean' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_avg_at_aeq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_avg_at_aeq[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d summer Tavg mean trend' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_avg_trend_at_aeq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_avg_at_aeq[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d summer Tavg variance' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_var_at_aeq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_var_at_aeq[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d summer Precip (cm)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_avg_at_aeq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_sum_at_aeq[0,:,:],landmask,UTM_zone,UTM_bounds,'YlGnBu','tight',titlestr,fname,1,1)
    titlestr = '%d-%d summer Precip trend' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_avg_trend_at_aeq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_sum_at_aeq[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d summer Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd0_avg_at_aeq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd0_at_aeq[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d summer moderate Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd10_avg_at_aeq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd10_at_aeq[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d summer heavy Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd25_avg_at_aeq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd25_at_aeq[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    #
    titlestr = '%d-%d autumn Tavg mean' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_avg_at_wsol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_avg_at_wsol[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d autumn Tavg mean trend' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_avg_trend_at_wsol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_avg_at_wsol[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d autumn Tavg variance' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_tavg_90d_var_at_wsol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_tavg_90d_var_at_wsol[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d autumn Precip (cm)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_avg_at_wsol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_sum_at_wsol[0,:,:],landmask,UTM_zone,UTM_bounds,'YlGnBu','tight',titlestr,fname,1,1)
    titlestr = '%d-%d autumn Precip trend' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_avg_trend_at_wsol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_sum_at_wsol[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d autumn Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd0_avg_at_wsol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd0_at_wsol[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d autumn moderate Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd10_avg_at_wsol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd10_at_wsol[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d autumn heavy Precip days' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_90d_nd25_avg_at_wsol_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_90d_nd25_at_wsol[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    #
    titlestr = '%d-%d mean CD at spring equinox' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_chill_d_avg_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_chill_d_at_veq[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of CD at spring equinox' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_chill_d_std_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_chill_d_at_veq[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of CD at spring equinox' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_chill_d_trend_at_veq_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_chill_d_at_veq[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean plateau CD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_plateau_chill_d_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_chill_d_at_ssol[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of plateau CD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_plateau_chill_d_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_chill_d_at_ssol[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of plateau CD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_plateau_chill_d_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_chill_d_at_ssol[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean winter intensity (CDD/CD)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_winter_intensity_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_intensity_winter[0,:,:],landmask,UTM_zone,UTM_bounds,'Blues','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of winter intensity (CDD/CD)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_winter_intensity_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_intensity_winter[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of winter intensity (CDD/CD)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_winter_intensity_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_intensity_winter[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','balanced',titlestr,fname,1,1)
    #
    titlestr = '%d-%d mean annual precipitation (cm)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_365d_avg_at_eoy_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_365d_at_eoy[0,:,:],landmask,UTM_zone,UTM_bounds,'YlGnBu','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of annual precipitation (cm)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_365d_std_at_eoy_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_365d_at_eoy[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of annual precipitation (cm)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_prcp_365d_trend_at_eoy_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_prcp_365d_at_eoy[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','balanced',titlestr,fname,1,1)
    #
    titlestr = '%d-%d mean DOY of last spring Tmin freeze' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_last_spring_tmin_frz_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_last_spring_tmin_frz[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of DOY of last spring Tmin freeze (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_last_spring_tmin_frz_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_last_spring_tmin_frz[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of DOY of last spring Tmin freeze (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_last_spring_tmin_frz_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_last_spring_tmin_frz[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow_r','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean GDD of last spring Tmin freeze' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_last_spring_tmin_frz_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_last_spring_tmin_frz[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of GDD of last spring Tmin freeze' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_last_spring_tmin_frz_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_last_spring_tmin_frz[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of GDD of last spring Tmin freeze' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_last_spring_tmin_frz_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_last_spring_tmin_frz[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean DOY of first autumn Tmin freeze' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_first_autumn_tmin_frz_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_first_autumn_tmin_frz[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of DOY of first autumn Tmin freeze (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_first_autumn_tmin_frz_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_first_autumn_tmin_frz[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of DOY of first autumn Tmin freeze (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_first_autumn_tmin_frz_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_first_autumn_tmin_frz[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean frost-free season (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_frost_free_season_days_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_frost_free_season_days[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of frost-free season (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_frost_free_season_days_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_frost_free_season_days[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of frost-free season (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_frost_free_season_days_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_frost_free_season_days[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean plateau beginning DOY' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_plateau_begin_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_plateau_begin[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of plateau beginning DOY (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_plateau_begin_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_plateau_begin[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of plateau beginning DOY (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_plateau_begin_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_plateau_begin[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean plateau beginning GDD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_plateau_begin_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_plateau_begin[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of plateau beginning GDD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_plateau_begin_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_plateau_begin[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of plateau beginning GDD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_plateau_begin_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_plateau_begin[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean plateau ending DOY' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_plateau_end_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_plateau_end[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of plateau ending DOY (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_plateau_end_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_plateau_end[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of plateau ending DOY (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_doy_plateau_end_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_doy_plateau_end[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean plateau ending GDD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_plateau_end_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_plateau_end[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of plateau ending GDD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_plateau_end_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_plateau_end[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of plateau ending GDD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_plateau_end_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_plateau_end[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean plateau length (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_days_plateau_length_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_days_plateau_length[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of plateau length (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_days_plateau_length_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_days_plateau_length[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of plateau length (days)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_days_plateau_length_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_days_plateau_length[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean plateau GDD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_plateau_length_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_plateau_length[0,:,:],landmask,UTM_zone,UTM_bounds,'Reds','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of plateau GDD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_plateau_length_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_plateau_length[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of plateau GDD' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_gdd_plateau_length_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_gdd_plateau_length[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
    titlestr = '%d-%d mean plateau intensity (GDD/day)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_plateau_intensity_avg_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_intensity_plateau[0,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','tight',titlestr,fname,1,1)
    titlestr = '%d-%d stdev of plateau intensity (GDD/day)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_plateau_intensity_std_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_intensity_plateau[1,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','zero',titlestr,fname,1,1)
    titlestr = '%d-%d trend of plateau intensity (GDD/day)' % (year_begin,year_end)
    fname = '%s/summary_maps/%d-%d_plateau_intensity_trend_map.png' % (path,year_begin,year_end)
    masked_map_plot_geo(stats_intensity_plateau[4,:,:],landmask,UTM_zone,UTM_bounds,'rainbow','balanced',titlestr,fname,1,1)
message(' ')
#
message('process_NCEI_09.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_09.py
