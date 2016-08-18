"""
Python module 'Plots.py'
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
               Garcia, M., and P.A. Townsend (in review): "Recent climatological trends and 
               potential influences on forest phenology around western Lake Superior, USA." 
               Submitted to J. Geophys. Res. Atmos. on 5 April 2016, revised 20 August 2016.
           See also 'README.md', 'CITATION.txt', and 'ACKNOWLEDGEMENTS.txt' for more information.

USAGE: insert 'from Plots import *' line near head of script, then call individual plotting
       routines as indicated

PURPOSE: Plotting T/P maps, daily/seasonal/annual time series, specially formatted series 

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The 'Date_Convert' module has no external requirements
              The 'UTM_Geo_Convert' module has its own requirements

INPUT: Provided by calling script

OUTPUT: New '.png' files in a location specified by variable <fname> in most cases

RUN TIME: very small (plot generation)
"""

import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from Date_Convert import *
from UTM_Geo_Convert import *

def message(str):
    print str
    sys.stdout.flush()
    return    

def p_map_plot(x,y,z,xmin,xmax,ymin,ymax,grid,UTMz,titlestr,fname,stations=1):
    plt.imshow(grid, extent=(xmin, xmax, ymin, ymax), origin='bottom', cmap=plt.cm.Blues)
    plt.clim(0,2.5)
    plt.colorbar()
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel('UTM %dN easting (m)' % UTMz)
    plt.yticks(rotation='vertical')
    plt.ylabel('UTM %dN northing (m)' % UTMz)
    plt.title(titlestr)
    if stations:
        plt.hold(True)
        plt.scatter(x, y, s=20, c=z, cmap=plt.cm.Blues)
        plt.clim(0,5)
        plt.hold(False)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()
    return

def t_map_plot(x,y,z,xmin,xmax,ymin,ymax,grid,UTMz,titlestr,fname,stations=1):
    plt.imshow(grid, extent=(xmin, xmax, ymin, ymax), origin='bottom', cmap=plt.cm.jet)
    plt.clim(-20,40)
    plt.colorbar()
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel('UTM %dN easting (m)' % UTMz)
    plt.yticks(rotation='vertical')
    plt.ylabel('UTM %dN northing (m)' % UTMz)
    plt.title(titlestr)
    if stations:
        plt.hold(True)
        plt.scatter(x, y, s=20, c=z, cmap=plt.cm.jet)
        plt.clim(-20,40)
        plt.hold(False)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()
    return

def masked_map_plot_utm(img,mask,UTM_zone,UTM_bounds,cmapname,cbar,titlestr,fname,notitle=1):
    UTM_W = UTM_bounds[0] / 1000.0
    UTM_N = UTM_bounds[1] / 1000.0
    UTM_E = UTM_bounds[2] / 1000.0
    UTM_S = UTM_bounds[3] / 1000.0    
    nrows, ncols = np.shape(img)
    stretch = (UTM_E - UTM_W) / (UTM_N - UTM_S)
    img_mod = np.where(mask == 1, img, np.nan)
    if cbar == 'tight':
        img_mod = np.where(img_mod == 0.0, np.nan, img_mod)
    plt.imshow(img_mod, extent=(UTM_W, UTM_E, UTM_S, UTM_N), aspect=stretch, origin='bottom', interpolation='nearest', cmap=plt.get_cmap(cmapname))
    if cbar == 'tight':
        img_mod = np.where(img_mod == 0.0, np.nan, img_mod)
        minval = np.floor(np.nanmin(img_mod))
        maxval = np.ceil(np.nanmax(img_mod))    
        plt.clim(minval, maxval)
    elif cbar == 'balanced':
        minval = np.nanmin(img_mod)
        maxval = np.nanmax(img_mod)
        bound = np.max([np.abs(minval),np.abs(maxval)])
        minval = -1 * bound
        maxval = bound
        plt.clim(minval, maxval)
    elif cbar == 'zero':
        minval = 0.0
        maxval = np.ceil(np.nanmax(img_mod))
        plt.clim(minval, maxval)
    else:
        pass
    plt.colorbar()
    plt.xlim([UTM_W, UTM_E])
    plt.ylim([UTM_S, UTM_N])
    plt.xlabel('UTM %dN easting (km)' % UTM_zone)
    plt.yticks(rotation='vertical')
    plt.ylabel('UTM %dN northing (km)' % UTM_zone)
    if notitle:
        fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
        plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()   
    return 

def masked_map_plot_geo(img,mask,UTM_zone,UTM_bounds,cmapname,cbar,titlestr,fname,filter=0,notitle=0):
    UTM_W = UTM_bounds[0]
    UTM_N = UTM_bounds[1]
    UTM_E = UTM_bounds[2]
    UTM_S = UTM_bounds[3]
    nrows, ncols = np.shape(img)
    stretch = (UTM_E - UTM_W) / (UTM_N - UTM_S)
    if filter:
        img = gaussian_filter(img, 3)
    img_mod = np.where(mask == 1, img, np.nan)
    if cbar == 'tight':
        img_mod = np.where(img_mod == 0.0, np.nan, img_mod)
    plt.imshow(img_mod, extent=(UTM_W, UTM_E, UTM_S, UTM_N), aspect=stretch, origin='bottom', interpolation='nearest', cmap=plt.get_cmap(cmapname))
    if cbar == 'tight':
        img_mod = np.where(img_mod == 0.0, np.nan, img_mod)
        minval = np.nanmin(img_mod)
        maxval = np.nanmax(img_mod)
        plt.clim(minval, maxval)
    elif cbar == 'balanced':
        minval = np.nanmin(img_mod)
        maxval = np.nanmax(img_mod)
        bound = np.max([np.abs(minval),np.abs(maxval)])
        minval = -1 * bound
        maxval = bound
        plt.clim(minval, maxval)
    elif cbar == 'zero':
        minval = 0.0
        maxval = np.ceil(np.nanmax(img_mod))
        plt.clim(minval, maxval)
    else:
        pass
    if cbar != 'none':
        plt.colorbar()
    plt.xlim([UTM_W, UTM_E])
    plt.ylim([UTM_S, UTM_N])
    LL_lon, LL_lat = utm_to_geographic(UTM_bounds[0], UTM_bounds[3], UTM_zone)
    UL_lon, UL_lat = utm_to_geographic(UTM_bounds[0], UTM_bounds[1], UTM_zone)
    UR_lon, UR_lat = utm_to_geographic(UTM_bounds[2], UTM_bounds[1], UTM_zone)
    LR_lon, LR_lat = utm_to_geographic(UTM_bounds[2], UTM_bounds[3], UTM_zone)
    # left_geog_ticks
    left_utm_ticks = []
    left_geog_ticks = np.arange(np.ceil(LL_lat),np.ceil(UL_lat))
    for lat in left_geog_ticks:
        lon = LL_lon + (((LL_lon - UL_lon) / (LL_lat - UL_lat)) * (lat - LL_lat))
        zone, east, north = geographic_to_utm(lon, lat, UTM_zone)
        left_utm_ticks.append(np.round(north,0))
    # bot_geog_ticks
    bot_utm_ticks = []
    bot_geog_ticks = np.arange(np.ceil(LL_lon),np.ceil(LR_lon))
    for lon in bot_geog_ticks:
        lat = LL_lat + (((LR_lat - LL_lat) / (LR_lon - LL_lon)) * (lon - LL_lon))
        zone, east, north = geographic_to_utm(lon, lat, UTM_zone)
        bot_utm_ticks.append(np.round(east,0))
    # print axis ticks
    plt.xticks(bot_utm_ticks,bot_geog_ticks)
    plt.yticks(left_utm_ticks,left_geog_ticks, rotation='vertical')
    if filter:
        fname = fname[:-4] + '_filter' + fname[-4:]
    if notitle:
        fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
        plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()   
    return 

def p_ts_daily_plot(avg_pts,min_pts,max_pts,titlestr,fname):
    doy_min = 1
    doy_max = 365
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy_first_dom = [1]
    for i in range(0,len(days_in_month)-1):
        doy_first_dom.append(doy_first_dom[i] + days_in_month[i])
    date_vequinox = 320; doy_vequinox = date_to_doy(2013,date_vequinox)
    date_ssolstice = 621; doy_ssolstice = date_to_doy(2013,date_ssolstice)
    date_aequinox = 922; doy_aequinox = date_to_doy(2013,date_aequinox)
    date_wsolstice = 1221; doy_wsolstice = date_to_doy(2013,date_wsolstice)
    ymin = 0.0
    ymax = 6.0
    plt.plot([doy_vequinox, doy_vequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_vequinox-1,4.0,'vernal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_ssolstice, doy_ssolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_ssolstice-1,4.0,'summer solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_aequinox, doy_aequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_aequinox-1,4.0,'autumnal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_wsolstice, doy_wsolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_wsolstice-1,4.0,'winter solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot(avg_pts, 'k-', linewidth=2)
    plt.plot(max_pts, 'b-', linewidth=2)
    plt.xlim([doy_min, doy_max])
    plt.xticks(doy_first_dom)
    plt.xlabel('Day of Year')
    plt.ylim([ymin, ymax])
    plt.yticks(rotation='vertical')
    plt.ylabel('Area-average Precipitation (cm)')
    fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
    plt.savefig(fname, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf() 
    return   

def p_ts_pentad_plot(avg_pts,titlestr,fname):
    doy_min = 1
    doy_max = 365
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy_first_dom = [1]
    for i in range(0,len(days_in_month)-1):
        doy_first_dom.append(doy_first_dom[i] + days_in_month[i])
    #
    # TO DO:
    #   accumulate daily precip to pentads
    #   plot pentad precip as column bar chart
    #
    plt.xlim([doy_min, doy_max])
    plt.xticks(doy_first_dom)
    plt.xlabel('Day of Year')
    ymin = 0.0
    ymax = 6.0
    plt.ylim([ymin, ymax])
    plt.yticks(rotation='vertical')
    plt.ylabel('Precipitation (cm)')
    fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
    plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf() 
    return   
    

def p_ts_daily_var_plot(mean_pts,std_vals,max_pts,titlestr,fname):
    doy_min = 1
    doy_max = 365
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy_first_dom = [1]
    for i in range(0,len(days_in_month)-1):
        doy_first_dom.append(doy_first_dom[i] + days_in_month[i])
    date_vequinox = 320; doy_vequinox = date_to_doy(2013,date_vequinox)
    date_ssolstice = 621; doy_ssolstice = date_to_doy(2013,date_ssolstice)
    date_aequinox = 922; doy_aequinox = date_to_doy(2013,date_aequinox)
    date_wsolstice = 1221; doy_wsolstice = date_to_doy(2013,date_wsolstice)
    ymin = 0.0
    ymax = 1.5
    zero_pts = np.zeros(np.shape(mean_pts))
    plt.plot([doy_vequinox, doy_vequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_vequinox-1,1.0,'vernal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_ssolstice, doy_ssolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_ssolstice-1,1.0,'summer solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_aequinox, doy_aequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_aequinox-1,1.0,'autumnal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_wsolstice, doy_wsolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_wsolstice-1,1.0,'winter solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    mean_plus_pts = mean_pts + std_vals
    mean_minus_pts = mean_pts - std_vals
    plt.plot(mean_pts, 'k-', linewidth=2)
    x = np.linspace(doy_min-1,doy_max,doy_max)
    plt.fill_between(x, mean_plus_pts, mean_minus_pts, color='k', alpha=0.3)
    # plt.plot(x, max_pts, 'k', marker='o', markeredgecolor='k', markersize=2, linewidth=0)
    plt.xlim([doy_min, doy_max])
    plt.xticks(doy_first_dom)
    plt.xlabel('Day of Year')
    plt.ylim([ymin, ymax])
    plt.yticks(rotation='vertical')
    plt.ylabel('Precipitation (cm)')
    fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
    plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()    
    return

def p_ts_annual_plot(winter,spring,summer,autumn,annual,annual_std,year_min,year_max,titlestr,fname):
    years = np.arange(year_min,year_max+1).astype(int)
    ymin = 0.0
    ymax = 120.0
    plt.errorbar(years, annual, yerr=annual_std, color='k', linestyle='-', marker='o', markeredgecolor='k', markersize=5, linewidth=2, elinewidth=1, label='Annual')
    plt.plot(years, winter, color='b', linestyle='-', marker='o', markeredgecolor='b', markersize=5, linewidth=2, label='Winter')
    plt.plot(years, spring, color='g', linestyle='-', marker='o', markeredgecolor='g', markersize=5, linewidth=2, label='Spring')
    plt.plot(years, summer, color='r', linestyle='-', marker='o', markeredgecolor='r', markersize=5, linewidth=2, label='Summer')
    plt.plot(years, autumn, color='y', linestyle='-', marker='o', markeredgecolor='y', markersize=5, linewidth=2, label='Autumn')
    plt.legend(fontsize=8, loc='upper right', bbox_to_anchor=(0.87, 0.99))
    plt.xlim([year_min-1, year_max+1])
    plt.xlabel('Year')
    plt.ylim([ymin, ymax])
    plt.yticks(rotation='vertical')
    plt.ylabel('Area-average Precipitation (cm)')
    fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
    plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()    
    return

def t_ts_daily_plot(avg_pts,min_pts,max_pts,titlestr,fname,tbase=0.0,fill=1):
    doy_min = 1
    doy_max = 365
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy_first_dom = [1]
    for i in range(0,len(days_in_month)-1):
        doy_first_dom.append(doy_first_dom[i] + days_in_month[i])
    date_vequinox = 320; doy_vequinox = date_to_doy(2013,date_vequinox)
    date_ssolstice = 621; doy_ssolstice = date_to_doy(2013,date_ssolstice)
    date_aequinox = 922; doy_aequinox = date_to_doy(2013,date_aequinox)
    date_wsolstice = 1221; doy_wsolstice = date_to_doy(2013,date_wsolstice)
    ymin = -45.0
    ymax = 40.0
    zero_pts = np.zeros(np.shape(avg_pts))
    plt.plot(zero_pts, 'k-', linewidth=1)
    plt.plot([doy_vequinox, doy_vequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_vequinox-1,20.0,'vernal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_ssolstice, doy_ssolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_ssolstice-1,-40.0,'summer solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_aequinox, doy_aequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_aequinox-1,-40.0,'autumnal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_wsolstice, doy_wsolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_wsolstice-1,20.0,'winter solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    if tbase > 0.0:
        tbase_pts = np.full_like(avg_pts, fill_value=tbase)
        if fill:
            x = np.linspace(doy_min-1,doy_max,doy_max)
            plt.fill_between(x, tbase_pts, avg_pts, where=(avg_pts > tbase), color='r', alpha=0.4)
            plt.fill_between(x, tbase_pts, avg_pts, where=(avg_pts < tbase), color='b', alpha=0.3)
            plt.fill_between(x, zero_pts, avg_pts, where=(avg_pts < 0.0), color='b', alpha=0.5)
        plt.plot(tbase_pts, 'k--', linewidth=1)
        pltstr = '%s %d%sC' % (r'$T_{base} = $',int(tbase),r'$^\circ$')
        plt.text(5,tbase+0.5,pltstr,fontsize=10,ha='left',va='bottom')
    plt.plot(avg_pts, 'k-', linewidth=2)
    plt.plot(min_pts, 'b-', linewidth=2)
    plt.plot(max_pts, 'r-', linewidth=2)
    plt.xlim([doy_min, doy_max])
    plt.xticks(doy_first_dom)
    plt.xlabel('Day of Year')
    plt.ylim([ymin, ymax])
    plt.yticks(rotation='vertical')
    plt.ylabel('Temperature (%sC)' % r'$^\circ$')
    fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
    plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()
    return    

def ttt_ts_daily_plot(avg_pts,min_pts,max_pts,m_plateau_beg,m_plateau_end,m_last_frost,m_first_frost,titlestr,fname,tbase=0.0,fill=1):
    doy_min = 1
    doy_max = 365
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy_first_dom = [1]
    for i in range(0,len(days_in_month)-1):
        doy_first_dom.append(doy_first_dom[i] + days_in_month[i])
    date_vequinox = 320; doy_vequinox = date_to_doy(2013,date_vequinox)
    date_ssolstice = 621; doy_ssolstice = date_to_doy(2013,date_ssolstice)
    date_aequinox = 922; doy_aequinox = date_to_doy(2013,date_aequinox)
    date_wsolstice = 1221; doy_wsolstice = date_to_doy(2013,date_wsolstice)
    ymin = -25.0
    ymax = 30.0
    zero_pts = np.zeros(np.shape(avg_pts))
    plt.plot(zero_pts, 'k-', linewidth=1)
    plt.plot([doy_vequinox, doy_vequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_vequinox-1,15.0,'vernal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_ssolstice, doy_ssolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_ssolstice-1,-23.0,'summer solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_aequinox, doy_aequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_aequinox-1,-23.0,'autumnal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_wsolstice, doy_wsolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_wsolstice-1,15.0,'winter solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    if tbase > 0.0:
        tbase_pts = np.full_like(avg_pts, fill_value=tbase)
        if fill:
            x = np.linspace(doy_min-1,doy_max,doy_max)
            plt.fill_between(x, tbase_pts, avg_pts, where=(avg_pts > tbase), color='r', alpha=0.5)
            plt.fill_between(x, tbase_pts, avg_pts, where=(avg_pts < tbase), color='b', alpha=0.3)
            plt.fill_between(x, zero_pts, avg_pts, where=(avg_pts < 0.0), color='b', alpha=0.5)
        plt.plot(tbase_pts, 'k--', linewidth=1)
        pltstr = '%s %d%sC' % (r'$T_{base} = $',int(tbase),r'$^\circ$')
        plt.text(10,tbase+0.5,pltstr,fontsize=10,ha='left',va='bottom')
    plt.plot(max_pts, 'r-', linewidth=2, label=r'$T_{max}$')
    plt.plot(avg_pts, 'k-', linewidth=2, label=r'$T_{avg}$')
    plt.plot(min_pts, 'b-', linewidth=2, label=r'$T_{min}$')
    plt.legend(fontsize=10, loc='upper left')
    plateau_begin = m_plateau_beg
    plateau_end = m_plateau_end
    plt.plot([plateau_begin, plateau_begin],[-2, 21], 'k-', linewidth=1)
    plt.text(plateau_begin,21.25,'mean plateau\nbegin date',fontsize=8,ha='center',va='bottom')
    plt.plot([plateau_end, plateau_end],[-2, 21], 'k-', linewidth=1)
    plt.text(plateau_end,21.25,'mean plateau\nend date',fontsize=8,ha='center',va='bottom')
    last_frost = m_last_frost
    first_frost = m_first_frost
    plt.plot([last_frost, last_frost],[-5, 7], 'k-', linewidth=1)
    plt.text(last_frost,-5.25,'mean last\nfrost date',fontsize=8,ha='center',va='top')
    plt.plot([first_frost, first_frost],[-5, 7], 'k-', linewidth=1)
    plt.text(first_frost,-5.25,'mean first\nfrost date',fontsize=8,ha='center',va='top')
    plt.text(200,8.5,'growing\ndegree days',fontsize=8,ha='center',va='center')
    plt.text(20,2.5,'chilling\ndays',fontsize=8,ha='center',va='center')
    plt.text(20,-2.5,'freezing\ndays',fontsize=8,ha='center',va='center')
    plt.xlim([doy_min, doy_max])
    plt.xticks(doy_first_dom)
    plt.xlabel('Day of Year')
    plt.ylim([ymin, ymax])
    plt.yticks(rotation='vertical')
    plt.ylabel('Temperature (%sC)' % r'$^\circ$')
    fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
    plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()
    return    

def shifted_cd_plot(avg_pts,min_pts,max_pts,m_plateau_beg,m_plateau_end,m_plateau_cd,titlestr,fname):
    shift_month = 7
    shift_day = 1
    shift_mmdd = shift_month * 100 + shift_day
    date_reset = shift_mmdd; doy_reset = date_to_doy(2013,date_reset)
    doy_min = 1
    doy_max = 365
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy_first_dom = [1]
    for i in range(0,len(days_in_month)-1):
        doy_first_dom.append(doy_first_dom[i] + days_in_month[i])
    doy_first_dom_shift = [doy_reset]
    for i in range(shift_month,12):
        doy_first_dom_shift.append(doy_first_dom[i])
    for i in range(0,shift_month):
        doy_first_dom_shift.append(doy_first_dom[i])
    days_in_month = [31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30]
    doy_first_dom = [1]
    for i in range(0,len(days_in_month)-1):
        doy_first_dom.append(doy_first_dom[i] + days_in_month[i])
    date_reset = shift_mmdd; doy_reset = date_to_doy(2013,date_reset)
    offset = doy_max - doy_reset
    date_vequinox = 320; doy_vequinox = date_to_doy(2013,date_vequinox) + offset
    date_ssolstice = 621; doy_ssolstice = date_to_doy(2013,date_ssolstice) + offset
    date_aequinox = 922; doy_aequinox = date_to_doy(2013,date_aequinox) + offset - doy_max
    date_wsolstice = 1221; doy_wsolstice = date_to_doy(2013,date_wsolstice) + offset - doy_max
    ymin = 0.0
    ymax = 200.0
    plt.plot([doy_vequinox, doy_vequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_vequinox-1,60.0,'vernal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_ssolstice, doy_ssolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_ssolstice-1,60.0,'summer solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_aequinox, doy_aequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_aequinox-1,60.0,'autumnal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_wsolstice, doy_wsolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_wsolstice-1,90.0,'winter solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    avg_pts_shift = []
    for i in range(doy_reset,doy_max):
        avg_pts_shift.append(avg_pts[i])
    for i in range(0,doy_reset-1):
        avg_pts_shift.append(avg_pts[i])
    plt.plot(avg_pts_shift, 'b-', linewidth=2)
    min_pts_shift = []
    for i in range(doy_reset,doy_max):
        min_pts_shift.append(min_pts[i])
    for i in range(0,doy_reset-1):
        min_pts_shift.append(min_pts[i])
    plt.plot(min_pts_shift, 'k--', linewidth=2)
    max_pts_shift = []
    for i in range(doy_reset,doy_max):
        max_pts_shift.append(max_pts[i])
    for i in range(0,doy_reset-1):
        max_pts_shift.append(max_pts[i])
    plt.plot(max_pts_shift, 'k--', linewidth=2)
    plateau_begin = m_plateau_beg + offset
    plateau_end = m_plateau_end + offset - doy_max
    plateau_cd = m_plateau_cd
    plt.plot([plateau_begin, plateau_begin],[plateau_cd-10, plateau_cd+10], 'k-', linewidth=1)
    plt.text(plateau_begin+2,plateau_cd-4,'mean plateau\nbegin date',fontsize=8,ha='left',va='top')
    plt.plot([plateau_end, plateau_end],[0, 20], 'k-', linewidth=1)
    plt.text(plateau_end-2,20,'mean plateau\nend date',fontsize=8,ha='right',va='top')
    plt.plot([offset-90, offset-20],[plateau_cd, plateau_cd], 'k-', linewidth=1)
    plt.text(offset-55,plateau_cd+1,'mean plateau CD',fontsize=8,ha='center',va='bottom')
    plt.xlim([doy_min, doy_max])
    plt.xticks(doy_first_dom,doy_first_dom_shift)
    plt.xlabel('Day of Year')
    plt.ylim([ymin, ymax])
    plt.yticks(rotation='vertical')
    plt.ylabel('Accumulated Chilling Days')
    fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
    plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()  
    return  

def gdd_plot(avg_pts,min_pts,max_pts,m_plateau_beg,m_plateau_end,titlestr,fname):
    doy_min = 1
    doy_max = 365
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy_first_dom = [1]
    for i in range(0,len(days_in_month)-1):
        doy_first_dom.append(doy_first_dom[i] + days_in_month[i])
    doy_first_dom_shift = doy_first_dom
    date_vequinox = 320; doy_vequinox = date_to_doy(2013,date_vequinox)
    date_ssolstice = 621; doy_ssolstice = date_to_doy(2013,date_ssolstice)
    date_aequinox = 922; doy_aequinox = date_to_doy(2013,date_aequinox)
    date_wsolstice = 1221; doy_wsolstice = date_to_doy(2013,date_wsolstice)
    ymin = 0.0
    ymax = 2100.0
    plt.plot([doy_vequinox, doy_vequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_vequinox-1,600.0,'vernal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_ssolstice, doy_ssolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_ssolstice-1,800.0,'summer solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_aequinox, doy_aequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_aequinox-1,600.0,'autumnal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_wsolstice, doy_wsolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_wsolstice-1,600.0,'winter solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot(avg_pts, 'r-', linewidth=2)
    plt.plot(min_pts, 'k--', linewidth=2)
    plt.plot(max_pts, 'k--', linewidth=2)
    plateau_begin = m_plateau_beg 
    plateau_end = m_plateau_end 
    plt.plot([plateau_begin, plateau_begin],[0, 250], 'k-', linewidth=1)
    plt.text(plateau_begin-2,200,'mean plateau\nbegin date',fontsize=8,ha='right',va='bottom')
    plt.plot([plateau_end, plateau_end],[1600, 1850], 'k-', linewidth=1)
    plt.text(plateau_end+2,1700,'mean plateau\nend date',fontsize=8,ha='left',va='top')
    plt.plot([doy_vequinox+10,doy_ssolstice-10],[np.max(avg_pts), np.max(avg_pts)], 'k-', linewidth=1)
    plt.text((doy_vequinox+doy_ssolstice)/2.0,np.max(avg_pts)+1,'mean warm season\nGDD accumulation',fontsize=8,ha='center',va='bottom')
    plt.xlim([doy_min, doy_max])
    plt.xticks(doy_first_dom) #,doy_first_dom_shift)
    plt.xlabel('Day of Year')
    plt.ylim([ymin, ymax])
    plt.yticks(rotation='vertical')
    plt.ylabel('Accumulated Growing Degree Days')
    fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
    plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()  
    return  

def t_ts_bydecade_plot(d1_pts,d2_pts,d3_pts,titlestr,fname):
    doy_min = 1
    doy_max = 365
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy_first_dom = [1]
    for i in range(0,len(days_in_month)-1):
        doy_first_dom.append(doy_first_dom[i] + days_in_month[i])
    date_vequinox = 320; doy_vequinox = date_to_doy(2013,date_vequinox)
    date_ssolstice = 621; doy_ssolstice = date_to_doy(2013,date_ssolstice)
    date_aequinox = 922; doy_aequinox = date_to_doy(2013,date_aequinox)
    date_wsolstice = 1221; doy_wsolstice = date_to_doy(2013,date_wsolstice)
    ymin = -30.0
    ymax = 35.0
    zero_pts = np.zeros(np.shape(d1_pts))
    plt.plot(zero_pts, 'k-', linewidth=1)
    plt.plot([doy_vequinox, doy_vequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_vequinox-1,15.0,'vernal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_ssolstice, doy_ssolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_ssolstice-1,-25.0,'summer solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_aequinox, doy_aequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_aequinox-1,-25.0,'autumnal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_wsolstice, doy_wsolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_wsolstice-1,15.0,'winter solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot(d1_pts, 'b-', linewidth=2, label='1984-1993')
    plt.plot(d2_pts, 'g-', linewidth=2, label='1994-2003')
    plt.plot(d3_pts, 'r-', linewidth=2, label='2004-2013')
    plt.legend(fontsize=8, loc='lower right')
    plt.xlim([doy_min, doy_max])
    plt.xticks(doy_first_dom)
    plt.xlabel('Day of Year')
    plt.ylim([ymin, ymax])
    plt.yticks(rotation='vertical')
    plt.ylabel('Temperature (%sC)' % r'$^\circ$')
    fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
    plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf() 
    return   

def p_ts_bydecade_plot(d1_pts,d2_pts,d3_pts,titlestr,fname):
    doy_min = 1
    doy_max = 365
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy_first_dom = [1]
    for i in range(0,len(days_in_month)-1):
        doy_first_dom.append(doy_first_dom[i] + days_in_month[i])
    date_vequinox = 320; doy_vequinox = date_to_doy(2013,date_vequinox)
    date_ssolstice = 621; doy_ssolstice = date_to_doy(2013,date_ssolstice)
    date_aequinox = 922; doy_aequinox = date_to_doy(2013,date_aequinox)
    date_wsolstice = 1221; doy_wsolstice = date_to_doy(2013,date_wsolstice)
    ymin = 0.0
    ymax = 1.0
    zero_pts = np.zeros(np.shape(d1_pts))
    plt.plot([doy_vequinox, doy_vequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_vequinox-1,0.7,'vernal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_ssolstice, doy_ssolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_ssolstice-1,0.7,'summer solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_aequinox, doy_aequinox],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_aequinox-1,0.7,'autumnal equinox',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot([doy_wsolstice, doy_wsolstice],[ymin, ymax], 'k:', linewidth=1)
    plt.text(doy_wsolstice-1,0.7,'winter solstice',fontsize=10,ha='right',va='bottom',rotation='vertical')
    plt.plot(d1_pts, 'b-', linewidth=2, label='1984-1993')
    plt.plot(d2_pts, 'g-', linewidth=2, label='1994-2003')
    plt.plot(d3_pts, 'r-', linewidth=2, label='2004-2013')
    plt.legend(fontsize=8, loc='upper left')
    plt.xlim([doy_min, doy_max])
    plt.xticks(doy_first_dom)
    plt.xlabel('Day of Year')
    plt.ylim([ymin, ymax])
    plt.yticks(rotation='vertical')
    plt.ylabel('Precipitation (cm)')
    fname_notitle = fname[:-4] + '_notitle' + fname[-4:]
    plt.savefig(fname_notitle, dpi=300, bbox_inches='tight')    
    plt.title(titlestr, fontsize=12)
    message('saving figure as %s' % fname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.clf()  
    return  

# end Plots.py
