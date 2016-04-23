"""
Python module 'Stats.py'
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

USAGE: insert 'from Stats import *' line near head of script, then call routines as indicated

PURPOSE: Basic statistical calculations that are used the same way more than once

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines

INPUT: Provided by calling script

OUTPUT: Returned to calling script

RUN TIME: negligible
"""

import numpy as np
from scipy.stats import pearsonr, linregress

def regress(ny,yvals):
    xvals = np.arange(0,ny)
    coeffs = np.polyfit(xvals, yvals, 1)
    corr, sig = pearsonr(xvals, yvals)
    return coeffs[0], corr, sig

def getstats(grids,mask,ny):
    grids_reduced = np.zeros((ny))
    for i in range(0,ny):
        grids_reduced[i] = np.nanmean(np.where(mask == 1, grids[i,:,:], np.nan))
    sts = np.zeros((7))
    sts[0] = np.mean(grids_reduced)
    sts[1] = np.std(grids_reduced)
    sts[2] = np.min(grids_reduced)
    sts[3] = np.max(grids_reduced)
    sts[4:7] = regress(ny, grids_reduced)
    return grids_reduced, sts

def indicator_stats(var,nvals):
    var_stats = np.zeros((12))
    # general time series stats
    var_stats[0] = np.min(var)
    var_stats[1] = np.max(var)
    var_stats[2] = np.median(var)    
    var_stats[3] = np.mean(var)
    var_stats[4] = np.std(var)
    # time series autocorrelation
    var_stats[5:7] = pearsonr(var[0:nvals-1], var[1:nvals]) 
    # time series trend
    years = np.arange(0,nvals).astype(int)
    var_stats[7:12] = linregress(years,var)
    return var_stats

# end Stats.py
