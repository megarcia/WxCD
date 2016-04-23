"""
Python module 'Date_Convert.py'
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

USAGE: insert 'from Date_Convert import *' line near head of script, then (for example)
           doy_X = date_to_doy(year, date_X) # date_X is an integer in mmdd format (without leading 0)
           date_X = doy_to_date(year, doy_X) # date_X will be an integer in mmdd format (without leading 0)

PURPOSE: Conversion between DOY and calendar date

DEPENDENCIES: None

INPUT: Provided by calling script

OUTPUT: Returned to calling script

RUN TIME: negligible
"""

def is_leapyear(y):
    leap = 0
    if y % 1000 == 0:
        leap = 1
    elif y % 4 == 0:
        if y % 100 != 0:
            leap = 1
    return leap

def date_to_doy(year, mmdd):
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if is_leapyear(year):
        days_in_month[1] += 1    
    mm = mmdd // 100
    dd = mmdd % 100
    if mm == 1:
        doy = dd
    else:
        doy = sum(days_in_month[:(mm - 1)]) + dd
    return doy

def doy_to_date(year, doy):
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if is_leapyear(year):
        days_in_month[1] += 1    
    accumulated_days = 0  
    for i in range(0,len(days_in_month)):
        next_accumulated_days = accumulated_days + days_in_month[i]
        if doy > next_accumulated_days:
            accumulated_days = next_accumulated_days
        else:
            mm = i + 1
            dd = doy - accumulated_days
            mmdd = (mm * 100) + dd
            break
    return mmdd

# end Date_Convert.py
