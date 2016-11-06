"""
Python module 'Date_Convert.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Conversion between DOY and calendar date

DEPENDENCIES: None

USAGE: insert 'from Date_Convert import *' near head of script, then
       (for example)
        # date_X is an integer in mmdd format (without leading 0)
        doy_X = date_to_doy(year, date_X)
        # date_X will be an integer in mmdd format (without leading 0)
        date_X = doy_to_date(year, doy_X)

INPUT: date/doy provided by calling script

OUTPUT: date/doy returned to calling script
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
    for i in range(0, len(days_in_month)):
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
