"""
Python script 'process_NCEI_03_dag_gen.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: QA/QC of daily meteorological station data in NOAA/NCEI datasets

DEPENDENCIES: none

USAGE: '$ python process_NCEI_03_dag_gen.py'
"""


wxcd_vars = ['preprocess',
             'chill_d', 'chill_dd', 'grow_dd', 'grow_dd_base0',
             'tavg_frz', 'tmax_frz', 'tmin_frz',
             'prcp_03d', 'tavg_03d', 'tmax_03d', 'tmin_03d', 'vpd_03d',
             'prcp_07d', 'tavg_07d', 'tmax_07d', 'tmin_07d', 'vpd_07d',
             'prcp_15d', 'tavg_15d', 'tmax_15d', 'tmin_15d', 'vpd_15d',
             'prcp_30d', 'tavg_30d', 'tmax_30d', 'tmin_30d', 'vpd_30d',
             'prcp_60d', 'tavg_60d', 'tmax_60d', 'tmin_60d', 'vpd_60d',
             'prcp_90d', 'tavg_90d', 'tmax_90d', 'tmin_90d', 'vpd_90d',
             'prcp_90d_nd0', 'prcp_90d_nd10', 'prcp_90d_nd25',
             'prcp_120d', 'prcp_180d', 'prcp_365d']
nvars = len(wxcd_vars)
#
jobs = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        'AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH', 'II', 'JJ', 'KK',
        'LL', 'MM', 'NN', 'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'UU', 'VV',
        'WW', 'XX', 'YY', 'ZZ']
#
years = range(1983, 2014)
nyears = len(years)
#
# create condor_submit_dag file
dag_fname = 'process_NCEI_03_dag.sub'
with open(dag_fname, 'w') as dag_f:
    # create JOB lines
    for j in range(nvars):
        for year in years:
            dag_f.write('JOB %s_%d process_NCEI_03_generic.sub \n' %
                        (jobs[j], year))
    # create VARS lines
    for j, wxcd_var in enumerate(wxcd_vars):
        if '_90d' in wxcd_var:
            memreq = '16GB'
        elif '_120d' in wxcd_var:
            memreq = '16GB'
        elif '_180d' in wxcd_var:
            memreq = '16GB'
        elif '_365d' in wxcd_var:
            memreq = '32GB'
        else:
            memreq = '8GB'
        for year in years:
            dag_f.write('VARS %s_%d var="%s" year="%d" mem="%s" \n' %
                        (jobs[j], year, wxcd_var, year, memreq))
    #
    # create ABORT lines
    for j in range(nvars):
        for year in years:
            dag_f.write('ABORT-DAG-ON %s_%d 1 RETURN 1 \n' % (jobs[j], year))
    #
    # create parent/child relationships within each variable to ensure
    # sequential execution
    for j in range(nvars - 1):
        for i in range(nyears):
            if i == 0:
                dag_f.write('PARENT %s_%d CHILD %s_%d %s_%d \n' %
                            (jobs[j], years[i], jobs[j], years[i + 1],
                             jobs[j + 1], years[i]))
            elif i < (nyears - 1):
                dag_f.write('PARENT %s_%d %s_%d CHILD %s_%d %s_%d \n' %
                            (jobs[j], years[i], jobs[j + 1], years[i - 1],
                             jobs[j], years[i + 1], jobs[j + 1], years[i]))
            else:
                dag_f.write('PARENT %s_%d %s_%d CHILD %s_%d \n' %
                            (jobs[j], years[i], jobs[j + 1], years[i - 1],
                             jobs[j + 1], years[i]))
print 'done!'

# end process_NCEI_03_dag_gen.py
