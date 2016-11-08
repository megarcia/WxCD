"""
Python script "setup.py"
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Verifies sample data, scripts, modules, documents, auxiliary files.
         Verifies availability of python dependencies used by various scripts.
         Uncompresses certain large example data files
         Builds directory structure for script output products.

DEPENDENCIES: all software package source dependencies are polled here

USAGE: '$ python setup.py'
"""


import os
import sys
import glob


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


txt_files = ['ACKNOWLEDGEMENTS.txt', 'CITATION.txt', 'DISCLAIMER.txt',
             'LICENSE_GnuGPLv3.txt']
md_files = ['README.md']
main_dirs = ['data', 'docs', 'htcondor', 'source', 'tools']
#
scripts = ['process_NCEI_00.py', 'process_NCEI_01.py',
           'process_NCEI_02a.py', 'process_NCEI_02b.py',
           'process_NCEI_03_chill_d.py', 'process_NCEI_03_chill_dd.py',
           'process_NCEI_03_grow_dd.py', 'process_NCEI_03_grow_dd_base0.py',
           'process_NCEI_03_prcp_03d.py', 'process_NCEI_03_prcp_07d.py',
           'process_NCEI_03_prcp_120d.py', 'process_NCEI_03_prcp_15d.py',
           'process_NCEI_03_prcp_180d.py', 'process_NCEI_03_prcp_30d.py',
           'process_NCEI_03_prcp_365d.py', 'process_NCEI_03_prcp_60d.py',
           'process_NCEI_03_prcp_90d.py', 'process_NCEI_03_prcp_90d_nd0.py',
           'process_NCEI_03_prcp_90d_nd10.py',
           'process_NCEI_03_prcp_90d_nd25.py',
           'process_NCEI_03_preprocess.py', 'process_NCEI_03_tavg_03d.py',
           'process_NCEI_03_tavg_07d.py', 'process_NCEI_03_tavg_15d.py',
           'process_NCEI_03_tavg_30d.py', 'process_NCEI_03_tavg_60d.py',
           'process_NCEI_03_tavg_90d.py', 'process_NCEI_03_tavg_frz.py',
           'process_NCEI_03_tmax_03d.py', 'process_NCEI_03_tmax_07d.py',
           'process_NCEI_03_tmax_15d.py', 'process_NCEI_03_tmax_30d.py',
           'process_NCEI_03_tmax_60d.py', 'process_NCEI_03_tmax_90d.py',
           'process_NCEI_03_tmax_frz.py', 'process_NCEI_03_tmin_03d.py',
           'process_NCEI_03_tmin_07d.py', 'process_NCEI_03_tmin_15d.py',
           'process_NCEI_03_tmin_30d.py', 'process_NCEI_03_tmin_60d.py',
           'process_NCEI_03_tmin_90d.py', 'process_NCEI_03_tmin_frz.py',
           'process_NCEI_03_vpd_03d.py', 'process_NCEI_03_vpd_07d.py',
           'process_NCEI_03_vpd_15d.py', 'process_NCEI_03_vpd_30d.py',
           'process_NCEI_03_vpd_60d.py', 'process_NCEI_03_vpd_90d.py',
           'process_NCEI_04a.py', 'process_NCEI_04b.py', 'process_NCEI_05.py',
           'process_NCEI_06.py', 'process_NCEI_07.py', 'process_NCEI_08.py',
           'process_NCEI_09.py', 'process_NCEI_10.py', 'process_NCEI_11.py',
           'process_NCEI_12.py', 'process_NCEI_13.py', 'process_NCEI_14.py',
           'process_NCEI_15.py']
#
modules = ['Date_Convert.py', 'Interpolation.py', 'Plots.py',
           'process_NCEI_03_aux.py', 'Read_Header_Files.py', 'Stats.py',
           'Teleconnections.py', 'UTM_Geo_Convert.py']
#
htcondor = ['process_NCEI_00.sh', 'process_NCEI_00.sub',
            'process_NCEI_01.sh', 'process_NCEI_01.sub',
            'process_NCEI_02a.sh', 'process_NCEI_02a.sub',
            'process_NCEI_02b.sh', 'process_NCEI_02b.sub',
            'process_NCEI_02b_dag.sub', 'process_NCEI_03_chill_d.sh',
            'process_NCEI_03_chill_dd.sh', 'process_NCEI_03_dag_gen.py',
            'process_NCEI_03_generic.sub', 'process_NCEI_03_grow_dd.sh',
            'process_NCEI_03_grow_dd_base0.sh', 'process_NCEI_03_prcp_03d.sh',
            'process_NCEI_03_prcp_07d.sh', 'process_NCEI_03_prcp_120d.sh',
            'process_NCEI_03_prcp_15d.sh', 'process_NCEI_03_prcp_180d.sh',
            'process_NCEI_03_prcp_30d.sh', 'process_NCEI_03_prcp_365d.sh',
            'process_NCEI_03_prcp_60d.sh', 'process_NCEI_03_prcp_90d.sh',
            'process_NCEI_03_prcp_90d_nd0.sh',
            'process_NCEI_03_prcp_90d_nd10.sh',
            'process_NCEI_03_prcp_90d_nd25.sh',
            'process_NCEI_03_preprocess.sh', 'process_NCEI_03_tavg_03d.sh',
            'process_NCEI_03_tavg_07d.sh', 'process_NCEI_03_tavg_15d.sh',
            'process_NCEI_03_tavg_30d.sh', 'process_NCEI_03_tavg_60d.sh',
            'process_NCEI_03_tavg_90d.sh', 'process_NCEI_03_tavg_frz.sh',
            'process_NCEI_03_tmax_03d.sh', 'process_NCEI_03_tmax_07d.sh',
            'process_NCEI_03_tmax_15d.sh', 'process_NCEI_03_tmax_30d.sh',
            'process_NCEI_03_tmax_60d.sh', 'process_NCEI_03_tmax_90d.sh',
            'process_NCEI_03_tmax_frz.sh', 'process_NCEI_03_tmin_03d.sh',
            'process_NCEI_03_tmin_07d.sh', 'process_NCEI_03_tmin_15d.sh',
            'process_NCEI_03_tmin_30d.sh', 'process_NCEI_03_tmin_60d.sh',
            'process_NCEI_03_tmin_90d.sh', 'process_NCEI_03_tmin_frz.sh',
            'process_NCEI_03_vpd_03d.sh', 'process_NCEI_03_vpd_07d.sh',
            'process_NCEI_03_vpd_15d.sh', 'process_NCEI_03_vpd_30d.sh',
            'process_NCEI_03_vpd_60d.sh', 'process_NCEI_03_vpd_90d.sh',
            'process_NCEI_04a.sh', 'process_NCEI_04a.sub',
            'process_NCEI_04b.sh', 'process_NCEI_04b.sub',
            'process_NCEI_05.sh', 'process_NCEI_05.sub',
            'process_NCEI_06.sh', 'process_NCEI_06.sub',
            'process_NCEI_07.sh', 'process_NCEI_07.sub',
            'process_NCEI_08.sh', 'process_NCEI_08.sub',
            'process_NCEI_09.sh', 'process_NCEI_09.sub']
#
dependencies = ['os', 'sys', 'datetime', 'glob', 'numpy', 'pandas', 'h5py',
                'matplotlib', 'matplotlib.pyplot', 'gdal', 'osgeo.osr',
                'scipy.interpolate', 'scipy.ndimage', 'scipy.stats',
                'mpl_toolkits', 'mpl_toolkits.basemap', 'pickle']
#
gz_data_files = ['EPA_L4_Ecoregions_WLS_UTM15N.bil.gz',
                 'NCEI_WLS_19830101-20151031.csv.gz',
                 'NLCD_2011_WLS_UTM15N.bil.gz']
#
data_files = ['EPA_L4_Ecoregions_WLS_polygonIDs.txt',
              'EPA_L4_Ecoregions_WLS_UTM15N.bil',
              'EPA_L4_Ecoregions_WLS_UTM15N.hdr',
              'NCEI_WLS_19830101-20151031.csv',
              'NCEP_CPC_AO_indices.csv',
              'NCEP_CPC_ENSO_indices.csv',
              'NCEP_CPC_NAO_indices.csv',
              'NCEP_CPC_PNA_indices.csv',
              'NLCD_2011_WLS_UTM15N.bil',
              'NLCD_2011_WLS_UTM15N.hdr',
              'NOAA_ESRL_AMO_indices.csv',
              'NOAA_ESRL_PDO_indices.csv',
              'NSIDC_MIFL_Superior_Ice.csv',
              'Query_locations_dates_sample.csv']
#
doc_files = ['How_to_get_NCEI_GHCND_data.txt',
             'NCEI_GHCND_documentation.pdf']
#
tools = ['query_NCEI_grids.py', 'orientation_maps.py']
#
add_dirs = ['analyses', 'grids', 'images']
#
analyses_dirs = ['annual_maps', 'cluster_maps', 'ecoregion_maps',
                 'figures', 'summary_maps']
#
os.system('rm .DS_Store')
os.system('rm */.DS_Store')
os.system('rm ._*')
os.system('rm */._*')
#
message('checking for auxiliary files that should accompany this software')
txts_present = glob.glob('*.txt')
mds_present = glob.glob('*.md')
absent = 0
for txt in txt_files:
    if txt in txts_present:
        message('- found auxiliary file \'%s\' as expected' % txt)
    else:
        message('- auxiliary file \'%s\' is absent' % txt)
        absent += 1
for md in md_files:
    if md in mds_present:
        message('- found auxiliary file \'%s\' as expected' % md)
    else:
        message('- auxiliary file \'%s\' is absent' % md)
        absent += 1
if absent > 0:
    message('- you don\'t need them to run things, but you do need them to \
            understand things')
    message('- you should probably download this package again from scratch')
    message('- exiting setup procedure')
    sys.exit(1)
message(' ')
#
message('checking for top-level directories that should already exist')
dirs_present = [d.replace('/', '') for d in glob.glob('*/')]
absent = 0
for dirname in main_dirs:
    if dirname in dirs_present:
        message('- found main directory \'%s\' as expected' % dirname)
    else:
        message('- main directory \'%s\' is absent' % dirname)
        absent += 1
if absent > 0:
    message('- you should download this package again from scratch')
    message('- exiting setup procedure')
    sys.exit(1)
message(' ')
#
message('checking for main scripts and modules that comprise this software')
src_present = glob.glob('source/*')
absent = 0
for srcfile in scripts:
    srcfile = 'source/%s' % srcfile
    if srcfile in src_present:
        message('- found script \'%s\' as expected' % srcfile)
    else:
        message('- script \'%s\' is absent' % srcfile)
        absent += 1
for srcfile in modules:
    srcfile = 'source/%s' % srcfile
    if srcfile in src_present:
        message('- found module \'%s\' as expected' % srcfile)
    else:
        message('- module \'%s\' is absent' % srcfile)
        absent += 1
if absent > 0:
    message('- you should download this package again from scratch')
    message('- exiting setup procedure')
    sys.exit(1)
message(' ')
#
message('checking for script-based tools that accompany this software')
src_present = glob.glob('tools/*')
absent = 0
for srcfile in tools:
    srcfile = 'tools/%s' % srcfile
    if srcfile in src_present:
        message('- found script \'%s\' as expected' % srcfile)
    else:
        message('- script \'%s\' is absent' % srcfile)
        absent += 1
if absent > 0:
    message('- if you need these tools, you should download this package \
            again from scratch')
message(' ')
#
message('checking for HTCondor example files that accompany this software')
src_present = glob.glob('htcondor/*')
absent = 0
for srcfile in htcondor:
    srcfile = 'htcondor/%s' % srcfile
    if srcfile in src_present:
        message('- found htcondor file \'%s\' as expected' % srcfile)
    else:
        message('- htcondor file \'%s\' is absent' % srcfile)
        absent += 1
if absent > 0:
    message('- if you need these files, you should download this package \
            again from scratch')
message(' ')
#
message('checking for essential python package dependencies for this software')
err = 0
#
try:
    import os
    message('- python dependency \'os\' is available')
except ImportError:
    message('- essential python dependency \'os\' is not available')
    err += 1
#
try:
    import sys
    message('- python dependency \'sys\' is available')
except ImportError:
    message('- essential python dependency \'sys\' is not available')
    err += 1
#
try:
    import datetime
    message('- python dependency \'datetime\' is available')
except ImportError:
    message('- essential python dependency \'datetime\' is not available')
    err += 1
#
try:
    import glob
    message('- python dependency \'glob\' is available')
except ImportError:
    message('- essential python dependency \'glob\' is not available')
    err += 1
#
try:
    import pickle
    message('- python dependency \'pickle\' is available')
except ImportError:
    message('- essential python dependency \'pickle\' is not available')
    err += 1
#
try:
    import numpy
    message('- python dependency \'numpy\' is available')
except ImportError:
    message('- essential python dependency \'numpy\' is not available')
    err += 1
#
try:
    import pandas
    message('- python dependency \'pandas\' is available')
except ImportError:
    message('- essential python dependency \'pandas\' is not available')
    err += 1
#
try:
    import h5py
    message('- python dependency \'h5py\' is available')
except ImportError:
    message('- essential python dependency \'h5py\' is not available')
    err += 1
#
try:
    import gdal
    message('- python dependency \'gdal\' is available')
except ImportError:
    message('- essential python dependency \'gdal\' is not available')
    err += 1
#
try:
    import osgeo.osr
    message('- python dependency \'osgeo.osr\' is available')
except ImportError:
    message('- essential python dependency \'osgeo.osr\' is not available')
    err += 1
#
try:
    import scipy.interpolate
    message('- python dependency \'scipy.interpolate\' is available')
except ImportError:
    message('- essential python dependency \'scipy.interpolate\' is not \
            available')
    err += 1
#
try:
    import scipy.ndimage
    message('- python dependency \'scipy.ndimage\' is available')
except ImportError:
    message('- essential python dependency \'scipy.ndimage\' is not available')
    err += 1
#
try:
    import scipy.stats
    message('- python dependency \'scipy.stats\' is available')
except ImportError:
    message('- essential python dependency \'scipy.stats\' is not available')
    err += 1
#
try:
    import matplotlib
    message('- python dependency \'matplotlib\' is available')
except ImportError:
    message('- essential python dependency \'matplotlib\' is not available')
    err += 1
#
try:
    import matplotlib.pyplot
    message('- python dependency \'matplotlib.pyplot\' is available')
except ImportError:
    message('- essential python dependency \'matplotlib.pyplot\' is not \
            available')
    err += 1
#
try:
    import mpl_toolkits
    message('- python dependency \'mpl_toolkits\' is available')
except ImportError:
    message('- essential python dependency \'mpl_toolkits\' is not available')
    err += 1
#
try:
    import mpl_toolkits.basemap
    message('- python dependency \'mpl_toolkits.basemap\' is available')
except ImportError:
    message('- essential python dependency \'mpl_toolkits.basemap\' is not \
            available')
    err += 1
#
if err > 0:
    message('- you need to install one or more additional python packages for \
            this software to work')
    message('- all of these packages are available via Anaconda (\'conda\') \
            and/or PyPI (\'pip\') repositories')
    message('- exiting setup procedure')
    sys.exit(1)
message(' ')
#
message('checking for example data files that should accompany this software')
gz_data_present = glob.glob('data/*.gz')
absent = 0
for gz_dfile in gz_data_files:
    gz_dfile_path = 'data/%s' % gz_dfile
    if gz_dfile_path in gz_data_present:
        message('- found compressed data file \'%s\' as expected' % gz_dfile)
        message('-- uncompressing \'%s\'' % gz_dfile)
        os.system('cd data')
        os.system('gunzip %s' % gz_dfile)
        os.system('cd ..')
    else:
        message('- compressed example data file \'%s\' is absent' % gz_dfile)
        absent += 1
if absent > 0:
    message('- you don\'t need these if you have your own data in the right \
            formats')
    message('- if you need the examples, you can find them at on GitHub at')
    message('      https://github.com/megarcia/WxCD')
#
data_present = glob.glob('data/*')
absent = 0
for dfile in data_files:
    dfile_path = 'data/%s' % dfile
    if dfile_path in data_present:
        message('- found data file \'%s\' as expected' % dfile)
    else:
        message('- example data file \'%s\' is absent' % dfile)
        absent += 1
if absent > 0:
    message('- you don\'t need these if you have your own data in the right \
            formats')
    message('- if you need the examples, you can find them at on GitHub at')
    message('      https://github.com/megarcia/WxCD')
message(' ')
#
message('checking for data documentation files that should accompany this \
        software')
docs_present = glob.glob('docs/*')
absent = 0
for dfile in doc_files:
    dfile = 'docs/%s' % dfile
    if dfile in docs_present:
        message('- found documentation file \'%s\' as expected' % dfile)
    else:
        message('- data documentation file \'%s\' is absent' % dfile)
        absent += 1
if absent > 0:
    message('- you don\'t need these if you have your own documentation')
    message('- if you need the examples, you can find them at on GitHub at')
    message('      https://github.com/megarcia/GT16_JGRA')
message(' ')
#
message('creating top-level and sub-directories that will be used for process \
        output')
for dirname in add_dirs:
    os.system('mkdir %s' % dirname)
    message('- made top-level directory \'%s\' ' % dirname)
for dirname in analyses_dirs:
    os.system('mkdir analyses/%s' % dirname)
    message('- made sub-directory \'analyses/%s\' ' % dirname)
message(' ')
#
message('copying source scripts and modules to top-level directory')
os.system('cp source/*.py .')
message('archiving original scripts and modules to \'source_orig\' directory')
os.system('mv source source_orig')
#
message('copying tools to top-level directory')
os.system('cp tools/*.py .')
message('archiving original tools scripts to \'tools_orig\' directory')
os.system('mv tools tools_orig')
message(' ')
#
message('all set!')
message(' ')
#
message('if you plan to use the HTCondor example files, you\'ll need to \
        move or copy them to')
message('    your top-level directory')
message(' ')
#
message('make sure to read the \'README.md\' file before you get started on \
        the scripts')
message(' ')
#
message('if you need help getting your own dataset of GHCND weather \
        observations, there is')
message('    a how-to document in the \'docs\' directory')
message(' ')
#
message('please send questions, bug reports, any other requests to \
        matt.e.garcia@gmail.com')
message('    (and include a helpfully descriptive subject line, if you could)')
message('or submit them through the Issues tab at the GitHub repository for \
        this package')
message(' ')
#
sys.exit(0)
