"""
Python script "setup.py"
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
Treat others as you would be treated. Pay it forward. Valar dohaeris.

USAGE: '$ python setup.py'

PURPOSE: Verifies sample data, scripts, modules, documents, auxiliary files.
         Verifies availability of python dependencies used by various scripts.
         Uncompresses certain large example data files
         Builds directory structure for script output products.

DEPENDENCIES: all software package source dependencies are polled here
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
main_dirs = ['data', 'docs', 'source', 'tools']
#
scripts = []
nscripts = 16
for i in range(0, nscripts):
    scripts.append('process_NCEI_%s.py' % str(i).zfill(2))
scripts.append('process_NCEI_02_preprocess.py')
scripts.append('process_NCEI_02_serial.py')
scripts.append('process_NCEI_03_checkpointed.py')
scripts.append('process_NCEI_03_serial_checkpointed.py')
scripts.append('process_NCEI_04a.py')
scripts.append('process_NCEI_04b.py')
#
modules = ['Date_Convert.py', 'Interpolation.py', 'Plots.py',
           'Read_Header_Files.py', 'Stats.py', 'Teleconnections.py',
           'UTM_Geo_Convert.py']
#
dependencies = ['os', 'sys', 'datetime', 'glob', 'numpy', 'pandas', 'h5py',
                'pp', 'matplotlib.pyplot', 'gdal', 'osgeo.osr',
                'scipy.interpolate', 'scipy.ndimage', 'scipy.stats']
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
tools = ['query_NCEI_grids.py']
#
add_dirs = ['analyses', 'grids', 'images']
#
analyses_dirs = ['annual_maps', 'cluster_maps', 'ecoregion_maps',
                 'figures', 'summary_maps']
message(' ')
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
message('checking for essential python package dependencies for this software')
err = 0
try:
    import os
    message('- python dependency \'os\' is available')
except ImportError:
    message('- essential python dependency \'os\' is not available')
    err += 1
try:
    import sys
    message('- python dependency \'sys\' is available')
except ImportError:
    message('- essential python dependency \'sys\' is not available')
    err += 1
try:
    import datetime
    message('- python dependency \'datetime\' is available')
except ImportError:
    message('- essential python dependency \'datetime\' is not available')
    err += 1
try:
    import glob
    message('- python dependency \'glob\' is available')
except ImportError:
    message('- essential python dependency \'glob\' is not available')
    err += 1
try:
    import numpy
    message('- python dependency \'numpy\' is available')
except ImportError:
    message('- essential python dependency \'numpy\' is not available')
    err += 1
try:
    import pandas
    message('- python dependency \'pandas\' is available')
except ImportError:
    message('- essential python dependency \'pandas\' is not available')
    err += 1
try:
    import h5py
    message('- python dependency \'h5py\' is available')
except ImportError:
    message('- essential python dependency \'h5py\' is not available')
    err += 1
try:
    import pp
    message('- python dependency \'pp\' is available')
except ImportError:
    message('- essential python dependency \'pp\' (ParallelPython) is not \
            available')
    message('  (it\'s essential for now, until I develop serial versions of \
            several routines)')
    err += 1
try:
    import gdal
    message('- python dependency \'gdal\' is available')
except ImportError:
    message('- essential python dependency \'gdal\' is not available')
    err += 1
try:
    import osgeo.osr
    message('- python dependency \'osgeo.osr\' is available')
except ImportError:
    message('- essential python dependency \'osgeo.osr\' is not available')
    err += 1
try:
    import scipy.interpolate
    message('- python dependency \'scipy.interpolate\' is available')
except ImportError:
    message('- essential python dependency \'scipy.interpolate\' is not \
            available')
    err += 1
try:
    import scipy.ndimage
    message('- python dependency \'scipy.ndimage\' is available')
except ImportError:
    message('- essential python dependency \'scipy.ndimage\' is not available')
    err += 1
try:
    import scipy.stats
    message('- python dependency \'scipy.stats\' is available')
except ImportError:
    message('- essential python dependency \'scipy.stats\' is not available')
    err += 1
try:
    import matplotlib.pyplot
    message('- python dependency \'matplotlib.pyplot\' is available')
except ImportError:
    message('- essential python dependency \'matplotlib.pyplot\' is not \
            available')
    err += 1
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
os.system('cp tools/*.py .')
message('archiving original tools scripts to \'tools_orig\' directory')
os.system('mv tools tools_orig')
message(' ')
#
message('all set!')
message(' ')
message('make sure to read the \'README.md\' file before you get started on \
        the scripts')
message(' ')
message('if you need help getting your own dataset of GHCND weather \
        observations, there is')
message('    a how-to document in the \'docs\' directory')
message(' ')
message('please send questions, bug reports, any other requests to \
        matt.e.garcia@gmail.com')
message('    (and include a helpfully descriptive subject line, if you could)')
message('or submit them through the Issues tab at the GitHub repository for \
        this package')
message(' ')
#
sys.exit(0)
