"""
Python script 'process_NCEI_08.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Generate land mask for maps, based on map of EPA Level-IV ecoregions

DEPENDENCIES: h5py, numpy, scipy.ndimage
              'Read_Header_Files' module has no external requirements
              'Plots' module has its own requirements

USAGE: '$ python process_NCEI_08.py EPA_L4_Ecoregions_WLS_UTM15N ./data'

INPUT: '.bil' raster map of ecoregions with corresponding '.hdr' header
       At least one '.h5' output file from process_NCEI_03.py
       (with the naming convention 'grids/[YYYYMMDD]_NCEI_grids_2.h5')

OUTPUT: 'clipped_ecoregions.h5' in data folder
"""


import sys
import datetime
import glob
import h5py as hdf
import numpy as np
import scipy.ndimage.interpolation
from Read_Header_Files import get_bil_hdr_info
from Plots import masked_map_plot_geo


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


message(' ')
message('process_NCEI_08.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 3:
    message('input warning: no input directory indicated, using ./data')
    path = './data'
else:
    path = sys.argv[2]
#
if len(sys.argv) < 2:
    message('input error: need prefix for ecoregion map .hdr file')
    sys.exit(1)
else:
    EPAfile = sys.argv[1]
EPAhdrfile = '%s/%s.hdr' % (path, EPAfile)
EPAbilfile = '%s/%s.bil' % (path, EPAfile)
#
# get working area size/shape/location from a weather derivatives files
wxlist = glob.glob('%s/../grids/*_NCEI_grids_2.h5' % path)
message('found %d weather derivative grid files' % len(wxlist))
message(' ')
UTM_bounds = []
message('extracting grid information from %s' % wxlist[0])
with hdf.File(wxlist[0], 'r') as h5infile:
    UTM_zone = np.copy(h5infile['grid/UTMzone'])
    NWeasting = np.copy(h5infile['grid/min_x'])
    UTM_bounds.append(float(NWeasting))
    NWnorthing = np.copy(h5infile['grid/max_y'])
    UTM_bounds.append(float(NWnorthing))
    SEeasting = np.copy(h5infile['grid/max_x'])
    UTM_bounds.append(float(SEeasting))
    SEnorthing = np.copy(h5infile['grid/min_y'])
    UTM_bounds.append(float(SEnorthing))
    dx = np.copy(h5infile['grid/dx'])
    dy = np.copy(h5infile['grid/dy'])
    ncols = np.copy(h5infile['grid/ncols'])
    nrows = np.copy(h5infile['grid/nrows'])
message('- UTM zone %d' % UTM_zone)
message('- Y from %d to %d (%d rows)' % (SEnorthing, NWnorthing, nrows))
message('- X from %d to %d (%d cols)' % (NWeasting, SEeasting, ncols))
message('- dx = dy = %d' % dx)
message(' ')
#
eco_corners = []
corners_tags = ['NW corner easting', 'NW corner northing', 'ncols', 'nrows',
                'pixel size', 'SE corner easting', 'SE corner northing']
eco_clipbounds = []
clipbounds_tags = ['W boundary', 'N boundary', 'E boundary', 'S boundary',
                   'W column', 'N row', 'E column', 'S row',
                   'ncols after clip', 'nrows after clip']
#
# get ecoregion map size/shape/location
message('extracting header information from %s' % EPAhdrfile)
UTMzone, eco_nrows, eco_ncols, eco_SEnorthing, eco_NWnorthing, eco_NWeasting, \
    eco_SEeasting, eco_dy, eco_dx = get_bil_hdr_info(EPAhdrfile)
eco_pixelsize = eco_dx
#
eco_corners.append(eco_NWeasting)   # j = 0
eco_corners.append(eco_NWnorthing)  # j = 1
eco_corners.append(eco_ncols)       # j = 2
eco_corners.append(eco_nrows)       # j = 3
eco_corners.append(eco_pixelsize)   # j = 4
eco_corners.append(eco_SEeasting)   # j = 5
eco_corners.append(eco_SEnorthing)  # j = 6
message('- ecoregion map corners: %s' % str(eco_corners))
#
eco_clipbounds.append(UTM_bounds[0])  # j = 0, W boundary
eco_clipbounds.append(UTM_bounds[1])  # j = 1, N boundary
eco_clipbounds.append(UTM_bounds[2])  # j = 2, E boundary
eco_clipbounds.append(UTM_bounds[3])  # j = 3, S boundary
enlarge = 0
#
Wcol = int(round(eco_clipbounds[0] - eco_corners[0]) / eco_pixelsize)
if Wcol < 0:
    message('-- ecoregion map boundary error: Wcol = %d < 0' % Wcol)
    eco_W = 0
    enlarge = 1
else:
    eco_W = Wcol
    message('eco_W column = %d' % eco_W)
eco_clipbounds.append(eco_W)  # j = 4
#
Nrow = int(round(eco_corners[1] - eco_clipbounds[1]) / eco_pixelsize)
if Nrow < 0:
    message('-- ecoregion map boundary error: Nrow = %d < 0' % Nrow)
    eco_N = 0
    enlarge = 1
else:
    eco_N = Nrow
    message('eco_N row = %d' % eco_N)
eco_clipbounds.append(eco_N)  # j = 5
#
Ecol = int(Wcol + (ncols * dx / eco_pixelsize))
if Ecol > eco_ncols:
    message('-- ecoregion map boundary error: Ecol = %d > eco_ncols = %d' %
            (Ecol, eco_ncols))
    eco_E = eco_W + eco_ncols
    enlarge = 1
else:
    eco_E = Ecol
    message('eco_E column = %d' % eco_E)
eco_clipbounds.append(Ecol)  # j = 6
#
Srow = int(Nrow + (nrows * dy / eco_pixelsize))
if Srow > eco_nrows:
    message('-- ecoregion map boundary error: Srow = %d > eco_nrows = %d' %
            (Srow, eco_nrows))
    eco_S = eco_N + eco_nrows
    enlarge = 1
else:
    eco_S = Srow
    message('eco_S row = %d' % eco_S)
eco_clipbounds.append(Srow)  # j = 7
#
eco_ncols_clip = Ecol - Wcol
message('eco_ncols_clip = %d' % eco_ncols_clip)
eco_clipbounds.append(eco_ncols_clip)  # j = 8
#
eco_nrows_clip = Srow - Nrow
message('eco_nrows_clip = %d' % eco_nrows_clip)
eco_clipbounds.append(eco_nrows_clip)  # j = 9
message('- ecoregion clip info: %s' % str(eco_clipbounds))
message(' ')
#
# clip ecoregion map to match climate derivatives domain
message('extracting ecoregion map from %s' % EPAbilfile)
with open(EPAbilfile, 'r') as ecobilfile:
    eco_raw = np.fromfile(file=ecobilfile, dtype=np.int16).reshape(eco_nrows,
                                                                   eco_ncols)
eco_clip = np.zeros((eco_nrows_clip, eco_ncols_clip))
if enlarge:
    message('- clipping/expanding ecoregion map to calculated boundaries')
    eco_clip[Nrow:Srow, Wcol:Ecol] = eco_raw[eco_N:eco_S, eco_W:eco_E]
else:
    message('- clipping ecoregion map to calculated boundaries')
    eco_clip[:, :] = eco_raw[eco_N:eco_S, eco_W:eco_E]
eco_clip = np.flipud(eco_clip)
eco_clip_rows, eco_clip_cols = np.shape(eco_clip)
message('- intermediate ecoregion map dimensions: %d rows, %d cols' %
        (eco_clip_rows, eco_clip_cols))
message('- generating full-resolution land mask')
landmask = np.where(eco_clip != 0, 1, 0)
message(' ')
#
outfile = '%s/clipped_ecoregions.h5' % path
message('saving clipped ecoregion map to %s' % outfile)
with hdf.File(outfile, 'w') as h5outfile:
    h5outfile.create_dataset('meta/filename', data=outfile)
    h5outfile.create_dataset('meta/created',
                             data=datetime.datetime.now().isoformat())
    h5outfile.create_dataset('meta/by',
                             data='M. Garcia, UWisconsin-Madison FWE')
    h5outfile.create_dataset('meta/last_updated',
                             data=datetime.datetime.now().isoformat())
    h5outfile.create_dataset('meta/sourcebil', data=EPAbilfile)
    h5outfile.create_dataset('meta/sourcebilhdr', data=EPAhdrfile)
    message('- 6 metadata items saved')
    h5outfile.create_dataset('grid/UTM_zone', data=UTM_zone)
    h5outfile.create_dataset('grid/corners_orig', data=eco_corners)
    h5outfile.create_dataset('grid/pixelsize_orig', data=eco_pixelsize)
    h5outfile.create_dataset('grid/clip_bounds', data=eco_clipbounds)
    message('- 4 grid items saved')
    h5outfile.create_dataset('eco_clip', data=eco_clip, dtype=np.float32,
                             compression='gzip')
    h5outfile.create_dataset('landmask', data=landmask, dtype=np.int8,
                             compression='gzip')
message(' ')
#
message('zoom operation to expand/reduce clipped ecoregion map resolution')
zoom_factor = float(eco_pixelsize) / float(dx)
message('- zoom factor = %.2f' % zoom_factor)
eco_clip_reduced = scipy.ndimage.interpolation.zoom(eco_clip, zoom_factor,
                                                    order=0)
eco_clip_reduced_rows, eco_clip_reduced_cols = np.shape(eco_clip_reduced)
message('- final ecoregion map dimensions: %d rows, %d cols' %
        (eco_clip_reduced_rows, eco_clip_reduced_cols))
message('- generating reduced-resolution land mask')
landmask_reduced = np.where(eco_clip_reduced != 0, 1, 0)
message(' ')
#
message('saving expanded/reduced ecoregion map to %s' % outfile)
with hdf.File(outfile, 'r+') as h5outfile:
    h5outfile.create_dataset('grid/pixelsize_reduced', data=dx)
    h5outfile.create_dataset('grid/zoom_factor', data=zoom_factor)
    message('- 2 grid metadata items saved')
    h5outfile.create_dataset('eco_clip_reduced', data=eco_clip_reduced,
                             dtype=np.float32, compression='gzip')
    h5outfile.create_dataset('landmask_reduced', data=landmask_reduced,
                             dtype=np.int8, compression='gzip')
    message('- 2 grids saved')
message(' ')
#
UTM_bounds = eco_clipbounds[:4]
titlestr = 'Ecoregions'
fname = '%s/ecoregions.png' % path
masked_map_plot_geo(eco_clip, landmask, UTM_zone, UTM_bounds, 'rainbow',
                    'none', titlestr, fname)
fname = '%s/ecoregions_reduced.png' % path
masked_map_plot_geo(eco_clip_reduced, landmask_reduced, UTM_zone, UTM_bounds,
                    'rainbow', 'none', titlestr, fname)
titlestr = 'Land Mask'
fname = '%s/landmask.png' % path
masked_map_plot_geo(landmask, landmask, UTM_zone, UTM_bounds, 'rainbow',
                    'none', titlestr, fname)
fname = '%s/landmask_reduced.png' % path
masked_map_plot_geo(landmask_reduced, landmask_reduced, UTM_zone, UTM_bounds,
                    'rainbow', 'none', titlestr, fname)
#
message('process_NCEI_08.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_08.py
