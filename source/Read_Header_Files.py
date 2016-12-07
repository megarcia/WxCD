"""
Python module 'Read_Header_Files.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2014-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Obtain grid information from '.hdr' files associated with '.bil'
         datasets (based on the '.hdr' file structure typical to ArcGIS output
         processes)

DEPENDENCIES: None

USAGE: insert 'from Read_Header_Files import *' line near head of script,
       then call individual routine(s) as indicated

INPUT: filename (with path) provided by calling script

OUTPUT: metadata elements returned to calling script
"""


"""
sample contents of Landsat header file

ENVI
description = {LEDAPS HDF File Imported into ENVI}
samples = 8041
lines   = 7201
bands   = 1
header offset = 0
file type = HDF Scientific Data
data type = 2
interleave = bsq
sensor type = Landsat
byte order = 0
map info = {UTM, 1.000, 1.000, 532485.000000, 5365815.000000, 30.000000, 30.000000, 15, North, WGS-84, units=Meters}

end of sample contents of Landsat header file
"""


def get_hdf_hdr_info(fname):
    with open(fname, 'r') as hdrfile:
        line = hdrfile.readline()  # ENVI
        line = hdrfile.readline()  # description
        line = hdrfile.readline()  # samples
        line = line.rstrip()
        parts = line.split('=')
        ncols = int(parts[1])
        line = hdrfile.readline()  # lines
        line = line.rstrip()
        parts = line.split('=')
        nrows = int(parts[1])
        line = hdrfile.readline()  # bands
        line = hdrfile.readline()  # header offset
        line = hdrfile.readline()  # file type
        line = hdrfile.readline()  # data type
        line = hdrfile.readline()  # interleave
        line = hdrfile.readline()  # sensor type
        line = hdrfile.readline()  # byte order
        line = hdrfile.readline()  # map info
        line = line.rstrip()
        parts = line.split('{')
        line = parts[1]
        parts = line.split(', ')
        proj = parts[0]
        NWeasting = float(parts[3])
        NWnorthing = float(parts[4])
        pixelsize = float(parts[5])
        zone = int(parts[7])
        hemi = parts[8]
        datum = parts[9]
        units = parts[10][:len(parts[10]) - 1]
    metadata = [proj, zone, hemi, datum, units]
    SEeasting = NWeasting + (ncols - 1) * pixelsize
    SEnorthing = NWnorthing - (nrows - 1) * pixelsize
    corners = [NWeasting, NWnorthing, ncols, nrows, pixelsize,
               SEeasting, SEnorthing]
    return metadata, corners


"""
sample contents of NLCD BIL-format header file (after reprojection in ArcGIS)

BYTEORDER      I
LAYOUT         BIL
NROWS          15657
NCOLS          18622
NBANDS         1
NBITS          8
BANDROWBYTES   18622
TOTALROWBYTES  18622
PIXELTYPE      UNSIGNEDINT
ULXMAP         340047.378402401
ULYMAP         5445784.00596
XDIM           30
YDIM           30
NODATA         255

end of sample contents of NLCD BIL-format header file
"""


def get_bil_hdr_info(fname):
    with open(fname, 'r') as lchdrfile:
        line = lchdrfile.readline()  # BYTEORDER
        line = lchdrfile.readline()  # LAYOUT
        line = lchdrfile.readline()  # NROWS
        line = line.rstrip()
        parts = line.split('          ')
        nr = int(parts[1])
        line = lchdrfile.readline()  # NCOLS
        line = line.rstrip()
        parts = line.split('          ')
        nc = int(parts[1])
        line = lchdrfile.readline()  # NBANDS
        line = lchdrfile.readline()  # NBITS
        line = lchdrfile.readline()  # BANDROWBYTES
        line = lchdrfile.readline()  # TOTALROWBYTES
        line = lchdrfile.readline()  # PIXELTYPE
        line = lchdrfile.readline()  # ULXMAP
        line = line.rstrip()
        parts = line.split('         ')
        minx = int(round(float(parts[1]), 0))
        line = lchdrfile.readline()  # ULYMAP
        line = line.rstrip()
        parts = line.split('         ')
        maxy = int(round(float(parts[1]), 0))
        line = lchdrfile.readline()  # XDIM
        line = line.rstrip()
        parts = line.split('           ')
        dx = int(parts[1])
        line = lchdrfile.readline()  # YDIM
        line = line.rstrip()
        parts = line.split('           ')
        dy = int(parts[1])
        line = lchdrfile.readline()  # NODATA
    fname_parts = fname.split('_')
    UTMz = 0
    for part in fname_parts:
        if part[:3] == 'UTM':
            UTMz = int(part[3:5])
        if UTMz != 0:
            break
    maxx = minx + nc * dx
    miny = maxy - nr * dy
    return UTMz, nr, nc, miny, maxy, minx, maxx, dy, dx

# end Read_Header_Files.py
