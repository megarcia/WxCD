"""
Python module 'Read_Header_Files.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
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


def get_bil_hdr_info(fname):
    lchdrfile = open(fname, 'r')
    line = lchdrfile.readline()  # BYTEORDER      I
    line = lchdrfile.readline()  # LAYOUT         BIL
    line = lchdrfile.readline()  # NROWS          15657
    line = line.rstrip()
    parts = line.split('          ')
    nr = int(parts[1])
    line = lchdrfile.readline()  # NCOLS          18622
    line = line.rstrip()
    parts = line.split('          ')
    nc = int(parts[1])
    line = lchdrfile.readline()  # NBANDS         1
    line = lchdrfile.readline()  # NBITS          8
    line = lchdrfile.readline()  # BANDROWBYTES   18622
    line = lchdrfile.readline()  # TOTALROWBYTES  18622
    line = lchdrfile.readline()  # PIXELTYPE      UNSIGNEDINT
    line = lchdrfile.readline()  # ULXMAP         340047.378402401
    line = line.rstrip()
    parts = line.split('         ')
    minx = int(round(float(parts[1]), 0))
    line = lchdrfile.readline()  # ULYMAP         5445784.00596
    line = line.rstrip()
    parts = line.split('         ')
    maxy = int(round(float(parts[1]), 0))
    line = lchdrfile.readline()  # XDIM           30
    line = line.rstrip()
    parts = line.split('           ')
    dx = int(parts[1])
    line = lchdrfile.readline()  # YDIM           30
    line = line.rstrip()
    parts = line.split('           ')
    dy = int(parts[1])
    line = lchdrfile.readline()  # NODATA         0
    lchdrfile.close()
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
