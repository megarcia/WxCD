"""
Python module 'Read_Header_Files.py'
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

USAGE: insert 'from Read_Header_Files import *' line near head of script, then call 
       individual routine(s) as indicated

PURPOSE: Getting grid information from '.hdr' files associated with '.bil' datasets 
         (based on the '.hdr' file structure typical to ArcGIS output processes)

DEPENDENCIES: None

INPUT: Provided by calling script

OUTPUT: Returned to calling script

RUN TIME: negligible
"""

def get_bil_hdr_info(fname):
    lchdrfile = open(fname,'r')
    line = lchdrfile.readline() # BYTEORDER      I
    line = lchdrfile.readline() # LAYOUT         BIL
    line = lchdrfile.readline() # NROWS          15657
    line = line.rstrip()
    parts = line.split('          ')
    nr = int(parts[1])
    line = lchdrfile.readline() # NCOLS          18622
    line = line.rstrip()
    parts = line.split('          ')
    nc = int(parts[1])
    line = lchdrfile.readline() # NBANDS         1
    line = lchdrfile.readline() # NBITS          8
    line = lchdrfile.readline() # BANDROWBYTES   18622
    line = lchdrfile.readline() # TOTALROWBYTES  18622
    line = lchdrfile.readline() # PIXELTYPE      UNSIGNEDINT
    line = lchdrfile.readline() # ULXMAP         340047.378402401
    line = line.rstrip()
    parts = line.split('         ')
    minx = int(round(float(parts[1]),0))
    line = lchdrfile.readline() # ULYMAP         5445784.00596
    line = line.rstrip()
    parts = line.split('         ')
    maxy = int(round(float(parts[1]),0))
    line = lchdrfile.readline() # XDIM           30
    line = line.rstrip()
    parts = line.split('           ')
    dx = int(parts[1])
    line = lchdrfile.readline() # YDIM           30
    line = line.rstrip()
    parts = line.split('           ')
    dy = int(parts[1])
    line = lchdrfile.readline() # NODATA         0
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
