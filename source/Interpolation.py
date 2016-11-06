"""
Python module 'Interpolation.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
portions based on python source code found at
http://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
http://stackoverflow.com/questions/12729228/simple-efficient-bilinear-interpolation-of-images-in-numpy-and-python
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: Interpolation from irregular station networks to regular grids

DEPENDENCIES: numpy, scipy

USAGE: insert 'import Interpolation' line near head of script
       see usage examples in 'process_NCEI_02.py'

INPUT: location/value arrays rovided by calling script

OUTPUT: grids returned to calling script
"""


import numpy
from scipy import interpolate


def distance_matrix_idw(x0, y0, x1, y1):
    obs = numpy.vstack((x0, y0)).T
    interp = numpy.vstack((x1, y1)).T
    # Make a distance matrix between pairwise observations
    # Based on example at http://stackoverflow.com/questions/1871536
    dx = numpy.subtract.outer(obs[:, 0], interp[:, 0])
    dy = numpy.subtract.outer(obs[:, 1], interp[:, 1])
    dist = numpy.hypot(dx, dy)
    return dist


def simple_idw(x, y, z, xi, yi, exp):
    dist = distance_matrix_idw(x, y, xi, yi)
    # weights are 1 / distance^exp
    weights = 1.0 / dist**exp
    # normalize weights
    weights /= weights.sum(axis=0)
    # dot-multiply weights for desired locations by original stations
    zi = numpy.dot(weights.T, z)
    return zi


def scipy_griddata(x, y, z, xi, yi, grid_method='cubic'):
    zi = interpolate.griddata((x, y), z, (xi, yi), method=grid_method)
    return zi


def scipy_bspline(x, y, z, xi, yi):
    sbsp = interpolate.bisplrep(x, y, z, s=0.5)
    zi = interpolate.bisplev(xi, yi, sbsp)
    return zi


def scipy_rbf(x, y, z, xi, yi):
    interp = interpolate.Rbf(x, y, z, epsilon=2)
    zi = interp(xi, yi)
    return zi


def bilinear(x0, dx, y0, dy, z, ix, iy):
    ix = numpy.asarray(ix)
    iy = numpy.asarray(iy)
    #
    x1 = x0 + dx
    y1 = y0 + dy
    #
    x0_idx = numpy.clip(x0, 0, len(x0) - 1)
    x1_idx = numpy.clip(x1, 0, len(x0) - 1)
    y0_idx = numpy.clip(y0, 0, len(y0) - 1)
    y1_idx = numpy.clip(y1, 0, len(y0) - 1)
    #
    ima = z[y0_idx][x0_idx]
    imb = z[y1_idx][x0_idx]
    imc = z[y0_idx][x1_idx]
    imd = z[y1_idx][x1_idx]
    #
    wgta = (x1 - ix) * (y1 - iy)
    wgtb = (x1 - ix) * (iy - y0)
    wgtc = (ix - x0) * (y1 - iy)
    wgtd = (ix - x0) * (iy - y0)
    #
    wgtsum = wgta + wgtb + wgtc + wgtd
    wgta /= wgtsum
    wgtb /= wgtsum
    wgtc /= wgtsum
    wgtd /= wgtsum
    #
    zi = wgta * ima + wgtb * imb + wgtc * imc + wgtd * imd
    return zi
