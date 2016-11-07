"""
Python module 'process_NCEI_03_aux.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: helper functions for process_NCEI_03_*.py

DEPENDENCIES: h5py, numpy, pickle

USAGE: insert 'from process_NCEI_03_aux import *' line near head of script
       see usage examples in 'process_NCEI_03_*.py'

INPUT: various, provided by calling script

OUTPUT: various, returned to calling script
"""


import sys
import pickle
import h5py as hdf
import numpy as np


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print char_string
    sys.stdout.flush()
    return


def grid_threshold_count(mm, dd, reset_mm, reset_dd, var_grid, grid_var_prev,
                         stns_all=0, stns=0):
    grid_var_new = np.where(var_grid < 0.0, 1.0, 0.0)
    if mm == reset_mm and dd == reset_dd:
        grid_var = grid_var_new
        if stns_all != 0:
            stns_all = stns
    else:
        grid_var = grid_var_prev + grid_var_new
        if stns_all != 0:
            stns_all = list(set(stns_all) | set(stns))
    if stns_all == 0:
        return grid_var
    else:
        return grid_var, stns_all


def grid_threshold_accumulate(mm, dd, reset_mm, reset_dd, var_grid,
                              grid_var_prev, stns_all, stns):
    grid_var_new = np.where(var_grid > 0.0, var_grid, 0.0)
    if mm == reset_mm and dd == reset_dd:
        grid_var = grid_var_new
        stns_all = stns
    else:
        grid_var = grid_var_prev + grid_var_new
        stns_all = list(set(stns_all) | set(stns))
    return grid_var, stns_all


def cube_sum(nd, var_cube, var_grid, stns_all, stns):
    var_cube[0:nd - 1, :, :] = var_cube[1:nd, :, :]
    var_cube[nd - 1, :, :] = var_grid[:, :]
    grid_var_sum = np.sum(var_cube, axis=0)
    stns_all_dims = np.shape(stns_all)
    if stns_all_dims[0] == nd:
        stns_all[0:nd - 1][:] = stns_all[1:nd][:]
        stns_all[nd - 1] = stns
    else:
        stns_all.append(stns)
    stns_all_list = []
    stns_all_dims = np.shape(stns_all)
    for i in range(stns_all_dims[0]):
        stns_all_list = list(set(stns_all_list) | set(stns_all[i]))
    return grid_var_sum, stns_all_list, var_cube, stns_all


def cube_sum_parts(nd, var_cube1, var_cube2, var_grid, stns_all, stns):
    var_cube1_dims = np.shape(var_cube1)
    nd1 = var_cube1_dims[0]
    var_cube2_dims = np.shape(var_cube2)
    nd2 = var_cube2_dims[0]
    var_cube1[0:nd1 - 1, :, :] = var_cube1[1:nd1, :, :]
    var_cube1[nd1 - 1, :, :] = var_cube2[0, :, :]
    var_cube2[0:nd2 - 1, :, :] = var_cube2[1:nd2, :, :]
    var_cube2[nd2 - 1, :, :] = var_grid[:, :]
    grid_var_sum = np.sum(var_cube1, axis=0) + np.sum(var_cube2, axis=0)
    stns_all_dims = np.shape(stns_all)
    if stns_all_dims[0] == nd:
        stns_all[0:nd - 1][:] = stns_all[1:nd][:]
        stns_all[nd - 1] = stns
    else:
        stns_all.append(stns)
    stns_all_list = []
    stns_all_dims = np.shape(stns_all)
    for i in range(stns_all_dims[0]):
        stns_all_list = list(set(stns_all_list) | set(stns_all[i]))
    return grid_var_sum, stns_all_list, var_cube1, var_cube2, stns_all


def cube_threshold_count(nd, var_cube, var_grid, thresh):
    thresh_var = np.where(var_grid > thresh, 1.0, 0.0)
    var_cube[0:nd - 1, :, :] = var_cube[1:nd, :, :]
    var_cube[nd - 1, :, :] = thresh_var[:, :]
    grid_var_count = np.sum(var_cube, axis=0)
    return grid_var_count, var_cube


def cube_mean_var(nd, var_cube, var_grid, stns_all, stns):
    var_cube[0:nd - 1, :, :] = var_cube[1:nd, :, :]
    var_cube[nd - 1, :, :] = var_grid[:, :]
    grid_var_mean = np.mean(var_cube, axis=0)
    grid_var_var = np.var(var_cube, axis=0)
    stns_all_dims = np.shape(stns_all)
    if stns_all_dims[0] == nd:
        stns_all[0:nd - 1][:] = stns_all[1:nd][:]
        stns_all[nd - 1] = stns
    else:
        stns_all.append(stns)
    stns_all_list = []
    stns_all_dims = np.shape(stns_all)
    for i in range(stns_all_dims[0]):
        stns_all_list = list(set(stns_all_list) | set(stns_all[i]))
    return grid_var_mean, grid_var_var, stns_all_list, var_cube, stns_all


def cube_mean_var_ns(nd, var_cube, var_grid):
    var_cube[0:nd - 1, :, :] = var_cube[1:nd, :, :]
    var_cube[nd - 1, :, :] = var_grid[:, :]
    grid_var_mean = np.mean(var_cube, axis=0)
    grid_var_var = np.var(var_cube, axis=0)
    return grid_var_mean, grid_var_var, var_cube


def write_to_file(h5file, gvar, gdata, svar=0, sdata=0):
    if gvar in h5file.keys():
        del h5file[gvar]
    h5file.create_dataset(gvar, data=gdata, dtype=np.float32,
                          compression='gzip')
    if svar:
        stnsdir = 'stns'
        datapath = '%s/%s' % (stnsdir, svar)
        if svar in h5file[stnsdir].keys():
            del h5file[datapath]
        h5file.create_dataset(datapath, data=sdata)
        message('- %s %s with %d stations' %
                (gvar, str(gdata.shape), len(sdata)))
    return


def write_to_file_2g(h5file, gvar1, gdata1, gvar2, gdata2, svar=0, sdata=0):
    if gvar1 in h5file.keys():
        del h5file[gvar1]
    if gvar2 in h5file.keys():
        del h5file[gvar2]
    h5file.create_dataset(gvar1, data=gdata1, dtype=np.float32,
                          compression='gzip')
    h5file.create_dataset(gvar2, data=gdata2, dtype=np.float32,
                          compression='gzip')
    if svar:
        stnsdir = 'stns'
        datapath = '%s/%s' % (stnsdir, svar)
        if svar in h5file[stnsdir].keys():
            del h5file[datapath]
        h5file.create_dataset(datapath, data=sdata)
        message('- %s %s and %s %s with %d stations' %
                (gvar1, str(gdata1.shape), gvar2,
                 str(gdata2.shape), len(sdata)))
    return


def write_stn_lists(path, year, var, contents):
    fname = '%s/%d_%s.pickle' % (path, year, var)
    f = open(fname, 'wb')
    pickle.dump(contents, f)
    f.close()
    return


def get_stn_lists(path, year, var):
    fname = '%s/%d_%s.pickle' % (path, year, var)
    f = open(fname, 'rb')
    contents = pickle.load(f)
    f.close()
    return contents

# end process_NCEI_03_aux.py
