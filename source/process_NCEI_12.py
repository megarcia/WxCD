"""
Python script 'process_NCEI_12.py'
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
               Garcia, M., and P.A. Townsend, in review: "Climatological trends influencing 
               forest phenology around western Lake Superior, USA." Submitted to J. Geophys. 
               Res. Atmos. on 5 April 2016.
           See also 'README.md', 'CITATION.txt', and 'ACKNOWLEDGEMENTS.txt' for more information.

USAGE: 'python process_NCEI_12.py 1984 2013 ./analyses'

PURPOSE: Use statistical test results to cluster ecoregions using a time-series similarity 
         measure that was defined in 'process_NCEI_11.py' and then examine the clusters for 
         subsets, overlaps, opportunities to merge/separate, singletons, etc. 
        
NOTE 1: The selected similarity measure uses p-values (either p or 1-p) from all three time 
        series statistical test results as analogous to Euclidean distance in x-y-z space. 
        See 'process_NCEI_11.py' and the reference listed above for more discussion on that 
        formulation.
 
NOTE 2: The clustering procedure is a multi-step process that may take some manual tweaking 
        on the user's part to get it right. There are many and more sophisticated ways to 
        cluster ecoregions according to arbitrary criteria. 

DEPENDENCIES: Some standard libraries/modules are required, as listed in 'import' lines
              The h5py module is required for handling of HDF5 files, used here for input

INPUT: Ecoregion time series statistical tests in '.h5' files from 'process_NCEI_11.py'

OUTPUT: Ecoregion clusters in '.txt' file

RUN TIME: ~5 sec

TO DO: Rewrite/rebuild this process for robustness and self-regulation with stable results
       -- ensure spatial contiguity of clustered ecoregions
       -- make <dist_lower_threshold> a dynamic variable, based on... ???
       -- use <dist_upper_threshold> to ensure dissimilar ecoregions don't get into same cluster
       -- use similarity as fuzzy likelihood of belonging to same cluster?
"""

import os, sys, datetime
import h5py as hdf
import numpy as np

def message(str):
    print str
    sys.stdout.flush()
    return    

def list_clusters(clusters,nclusters):
    for i in range(0,nclusters):
        message('cluster %d contains %d ecoregions: %s' % (i+1,len(clusters[i]),str(sorted(clusters[i]))))
    return

message(' ')
message('process_NCEI_12.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
# 
if len(sys.argv) < 4:
    message('input warning: no input analyses directory path indicated, using ./analyses')
    path = './analyses'
else:
    path = sys.argv[3]
# 
if len(sys.argv) < 3:
    year_begin = 1984
    year_end = 2013 
else:
    year_begin = int(sys.argv[1])      
    year_end = int(sys.argv[2])      
#
tsfname = '%s/%d-%d_ecoregion_timeseries.h5' % (path,year_begin,year_end)
message('getting ecoregion designations and variable names from %s' % tsfname)
h5infile = hdf.File(tsfname,'r')
ecoregion_IDs = np.copy(h5infile['ecoregion_IDs'])
necoregions = len(ecoregion_IDs)
varnames = np.copy(h5infile['varnames'])
nvars = len(varnames)
h5infile.close()
message(' ')
#
var_dist_matrix = np.zeros((nvars-1,necoregions,necoregions))
message('reading variable-specific distance matrices from %s' % tsfname)
h5infile = hdf.File(tsfname,'r')
for i in range(1,nvars):
    var = '%s_dist_matrix' % varnames[i]
    var_dist_matrix[i-1,:,:] = np.copy(h5infile[var])
    message('- %s' % varnames[i])
h5infile.close()
message(' ')
#
message('analyzing overall paired ecoregion distance measures')
dist_tuples = []
overall_dist_matrix = np.mean(var_dist_matrix, axis=0)
overall_dist_matrix_nan = np.where(overall_dist_matrix == 0.0, np.nan, overall_dist_matrix)
dist_mean = np.nanmean(overall_dist_matrix_nan)
dist_std = np.nanstd(overall_dist_matrix_nan)
dist_lower_threshold = dist_mean - 1.7 * dist_std
dist_upper_threshold = dist_mean + 1.0 * dist_std
message('- distance overall mean = %.3f  std = %.3f  lower threshold = %.3f  upper threshold = %.3f' % (dist_mean,dist_std,dist_lower_threshold,dist_upper_threshold))
for i in range(0,necoregions):
    for j in range(i,necoregions):
        if (overall_dist_matrix[i,j] < dist_lower_threshold) and (j != i):
            message('- similarity: ecoregions %s and %s are overall similar with D = %.3f' % (ecoregion_IDs[i],ecoregion_IDs[j],overall_dist_matrix[i,j]))
            dist_tuples.append((overall_dist_matrix[i,j],ecoregion_IDs[i],ecoregion_IDs[j]))
            #
message(' ')
#
dist_sorted = sorted(dist_tuples, key=lambda distance: distance[0])
message('sorted similar pairs')
for pair in dist_sorted:
    message('  %s' % str(pair))
message(' ')
#
message('clustering according to pair similarity')
clusters = []
nclusters = 0
for pair in dist_sorted:
    if clusters == 0:
        clusters.append([pair[1],pair[2]])
        nclusters = 1
    else:
        appended1 = 0
        appended2 = 0
        for cluster in clusters:
            if (pair[1] in cluster) and (pair[2] in cluster):
                pass
            elif (pair[1] not in cluster) and (pair[2] in cluster):
                cluster.append(pair[1])
                appended1 += 1
            elif (pair[1] in cluster) and (pair[2] not in cluster):
                cluster.append(pair[2])
                appended2 += 1
        if appended1 > 0:
            if appended1 > 1:
                message('- ecoregion %s appended to more than one cluster' % pair[1])
        if appended2 > 0:
            if appended2 > 1:
                message('- ecoregion %s appended to more than one cluster' % pair[2])
        if (appended1 == 0) and (appended2 == 0):
            clusters.append([pair[1],pair[2]])
            nclusters += 1
message(' ')
list_clusters(clusters,nclusters)
message(' ')
#
message('finding identical clusters')
to_remove = []
for i in range(0,nclusters):
    for j in range(i,nclusters):
        if j != i:
            c1len = len(clusters[i])
            c2len = len(clusters[j])
            culen = len(list(set(clusters[j]) | set(clusters[i])))
            if (c1len == culen) and (c2len == culen):
                message('- cluster %d is identical to cluster %d' % (i+1,j+1))
                to_remove.append(j)
message(' ')
to_remove = sorted(set(to_remove))
message('%d clusters to be removed: %s' % (len(to_remove),str(to_remove)))
for i in reversed(to_remove):
    del clusters[i]
    nclusters -= 1
message(' ')
list_clusters(clusters,nclusters)
message(' ')
#
message('finding subset clusters')
to_remove = []
for i in range(0,nclusters):
    for j in range(i,nclusters):
        if j != i:
            c1len = len(clusters[i])
            c2len = len(clusters[j])
            culen = len(list(set(clusters[j]) | set(clusters[i])))
            if (c1len == culen) and (c2len < culen):
                message('- cluster %d is a subset of cluster %d' % (j+1,i+1))
                to_remove.append(j)
            elif (c2len == culen) and (c1len < culen):
                message('- cluster %d is a subset of cluster %d' % (i+1,j+1))
                to_remove.append(i)
            # elif (c1len == culen) and (c2len == culen):
            #     message('- cluster %d is identical to cluster %d' % (i+1,j+1))
            #    to_remove.append(j)
message(' ')
to_remove = sorted(set(to_remove))
message('%d clusters to be removed: %s' % (len(to_remove),str(to_remove)))
for i in reversed(to_remove):
    del clusters[i]
    nclusters -= 1
message(' ')
list_clusters(clusters,nclusters)
message(' ')
#
message('finding combination possibilities based on overlap')
to_combine = []
for i in range(0,nclusters):
    for j in range(i,nclusters):
        if j != i:
            noverlap = len(list(set(clusters[i]) & set(clusters[j])))
            message('clusters %d and %d overlap on %d ecoregions' % (i+1,j+1,noverlap))
            if (noverlap > len(clusters[i])/2) and (noverlap > len(clusters[j])/2):
                to_combine.append([i, j])            
message('%d cluster pairs to be combined: %s' % (len(to_combine),str(to_combine)))
for pair in to_combine:
    clusters[pair[0]] = list(set(clusters[pair[0]]) | set(clusters[pair[1]]))
    del clusters[pair[1]]
    nclusters -= 1
message(' ')
list_clusters(clusters,nclusters)
message(' ')
#
message('finding split possibilities based on overlap')
for i in range(0,nclusters):
    for j in range(i,nclusters):
        if j != i:
            overlap = list(set(clusters[i]) & set(clusters[j]))
            noverlap = len(overlap)
            message('clusters %d and %d overlap on %d ecoregions' % (i+1,j+1,noverlap))
            if noverlap != 0:
                message('- splitting %s from clusters %d and %d' % (str(overlap),i+1,j+1))
                clusters[i] = list(set(clusters[i]) - set(overlap))
                clusters[j] = list(set(clusters[j]) - set(overlap))
                clusters.append(overlap)
                nclusters += 1
to_remove = []
for i in range(0,nclusters):
    if len(clusters[i]) == 0:
        to_remove.append(i)
for i in reversed(sorted(set(to_remove))):
    del clusters[i]
    nclusters -= 1
message(' ')
list_clusters(clusters,nclusters)
message(' ')
#
cluster_set = []
for i in range(0,nclusters):
    cluster_set = list(set(cluster_set) | set(clusters[i]))
message('full cluster set contains %d ecoregions: %s' % (len(cluster_set),str(sorted(cluster_set))))
missing = list(set(ecoregion_IDs) - set(cluster_set))
message('- missing %d ecoregions from total set: %s' % (len(missing),str(sorted(missing))))
message('-- adding each as singleton clusters')
if len(missing) > 0:
    for i in missing:
        clusters.append([i])
        nclusters += 1
message(' ')
list_clusters(clusters,nclusters)
message(' ')
#
clustersfname = '%s/%d-%d_ecoregion_clusters.txt' % (path,year_begin,year_end)
clustersf = open(clustersfname,'w')
for i in range(0,nclusters):
    line = '%d: ' % (i+1)
    if len(clusters[i]) == 1:
        line += '%s \n' % clusters[i][0]
    else:
        cluster_ecoregions = sorted(clusters[i])
        for j in range(0,len(cluster_ecoregions)-1):
            line += '%s, ' % cluster_ecoregions[j]
        line += '%s \n' % cluster_ecoregions[j+1]
    clustersf.write(line)
clustersf.close()
message('wrote ecoregion clusters to %s' % clustersfname)
message(' ')
#
message('process_NCEI_12.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end process_NCEI_12.py
