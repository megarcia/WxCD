#!/bin/bash

tar -xzf python.tar.gz
export PATH=miniconda2/bin:$PATH
python process_NCEI_02b.py NLCD_2011_WLS_UTM15N NCEI_WLS_$1 /mnt/gluster/megarcia/WLS_Climatology/grids 500 RBF 1
