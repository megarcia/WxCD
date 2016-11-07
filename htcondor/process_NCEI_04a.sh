#!/bin/bash

tar -xzf python.tar.gz
export PATH=miniconda2/bin:$PATH
python process_NCEI_04a.py $1 $2 /mnt/gluster/megarcia/WLS_Climatology/grids
