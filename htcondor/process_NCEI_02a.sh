#!/bin/bash

tar -xzf python.tar.gz
export PATH=miniconda2/bin:$PATH
python process_NCEI_02a.py NCEI_WLS_19830101-20131231 /mnt/gluster/megarcia/WLS_Climatology/data
