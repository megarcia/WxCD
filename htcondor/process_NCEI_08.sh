#!/bin/bash

tar -xzf python.tar.gz
export PATH=miniconda2/bin:$PATH
python process_NCEI_08.py EPA_L4_Ecoregions_WLS_UTM15N /mnt/gluster/megarcia/WLS_Climatology/data
