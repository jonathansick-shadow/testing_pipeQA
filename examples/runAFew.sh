#!/usr/bin/env bash

# lsstSim
LSST="pipeQa.py -e -d -C lsstSim -b ccd -v 890104911 -r 1,1 -c 1,1 krughoff_S12_lsstsim_u_krughoff_2012_0706_183555"

# coadd
COADD="pipeQa.py -e -d -z 'all' -T goodSeeing -C coadd -b ccd -k -v '0-250,0-r' daues_S12_sdss_u_s2012prod_sd1000"

# sdss
SDSS="pipeQa.py -e -F -d -C sdss -b ccd -k -v '1056-200' -c g daues_S12_sdss_u_s2012prod_sd1000 "


for CMD in "$LSST $COADD $SDSS"
do
    echo '#############################################################################'
    echo $CMD
    echo '#############################################################################'
    $CMD
done
