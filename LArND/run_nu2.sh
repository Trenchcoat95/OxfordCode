#! /usr/bin/env bash

####################

HORN=$1
NPER=$2
FIRST=$3
TEST=$4

if [ "${HORN}" != "FHC" ] && [ "${HORN}" != "RHC" ]; then
    echo "Invalid beam mode ${HORN}"
    echo "Must be FHC or RHC"
    kill -INT $$
fi

if [ "${NPER}" = "" ]; then
    echo "Number of events per job not specified, using 1000"
    NPER=1000
fi

if [ "${FIRST}" = "" ]; then
    echo "First run number not specified, using 0"
    FIRST=0
fi

CP="ifdh cp"
#CP="cp"
if [ "${TEST}" = "test" ]; then
    echo "In TEST mode, assuming interactive running"
    PROCESS=0
fi

MODE="neutrino"
if [ "${HORN}" = "RHC" ]; then
MODE="antineutrino"
fi

echo "Running gevgen for ${NPER} events in ${HORN} mode"
RNDSEED=$((${PROCESS}+${FIRST}))
#NEVENTS="-n ${NPER}"      # No. of events, -e XE16 for POT

GEOMETRY="nd_hall_only_kloe"
TOPVOL="volMainDet_3DST"
SPOT=10
#SEEDING=$(($RANDOM%$SPOT+314159*$PROCESS ))
#NEVENTS="-e 5E16"
NEVENTS="-n 200"
#NEVENTS = 200
PROD="PROD1"

USERDIR="/pnfs/dune/persistent/users/mtanaz/test_12"
OUTDIR="/pnfs/dune/persistent/users/mtanaz/test_12"

OUTFLAG="3DST"

####################

## Setup UPS and required products

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup genie        v2_12_10c  -q e15:prof
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
#setup dk2nu        v01_05_01b -q e15:prof
setup ifdhc

####################

## Run GENIE 

export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux_DUNE_ND.xml"

gevgen_fnal \
    -f /pnfs/dune/persistent/users/mtanaz/gsimple_subdetectors/${MODE}/gsimple*.root,DUNEND_ND_NEAR_FID1 \
    -g ${GEOMETRY}.gdml \
    -t ${TOPVOL} \
    -L cm -D g_cm3 \
    ${NEVENTS} \
    --seed ${RNDSEED} \
    -r ${RNDSEED} \
    -o ${MODE} \
    --message-thresholds Messenger_production.xml \
    --cross-sections ${GENIEXSECPATH}/gxspl-FNALsmall.xml \
    --event-record-print-level 0 \
    --event-generator-list Default+CCMEC