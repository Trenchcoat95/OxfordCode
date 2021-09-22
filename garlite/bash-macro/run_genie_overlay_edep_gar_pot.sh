#!/bin/bash

POT=150E13  #1 spill worth is 7.5E13
SPILLS=20
SEED=${RANDOM}
MODE="neutrino"
PROCESSO="Default+CCMEC"
FCLLOC="/pnfs/dune/scratch/users/battisti/fcl_bare"
GDML_SHORT="nd_hall_mpd_lar_only"
GEONAME="nd_hall_mpd_lar_DayOne_SPY_v3_wMuID"
GEO="/pnfs/dune/scratch/users/battisti/garlite/Geom/${GDML_SHORT}/${GEONAME}.gdml"
VOLUME="volArgonCubeDetector"
PROFILE="spill_profile.root"
PROFILEDIR="/pnfs/dune/scratch/users/battisti/garlite/Spill/production/Profile/${PROFILE}"
FLUXFILES="/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/Flux/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017/neutrino/gsimpleND/gsimple*.root"
LOC="DUNE_ND_HALL"

cd ${_CONDOR_SCRATCH_DIR}

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh

setup genie        v2_12_10c  -q e15:prof
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup dk2nugenie   v01_06_01f -q e15:prof

setup ifdhc
export IFDH_GRIDFTP_EXTRA="-st 10" #set ifdh cp stall timeout to 10 sec
export IFDH_CP_MAXRETRIES=3

ifdh cp ${GEO} ${GEONAME}.gdml
ifdh cp ${PROFILEDIR} ${PROFILE}
ifdh cp ${FCLLOC}/conversionedepjob.fcl conversionedepjob${SEED}.fcl
ifdh cp ${FCLLOC}/digi_reco_lite_spill.fcl digi_reco_lite_spill${SEED}.fcl
ifdh cp ${FCLLOC}/anajob_dayone_spill.fcl anajob_dayone_spill${SEED}.fcl


gevgen_fnal \
    -f ${FLUXFILES},${LOC} \
    -g ${GEONAME}.gdml \
    -t ${VOLUME} \
    -L cm \
    -D g_cm3 \
    -e ${POT} \
    --SEED ${SEED} \
    -r ${SEED} \
    -o ${MODE}.${GEONAME}.${VOLUME} \
    --message-thresholds Messenger_production.xml \
    --cross-sections ${GENIEXSECPATH}/gxspl-FNALsmall.xml \
    --event-record-print-level 0 \
    --event-generator-list ${PROCESSO}


ifdh cp -D ${MODE}.${GEONAME}.${VOLUME}.${SEED}.ghep.root  /pnfs/dune/scratch/users/battisti/garlite/Spill/production/Genie/

#######Now the overlay


${CONDOR_DIR_INPUT}/OverlayGenie/build/bin/overlay_genie \
    --source=${MODE}.${GEONAME}.${VOLUME}.${SEED}.ghep.root,gtree,2718,10,y,poisson,50 \
    --output=${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.ghep.root,gtree \
    --time_hist=${PROFILE},spill_profile \
    --nspills=${SPILLS}

ifdh cp -D ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.ghep.root  /pnfs/dune/scratch/users/battisti/garlite/Spill/production/Overlay/

#############Convert to gtracker

gntpc -i ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.ghep.root -f rootracker \
      -o ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.gtracker.root

ifdh cp -D ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.gtracker.root  /pnfs/dune/scratch/users/battisti/garlite/Spill/production/Overlay/

#############Propagate with edepsim

###clean setups

unsetup genie        v2_12_10c  -q e15:prof
unsetup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
unsetup genie_phyopt v2_12_10   -q dkcharmtau
unsetup dk2nugenie   v01_06_01f -q e15:prof

###create mac file

printf "/generator/kinematics/set rooTracker" > edepmac${SEED}.mac
printf "\n/generator/kinematics/rooTracker/input ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.gtracker.root" >> edepmac${SEED}.mac
printf "\n/generator/kinematics/rooTracker/tree gRooTracker" >> edepmac${SEED}.mac
printf "\n/generator/position/set free \n/generator/time/set free \n/generator/count/fixed/number 10000 \n/generator/count/set fixed \n/generator/add" >> edepmac${SEED}.mac

#ifdh cp -D edepmac${SEED}.mac  /pnfs/dune/scratch/users/battisti/garlite/Spill/production/Edep/


###setup and run edepsim

setup edepsim v3_0_1 -q e19:prof

edep-sim -C -g ${GEONAME}.gdml \
            -o ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.edep.root \
            -u \
            -e ${SPILLS} \
            edepmac${SEED}.mac

ifdh cp -D ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.edep.root  /pnfs/dune/scratch/users/battisti/garlite/Spill/production/Edep/

unsetup edepsim v3_0_1 -q e19:prof

#######garsoft

#####garsoft setup

#cd ${_CONDOR_SCRATCH_DIR}

#source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
#source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh
unsetup ifdhc

export MRB_PROJECT=garsoft
export MRB_PROJECT_VERSION=develop
export MRB_QUALS=e20:prof
export MRB_TOP=${CONDOR_DIR_INPUT}/garsoft
export MRB_SOURCE=${CONDOR_DIR_INPUT}/garsoft/srcs
export MRB_BUILDDIR=${CONDOR_DIR_INPUT}/garsoft/build_slf7.x86_64
export MRB_INSTALL=${CONDOR_DIR_INPUT}/garsoft/localProducts_garsoft_develop_e20_prof

export PRODUCTS=${CONDOR_DIR_INPUT}/garsoft/localProducts_garsoft_develop_e20_prof:${CONDOR_DIR_INPUT}/garsoft/localProducts_garsoft_develop_e20_prof:/cvmfs/dune.opensciencegrid.org/products/dune:/cvmfs/larsoft.opensciencegrid.org/products:/cvmfs/fermilab.opensciencegrid.org/products/common/db

cd ${MRB_BUILDDIR}

#mrbsetenv
mrbslp

cd ${_CONDOR_SCRATCH_DIR}


#####edepconvert


printf "\nphysics.producers.edepconvert.EDepSimFile: \"${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.edep.root\"" >> conversionedepjob${SEED}.fcl
printf "\nphysics.producers.edepconvert.GhepFile: \"${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.ghep.root\"" >> conversionedepjob${SEED}.fcl
printf "\nphysics.producers.edepconvert.OverlayFile: true" >> conversionedepjob${SEED}.fcl
printf "\nservices.GeometryGAr.GDML: \"${GEONAME}.gdml\"" >> conversionedepjob${SEED}.fcl
printf "\nservices.GeometryGAr.ROOT: \"${GEONAME}.gdml\"" >> conversionedepjob${SEED}.fcl

#ifdh cp -D conversionedepjob${SEED}.fcl  ${FCLLOC}/

art -n ${SPILLS} -c conversionedepjob${SEED}.fcl -o ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.edepconv.root

ifdh cp -D ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.edepconv.root  /pnfs/dune/scratch/users/battisti/garlite/Spill/conversion/

########digireco

printf "\nphysics.producers.dayone.InputLabel: \"edepconvert\"" >> digi_reco_lite_spill${SEED}.fcl
printf "\nservices.GeometryGAr.GDML: \"${GEONAME}.gdml\"" >> digi_reco_lite_spill${SEED}.fcl
printf "\nservices.GeometryGAr.ROOT: \"${GEONAME}.gdml\"" >> digi_reco_lite_spill${SEED}.fcl

#ifdh cp -D digi_reco_lite_spill${SEED}.fcl  ${FCLLOC}/

art -n ${SPILLS} -c digi_reco_lite_spill${SEED}.fcl ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.edepconv.root -o ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.reco.root

ifdh cp -D ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.reco.root  /pnfs/dune/scratch/users/battisti/garlite/Reco/

##########ana

printf "\nphysics.analyzers.anatree.GEANTLabel: \"edepconvert\"" >> anajob_dayone_spill${SEED}.fcl
printf "\nservices.GeometryGAr.GDML: \"${GEONAME}.gdml\"" >> anajob_dayone_spill${SEED}.fcl
printf "\nservices.GeometryGAr.ROOT: \"${GEONAME}.gdml\"" >> anajob_dayone_spill${SEED}.fcl

#ifdh cp -D anajob_dayone_spill${SEED}.fcl  ${FCLLOC}/

art -n ${SPILLS} -c anajob_dayone_spill${SEED}.fcl ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.reco.root -T ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.ana.root

ifdh cp -D ${MODE}.${GEONAME}.${VOLUME}.${SEED}.Overlay.ana.root  /pnfs/dune/scratch/users/battisti/garlite/Ana/

