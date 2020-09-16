#!/bin/bash

export FW_SEARCH_PATH=.:${FW_SEARCH_PATH}

for file in $(ls 40.fcl)
do
art -n 1000 -c $file -o gen${file/fcl/root}
art -c readoutsimjob.fcl -o sim${file/fcl/root} gen${file/fcl/root}
art -c recojob.fcl -o reco${file/fcl/root} sim${file/fcl/root}
art -c anajob.fcl reco${file/fcl/root}
mv anatree.root anatree${file/fcl/root}
done
