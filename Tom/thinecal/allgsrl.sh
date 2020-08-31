#!/bin/sh

export FW_SEARCH_PATH=.:${FW_SEARCH_PATH}
art -n 1000 -c lem.fcl -o lemgen.root
art -c readoutsimjob.fcl -o lemsim.root lemgen.root
art -c recopatreccheatjob.fcl -o lemreco.root lemsim.root
art -c anajob.fcl lemreco.root
mv anatree.root lemanatree.root

