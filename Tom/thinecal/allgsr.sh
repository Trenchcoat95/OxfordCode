#!/bin/sh

export FW_SEARCH_PATH=.:${FW_SEARCH_PATH}
art -n 1000 -c lem.fcl -o lemgen.root
art -c readoutsimjob.fcl -o lemsim.root lemgen.root
art -c recojob.fcl -o lemreco.root lemsim.root
art -c anajob.fcl lemreco.root
mv anatree.root lemanatree.root

art -n 1000 -c hem.fcl -o hemgen.root
art -c readoutsimjob.fcl -o hemsim.root hemgen.root
art -c recojob.fcl -o hemreco.root hemsim.root
art -c anajob.fcl hemreco.root
mv anatree.root hemanatree.root
