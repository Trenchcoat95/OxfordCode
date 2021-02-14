#!/bin/sh

export FW_SEARCH_PATH=.:${FW_SEARCH_PATH}

art -n 1000 -c hem.fcl -o hemgen.root
art -c readoutsimjob.fcl -o hemsim.root hemgen.root
art -c recopatreccheatjob.fcl -o hemreco.root hemsim.root
art -c anajob.fcl hemreco.root
mv anatree.root hemanatree.root
