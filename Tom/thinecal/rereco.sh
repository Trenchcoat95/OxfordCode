#!/bin/sh

export FW_SEARCH_PATH=.:${FW_SEARCH_PATH}

art -c recojob.fcl -o hemreco.root hemsim.root
art -c anajob.fcl hemreco.root
mv anatree.root hemanatree.root

art -c recojob.fcl -o lemreco.root lemsim.root
art -c anajob.fcl lemreco.root
mv anatree.root lemanatree.root

