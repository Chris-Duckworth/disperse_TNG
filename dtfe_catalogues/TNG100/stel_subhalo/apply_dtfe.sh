#!/bin/bash

FILES=$DISPERSE/disperse_TNG/tracer_catalogues/TNG100/stel_subhalo/*.ascii

for file in $FILES
do
    name="$(basename $file)"
    echo $name
    name+=".NDnet_bash_readout.txt"
    #echo ./readout/$name
    $DISPERSE/bin/delaunay_3D $file > ./readout/$name
done
