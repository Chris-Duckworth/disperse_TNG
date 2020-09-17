#!/bin/bash

FILES=$DISPERSE/disperse_TNG/morse_smale_catalogues/TNG300/stel_subhalo/*.NDskl

for file in $FILES
do
    name="$(basename $file)"
    echo $name
    name+=".segs_bash_readout.txt"
    echo ./readout/$name
    $DISPERSE/bin/skelconv $file -breakdown -smooth 1 -to segs_ascii > ./readout/$name
done
