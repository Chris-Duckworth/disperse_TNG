#!/bin/bash

FILES=$DISPERSE/disperse_TNG/morse_smale_catalogues/TNG100/stel_subhalo/*NDskl

for file in $FILES
do
    name="$(basename $file)"
    echo $name
    name+=".crits_bash_readout.txt"
    echo ./readout/$name
    $DISPERSE/bin/skelconv $file -breakdown -smooth 1 -to crits_ascii > ./readout/$name
done
