#!/bin/bash

FILES=$DISPERSE/disperse_TNG/dtfe_catalogues/TNG100/stel_subhalo/*.NDnet

for file in $FILES
do
    name="$(basename $file)"
    echo $name
    name+=".mse_bash_readout.txt"
    echo ./readout/$name
    $DISPERSE/bin/mse $file -nsig 4 -upSkl > ./readout/$name
done
