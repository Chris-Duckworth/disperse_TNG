#!/bin/bash

FILES=$DISPERSE/disperse_TNG/tracer_catalogues/TNG100/stel_subhalo/*.ascii

for file in $FILES
do
    echo $(basename $file)
done
