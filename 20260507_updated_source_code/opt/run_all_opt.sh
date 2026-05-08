#!/bin/bash

for i in {1..50}; do
    echo "Running seed $i"


    ./program input_opt.i $i

done
