#!/bin/bash

for i in {1..50}; do
    echo "Running seed $i"


    ./program input_sub.i $i

done
