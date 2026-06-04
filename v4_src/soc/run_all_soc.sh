#!/bin/bash

for i in {1..100}; do
    echo "Running seed $i"


./program input_soc.i $i

done
