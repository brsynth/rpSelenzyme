#!/bin/bash

echo "Test in normal mode"
python3 ../src/main.py \
  in/rhea15870.rxn \
  /home/data \
  out

# echo "Test in db mode"
# python3 ../src/main.py \
#   -sm db \
#   in/rhea15870.rxn \
#   /home/data \
#   out
