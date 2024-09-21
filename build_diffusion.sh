#!/bin/bash


. ./env_JUWELS_Stage2024.sh

rm -rf build
mkdir build

cmake -S . -B build
cmake --build ./build -- VERBOSE=1

mv build/src/diffusion bin/
