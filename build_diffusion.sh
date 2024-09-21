#!/bin/bash


. ./env_JUWELS_Stage2024.sh
rm -rf build

cmake -S src -B build
cmake --build ./build -- VERBOSE=1

