#!/bin/bash

rm -r build

cmake -S src -B build
cmake --build ./build -- VERBOSE=1

