#!/usr/bin/env bash

mkdir build
cd build
cmake ..
make -j4
cp ripser ..
cd ..
./ripser --benchmark_format=json --benchmark_out=result.json --benchmark_repetitions=3
