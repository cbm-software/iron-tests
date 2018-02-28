#!/bin/bash

echo "compiling and running example $(pwd)"

folder=$1

echo "  running $folder"
echo "    2D tests"
mkdir -p results/current_run/d2_i1_s0 && ./bin/example 2 1 0
mkdir -p results/current_run/d2_i2_s0 && ./bin/example 2 2 0
mkdir -p results/current_run/d2_i1_s1 && ./bin/example 2 1 1
mkdir -p results/current_run/d2_i2_s1 && ./bin/example 2 2 1
echo "    3D tests"
mkdir -p results/current_run/d3_i1_s0 && ./bin/example 3 1 0
mkdir -p results/current_run/d3_i2_s0 && ./bin/example 3 2 0
mkdir -p results/current_run/d3_i1_s1 && ./bin/example 3 1 1
mkdir -p results/current_run/d3_i2_s1 && ./bin/example 3 2 1
