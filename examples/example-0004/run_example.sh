#!/bin/bash

echo "compiling and running example $(pwd)"

folder=$1

echo "  running $folder"
echo "    2D tests"
mkdir -p results/current_run/l4x2x0_n4x2x0_i1_s0 && ./bin/example 4 2 0 1 0
mkdir -p results/current_run/l4x2x0_n8x4x0_i1_s0 && ./bin/example 8 4 0 1 0
mkdir -p results/current_run/l4x2x0_n2x1x0_i2_s0 && ./bin/example 2 1 0 2 0
mkdir -p results/current_run/l4x2x0_n4x2x0_i2_s0 && ./bin/example 4 2 0 2 0
mkdir -p results/current_run/l4x2x0_n8x4x0_i2_s0 && ./bin/example 8 4 0 2 0
mkdir -p results/current_run/l4x2x0_n4x2x0_i1_s1 && ./bin/example 4 2 0 1 1
mkdir -p results/current_run/l4x2x0_n8x4x0_i1_s1 && ./bin/example 8 4 0 1 1
mkdir -p results/current_run/l4x2x0_n2x1x0_i2_s1 && ./bin/example 2 1 0 2 1
mkdir -p results/current_run/l4x2x0_n4x2x0_i2_s1 && ./bin/example 4 2 0 2 1
mkdir -p results/current_run/l4x2x0_n8x4x0_i2_s1 && ./bin/example 8 4 0 2 1
mkdir -p results/current_run/l4x2x0_n100x50x0_i1_s0 && ./bin/example 100 50 0 1 0
mkdir -p results/current_run/l4x2x0_n100x50x0_i1_s1 && ./bin/example 100 50 0 1 1
mkdir -p results/current_run/l4x2x0_n100x50x0_i2_s0 && ./bin/example 100 50 0 2 0
mkdir -p results/current_run/l4x2x0_n100x50x0_i2_s1 && ./bin/example 100 50 0 2 1
