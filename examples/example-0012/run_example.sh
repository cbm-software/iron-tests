#!/bin/bash

echo "compiling and running example $(pwd)"

folder=$1

echo "  running $folder"
echo "    2D tests"
mkdir -p results/current_run/l2x1x0_n2x1x0_i1_s0 && ./bin/example 2.0 1.0 0.0 2 1 0 1 0 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n4x2x0_i1_s0 && ./bin/example 2.0 1.0 0.0 4 2 0 1 0 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n8x4x0_i1_s0 && ./bin/example 2.0 1.0 0.0 8 4 0 1 0 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n2x1x0_i2_s0 && ./bin/example 2.0 1.0 0.0 2 1 0 2 0 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n4x2x0_i2_s0 && ./bin/example 2.0 1.0 0.0 4 2 0 2 0 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n8x4x0_i2_s0 && ./bin/example 2.0 1.0 0.0 8 4 0 2 0 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n2x1x0_i1_s1 && ./bin/example 2.0 1.0 0.0 2 1 0 1 1 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n4x2x0_i1_s1 && ./bin/example 2.0 1.0 0.0 4 2 0 1 1 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n8x4x0_i1_s1 && ./bin/example 2.0 1.0 0.0 8 4 0 1 1 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n2x1x0_i2_s1 && ./bin/example 2.0 1.0 0.0 2 1 0 2 1 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n4x2x0_i2_s1 && ./bin/example 2.0 1.0 0.0 4 2 0 2 1 2 3 0 0 0 0
mkdir -p results/current_run/l2x1x0_n8x4x0_i2_s1 && ./bin/example 2.0 1.0 0.0 8 4 0 2 1 2 3 0 0 0 0
echo "    3D tests"
mkdir -p results/current_run/l2x1x1_n2x1x1_i1_s0 && ./bin/example 2.0 1.0 1.0 2 1 1 1 0 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n4x2x2_i1_s0 && ./bin/example 2.0 1.0 1.0 4 2 2 1 0 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n8x4x4_i1_s0 && ./bin/example 2.0 1.0 1.0 8 4 4 1 0 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n2x1x1_i2_s0 && ./bin/example 2.0 1.0 1.0 2 1 1 2 0 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n4x2x2_i2_s0 && ./bin/example 2.0 1.0 1.0 4 2 2 2 0 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n8x4x4_i2_s0 && ./bin/example 2.0 1.0 1.0 8 4 4 2 0 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n2x1x1_i1_s1 && ./bin/example 2.0 1.0 1.0 2 1 1 1 1 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n4x2x2_i1_s1 && ./bin/example 2.0 1.0 1.0 4 2 2 1 1 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n8x4x4_i1_s1 && ./bin/example 2.0 1.0 1.0 8 4 4 1 1 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n2x1x1_i2_s1 && ./bin/example 2.0 1.0 1.0 2 1 1 2 1 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n4x2x2_i2_s1 && ./bin/example 2.0 1.0 1.0 4 2 2 2 1 2 3 7 0 0 0
mkdir -p results/current_run/l2x1x1_n8x4x4_i2_s1 && ./bin/example 2.0 1.0 1.0 8 4 4 2 1 2 3 7 0 0 0
