#!/bin/bash

echo "compiling and running example $(pwd)"

folder=$1
kind=$2

echo "  running $folder"
#mkdir -p results/current_run && ./bin/example
#mkdir -p results/current_run && ./bin/example 160 120 0 32 24 0 1 0 10000 0.3 0.025
# These are the linear elasticity tests for simple shear (displacement driven)
# 2D Tests
mkdir -p results/current_run/l160x120x000_n08x06x00_i1_s0 && ./bin/example 160 120 0 8 6 0 1 0 10000 0.3 0.05
mkdir -p results/current_run/l160x120x000_n16x12x00_i1_s0 && ./bin/example 160 120 0 16 12 0 1 0 10000 0.3 0.05
mkdir -p results/current_run/l160x120x000_n32x24x00_i1_s0 && ./bin/example 160 120 0 32 24 0 1 0 10000 0.3 0.05
mkdir -p results/current_run/l160x120x000_n08x06x00_i2_s0 && ./bin/example 160 120 0 8 6 0 2 0 10000 0.3 0.05
mkdir -p results/current_run/l160x120x000_n16x12x00_i2_s0 && ./bin/example 160 120 0 16 12 0 2 0 10000 0.3 0.05
mkdir -p results/current_run/l160x120x000_n32x24x00_i2_s0 && ./bin/example 160 120 0 32 24 0 2 0 10000 0.3 0.05
# 3D Tests
#echo "    3D tests"
mkdir -p results/current_run/l160x120x120_n08x06x06_i1_s0 && ./bin/example 160 120 120 8 6 6 1 0 10000 0.3 0.05
mkdir -p results/current_run/l160x120x120_n16x12x12_i1_s0 && ./bin/example 160 120 120 16 12 12 1 0 10000 0.3 0.05
if [ "$kind" == "big" ]; then
    mkdir -p results/current_run/l160x120x120_n32x24x24_i1_s0 && ./bin/example 160 120 120 32 24 24 1 0 10000 0.3 0.05
    mkdir -p results/current_run/l160x120x120_n08x06x06_i2_s0 && ./bin/example 160 120 120 8 6 6 2 0 10000 0.3 0.05
    mkdir -p results/current_run/l160x120x120_n16x12x12_i2_s0 && ./bin/example 160 120 120 16 12 12 2 0 10000 0.3 0.05
    mkdir -p results/current_run/l160x120x120_n32x24x24_i2_s0 && ./bin/example 160 120 120 32 24 24 2 0 10000 0.3 0.05
fi
#
# Write out x-displacements into a Text-File (-> Gives a list as ASCII-File [1,last_node])
#
echo "    Processing files"
# 2d
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n08x06x00_i1_s0/ && cd results/current_run/l160x120x000_n08x06x00_i1_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n08x06x00_i1_s0_shear_disp.txt && cd ../../../
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n16x12x00_i1_s0/ && cd results/current_run/l160x120x000_n16x12x00_i1_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n16x12x00_i1_s0_shear_disp.txt && cd ../../../
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n32x24x00_i1_s0/ && cd results/current_run/l160x120x000_n32x24x00_i1_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n32x24x00_i1_s0_shear_disp.txt && cd ../../../
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n08x06x00_i2_s0/ && cd results/current_run/l160x120x000_n08x06x00_i2_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n08x06x00_i2_s0_shear_disp.txt && cd ../../../
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n16x12x00_i2_s0/ && cd results/current_run/l160x120x000_n16x12x00_i2_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n16x12x00_i2_s0_shear_disp.txt && cd ../../../
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n32x24x00_i2_s0/ && cd results/current_run/l160x120x000_n32x24x00_i2_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n32x24x00_i2_s0_shear_disp.txt && cd ../../../
# 3d
cp results/get_displacement_3d.pl results/current_run/l160x120x120_n08x06x06_i1_s0/ && cd results/current_run/l160x120x120_n08x06x06_i1_s0/ && ./get_displacement_3d.pl && mv my_displacement_5.txt l160x120x120_n08x06x06_i1_s0_shear_disp.txt && cd ../../../
cp results/get_displacement_3d.pl results/current_run/l160x120x120_n16x12x12_i1_s0/ && cd results/current_run/l160x120x120_n16x12x12_i1_s0/ && ./get_displacement_3d.pl && mv my_displacement_5.txt l160x120x120_n16x12x12_i1_s0_shear_disp.txt && cd ../../../
if [ "$kind" == "big" ]; then
    cp results/get_displacement_3d.pl results/current_run/l160x120x120_n32x24x24_i1_s0/ && cd results/current_run/l160x120x120_n32x24x24_i1_s0/ && ./get_displacement_3d.pl && mv my_displacement_5.txt l160x120x120_n32x24x24_i1_s0_shear_disp.txt && cd ../../../
    cp results/get_displacement_3d.pl results/current_run/l160x120x120_n08x06x06_i2_s0/ && cd results/current_run/l160x120x120_n08x06x06_i2_s0/ && ./get_displacement_3d.pl && mv my_displacement_5.txt l160x120x120_n08x06x06_i2_s0_shear_disp.txt && cd ../../../
    cp results/get_displacement_3d.pl results/current_run/l160x120x120_n16x12x12_i2_s0/ && cd results/current_run/l160x120x120_n16x12x12_i2_s0/ && ./get_displacement_3d.pl && mv my_displacement_5.txt l160x120x120_n16x12x12_i2_s0_shear_disp.txt && cd ../../../
    cp results/get_displacement_3d.pl results/current_run/l160x120x120_n32x24x24_i2_s0/ && cd results/current_run/l160x120x120_n32x24x24_i2_s0/ && ./get_displacement_3d.pl && mv my_displacement_5.txt l160x120x120_n32x24x24_i2_s0_shear_disp.txt && cd ../../../
fi
