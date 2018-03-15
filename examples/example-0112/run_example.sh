#!/bin/bash

echo "compiling and running example $(pwd)"

folder=$1
kind=$2

echo "  running $folder"
#mkdir -p results/current_run && ./bin/example
#mkdir -p results/current_run && ./bin/example 160 120 0 32 24 0 1 0 10000 0.3 0.025
# These are the linear elasticity tests for simple shear (force driven)
# 2D Tests
echo "    2D tests"
mkdir -p results/current_run/l160x120x000_n08x06x00_i1_s0 && ./bin/example 160 120 0 8 6 0 1 0 10000 0.3 60000
mkdir -p results/current_run/l160x120x000_n16x12x00_i1_s0 && ./bin/example 160 120 0 16 12 0 1 0 10000 0.3 60000
mkdir -p results/current_run/l160x120x000_n32x24x00_i1_s0 && ./bin/example 160 120 0 32 24 0 1 0 10000 0.3 60000
mkdir -p results/current_run/l160x120x000_n08x06x00_i2_s0 && ./bin/example 160 120 0 8 6 0 2 0 10000 0.3 60000
mkdir -p results/current_run/l160x120x000_n16x12x00_i2_s0 && ./bin/example 160 120 0 16 12 0 2 0 10000 0.3 60000
mkdir -p results/current_run/l160x120x000_n32x24x00_i2_s0 && ./bin/example 160 120 0 32 24 0 2 0 10000 0.3 60000
# 3D Tests
echo "    3D tests"
mkdir -p results/current_run/l160x120x120_n08x06x06_i1_s0 && ./bin/example 160 120 120 8 6 6 1 0 10000 0.3 7200000
mkdir -p results/current_run/l160x120x120_n16x12x12_i1_s0 && ./bin/example 160 120 120 16 12 12 1 0 10000 0.3 7200000
if [ "$kind" == "big" ]; then
    mkdir -p results/current_run/l160x120x120_n32x24x24_i1_s0 && ./bin/example 160 120 120 32 24 24 1 0 10000 0.3 7200000
    echo ">>>WARNING: Interpolation type 2 not tested for 3D"
fi
#
# Write out x-displacements into a Text-File (-> Gives a list as ASCII-File [1,last_node])
#
# 2d tests
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n08x06x00_i1_s0/ && cd results/current_run/l160x120x000_n08x06x00_i1_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n08x06x00_i1_s0_shear_force.txt &&
cd ../../../
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n16x12x00_i1_s0/ && cd results/current_run/l160x120x000_n16x12x00_i1_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n16x12x00_i1_s0_shear_force.txt &&
cd ../../../
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n32x24x00_i1_s0/ && cd results/current_run/l160x120x000_n32x24x00_i1_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n32x24x00_i1_s0_shear_force.txt &&
cd ../../../
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n08x06x00_i2_s0/ && cd results/current_run/l160x120x000_n08x06x00_i2_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n08x06x00_i2_s0_shear_force.txt &&
cd ../../../
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n16x12x00_i2_s0/ && cd results/current_run/l160x120x000_n16x12x00_i2_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n16x12x00_i2_s0_shear_force.txt &&
cd ../../../
cp results/get_displacement_2d.pl results/current_run/l160x120x000_n32x24x00_i2_s0/ && cd results/current_run/l160x120x000_n32x24x00_i2_s0/ && ./get_displacement_2d.pl && mv my_displacement_5.txt l160x120x000_n32x24x00_i2_s0_shear_force.txt &&
cd ../../../
# 3d tests
cp results/get_displacement_3d.pl results/current_run/l160x120x120_n08x06x06_i1_s0/ && cd results/current_run/l160x120x120_n08x06x06_i1_s0/ && ./get_displacement_3d.pl && mv my_displacement_5.txt l160x120x120_n08x06x00_i1_s0_shear_force.txt &&
cd ../../../
cp results/get_displacement_3d.pl results/current_run/l160x120x120_n16x12x12_i1_s0/ && cd results/current_run/l160x120x120_n16x12x12_i1_s0/ && ./get_displacement_3d.pl && mv my_displacement_5.txt l160x120x120_n16x12x12_i1_s0_shear_force.txt &&
cd ../../../
if [ "$kind" == "big" ]; then
    cp results/get_displacement_3d.pl results/current_run/l160x120x120_n32x24x24_i1_s0/ && cd results/current_run/l160x120x120_n32x24x24_i1_s0/ && ./get_displacement_3d.pl && mv my_displacement_5.txt l160x120x120_n32x24x24_i1_s0_shear_force.txt &&
    cd ../../../
fi
