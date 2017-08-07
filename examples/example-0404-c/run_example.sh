#!/bin/bash
# test cases for example-0404 (1D problem with Hodgkin-Huxley)

echo "compiling and running example $(pwd)"

folder=$1   # this should be 'debug' or 'release'
kind=$2     # this should be 'big' or <anything else>='fast'

mkdir -p $folder

echo "  compiling $folder"
cd $folder
cmake -DCMAKE_BUILD_TYPE=$folder -DOPENCMISS_BUILD_TYPE=$folder ..
make
cd ..
echo "  running $folder"
# <number elements X> <interpolation type> <solver type> <PDE step size> <stop time> <output frequency> <CellML Model URL> <slow-twitch> <ODE time-step>

t_end=4.0
# for 'big' target set t_end to 4
if [ "$kind" == "big" ]; then
  t_end=10.0
fi

mkdir -p results/current_run/l1x1_n64_i2_s0_01 && ./$folder/src/example 64 2 0 0.01 $t_end 5 hodgkin_huxley_1952.cellml F 0.01 && mv *.ex* results/current_run/l1x1_n64_i2_s0_01
mkdir -p results/current_run/l1x1_n64_i2_s0_005 && ./$folder/src/example 64 2 0 0.005 $t_end 10 hodgkin_huxley_1952.cellml F 0.005 && mv *.ex* results/current_run/l1x1_n64_i2_s0_005
mkdir -p results/current_run/l1x1_n64_i2_s0_001 && ./$folder/src/example 64 2 0 0.001 $t_end 50 hodgkin_huxley_1952.cellml F 0.001 && mv *.ex* results/current_run/l1x1_n64_i2_s0_001
mkdir -p results/current_run/l1x1_n64_i2_s0_0005 && ./$folder/src/example 64 2 0 0.0005 $t_end 100 hodgkin_huxley_1952.cellml F 0.0005 && mv *.ex* results/current_run/l1x1_n64_i2_s0_0005
mkdir -p results/current_run/l1x1_n64_i2_s0_00025 && ./$folder/src/example 64 2 0 0.00025 $t_end 200 hodgkin_huxley_1952.cellml F 0.00025 && mv *.ex* results/current_run/l1x1_n64_i2_s0_00025


