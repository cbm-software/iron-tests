#!/bin/bash
# test cases for example-0402 (2D problem with Noble's 1998 improved guinea-pig ventricular cell model)

echo "compiling and running example $(pwd)"

folder=$1   # this should be 'debug' or 'release'
kind=$2     # this should be 'big' or <anything else>='fast'

echo "  running $folder"
# <number elements X> <number elements Y> <interpolation type> <solver type> <PDE step size> <stop time> <output frequency> <CellML Model URL> <slow-twitch> <ODE time-step>

output_frequency=10
t_end=1.0

# for 'big' target output every timestep and set t_end to 3
if [ "$kind" == "big" ]; then
  output_frequency=1
  t_end=3.0
fi

mkdir -p results/current_run/l1x1_n24x24_i1_s0_p1 && ./bin/example 24 24 1 0 0.005 $t_end $output_frequency n98.xml F 0.0001 && mv *.ex* results/current_run/l1x1_n24x24_i1_s0_p1
mkdir -p results/current_run/l1x1_n24x24_i1_s1_p1 && ./bin/example 24 24 1 0 0.005 $t_end $output_frequency n98.xml F 0.005  && mv *.ex* results/current_run/l1x1_n24x24_i1_s1_p1
mkdir -p results/current_run/l1x1_n10x10_i1_s0_p1 && ./bin/example 10 10 1 0 0.005 $t_end $output_frequency n98.xml F 0.0001 && mv *.ex* results/current_run/l1x1_n10x10_i1_s0_p1
mkdir -p results/current_run/l1x1_n24x24_i1_s0_p2 && mpirun -n 2 ./bin/example 24 24 1 0 0.005 $t_end $output_frequency n98.xml F 0.0001 && mv *.ex* results/current_run/l1x1_n24x24_i1_s0_p2
mkdir -p results/current_run/l1x1_n24x24_i1_s0_p8 && mpirun -n 8 ./bin/example 24 24 1 0 0.005 $t_end $output_frequency n98.xml F 0.0001 && mv *.ex* results/current_run/l1x1_n24x24_i1_s0_p8
mkdir -p results/current_run/l1x1_n2x2_i1_s0_p1 && ./bin/example 2 2 1 0 0.005 $t_end $output_frequency n98.xml F 0.0001 && mv *.ex* results/current_run/l1x1_n2x2_i1_s0_p1
mkdir -p results/current_run/l1x1_n2x2_i1_s0_p2 && mpirun -n 2 ./bin/example 2 2 1 0 0.005 $t_end $output_frequency n98.xml F 0.0001 && mv *.ex* results/current_run/l1x1_n2x2_i1_s0_p2
