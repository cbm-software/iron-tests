#!/bin/bash

echo "compiling and running example $(pwd)"

echo "  running 3D tests"
# arguments:
# width, height, length, NumberGlobalXElements, NumberGlobalYElements, NumberGlobalZElements, SolverIsDirect, JACOBIAN_FD ,MooneyRivlin1, MooneyRivlin2, useGeneratedMesh, BCDISP_MAX, bcType
#
# generated mesh - direct solver
./bin/example 2.0 1.0 1.0  4 2 2 1 1 17.85 0.0 1 0.2 0
./bin/example 2.0 1.0 1.0  4 2 2 1 1 17.85 0.0 1 0.2 1
./bin/example 2.0 1.0 1.0  8 4 4 1 1 17.85 0.0 1 0.2 0
./bin/example 2.0 1.0 1.0  8 4 4 1 1 17.85 0.0 1 0.2 1
./bin/example 2.0 1.0 1.0 16 8 8 1 1 17.85 0.0 1 0.2 0
./bin/example 2.0 1.0 1.0 16 8 8 1 1 17.85 0.0 1 0.2 1
# generated mesh - iterative solver
./bin/example 2.0 1.0 1.0  4 2 2 0 1 17.85 0.0 1 0.2 0
./bin/example 2.0 1.0 1.0  4 2 2 0 1 17.85 0.0 1 0.2 1
./bin/example 2.0 1.0 1.0  8 4 4 0 1 17.85 0.0 1 0.2 0
./bin/example 2.0 1.0 1.0  8 4 4 0 1 17.85 0.0 1 0.2 1
./bin/example 2.0 1.0 1.0 16 8 8 0 1 17.85 0.0 1 0.2 0
./bin/example 2.0 1.0 1.0 16 8 8 0 1 17.85 0.0 1 0.2 1

# user-defined mesh - direct solver
./bin/example 2.0 1.0 1.0  4 2 2 1 1 17.85 0.0 0 0.2 0
./bin/example 2.0 1.0 1.0  4 2 2 1 1 17.85 0.0 0 0.2 1
./bin/example 2.0 1.0 1.0  8 4 4 1 1 17.85 0.0 0 0.2 0
./bin/example 2.0 1.0 1.0  8 4 4 1 1 17.85 0.0 0 0.2 1
./bin/example 2.0 1.0 1.0 16 8 8 1 1 17.85 0.0 0 0.2 0
./bin/example 2.0 1.0 1.0 16 8 8 1 1 17.85 0.0 0 0.2 1
# user-defined mesh - iterative solver
./bin/example 2.0 1.0 1.0  4 2 2 0 1 17.85 0.0 0 0.2 0
./bin/example 2.0 1.0 1.0  4 2 2 0 1 17.85 0.0 0 0.2 1
./bin/example 2.0 1.0 1.0  8 4 4 0 1 17.85 0.0 0 0.2 0
./bin/example 2.0 1.0 1.0  8 4 4 0 1 17.85 0.0 0 0.2 1
./bin/example 2.0 1.0 1.0 16 8 8 0 1 17.85 0.0 0 0.2 0
./bin/example 2.0 1.0 1.0 16 8 8 0 1 17.85 0.0 0 0.2 1
