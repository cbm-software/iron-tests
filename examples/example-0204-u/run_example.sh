#!/bin/bash

echo "running example $(pwd)"

# arguments: SolverIsDirect, JACOBIAN_FD, MooneyRivlin1, MooneyRivlin2, BCDISP_MAX, NumberOfLoadIncrements, bcType

# direct solver
./bin/example ./src/cheart/meshes/FinalModel_quad_FE 1 1 17.85 0.0 2.0 1 0
./bin/example ./src/cheart/meshes/FinalModel_quad_FE 1 0 17.85 0.0 2.0 1 0

# iterative solver
./bin/example ./src/cheart/meshes/FinalModel_quad_FE 0 1 17.85 0.0 2.0 1 0
./bin/example ./src/cheart/meshes/FinalModel_quad_FE 0 0 17.85 0.0 2.0 1 0
