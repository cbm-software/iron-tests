#!/bin/bash

# run 2D
cheartsolver.out l2x1x0_n2x1x0_i1.P
cheartsolver.out l2x1x0_n4x2x0_i1.P
cheartsolver.out l2x1x0_n8x4x0_i1.P
cheartsolver.out l2x1x0_n2x1x0_i2.P
cheartsolver.out l2x1x0_n4x2x0_i2.P
cheartsolver.out l2x1x0_n8x4x0_i2.P
# run 3D
cheartsolver.out l2x1x1_n2x1x1_i1.P
cheartsolver.out l2x1x1_n4x2x2_i1.P
cheartsolver.out l2x1x1_n8x4x4_i1.P
cheartsolver.out l2x1x1_n2x1x1_i2.P
cheartsolver.out l2x1x1_n4x2x2_i2.P
cheartsolver.out l2x1x1_n8x4x4_i2.P

rm -f gmon.out
