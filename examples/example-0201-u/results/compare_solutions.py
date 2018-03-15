#!/bin/python
#
# script to compare nodal field values between OpenCMISS-iron and CHeart
################################################################################
import numpy as np
import sys
import os
################################################################################
print("Numpy version: "+np.version.version)
################################################################################
NumberOfTests       = 0
NumberOfFailedTests = 0
failedtests_file = open("failed.tests", "w")
################################################################################
# 3D
################################################################################
print("  Testing 3D")
################################################################################
# exnode header size
skiprows      = 20
# size of block containing node ID, coordinates, solution values in exnode file
blocksize     = 15
# tolerance used to check numerical solution
tolu          = 1.0e-11

# stretch for analytic solution
lambd         = 1.2
mus           = 35.7
presval       = - (mus / 3.0 / lambd - mus * lambd * lambd / 3.0)

# number of elements in each coordinate direction for different refinement levels
NumberOfElements = np.zeros((3,3), dtype=int)
NumberOfElements[0][0] = 2
NumberOfElements[0][1] = 1
NumberOfElements[0][2] = 1
NumberOfElements[1][0] = 4
NumberOfElements[1][1] = 2
NumberOfElements[1][2] = 2
NumberOfElements[2][0] = 8
NumberOfElements[2][1] = 4
NumberOfElements[2][2] = 4

# width, height, length, NumberGlobalXElements, NumberGlobalYElements, NumberGlobalZElements, SolverIsDirect, JACOBIAN_FD ,MooneyRivlin1, MooneyRivlin2, useGeneratedMesh, BCDISP_MAX, bcType
# for all mesh types
for m in np.arange(0, 2, 1):
    # for all solver types
    for s in np.arange(0, 2, 1):
        # for all Jacobian types
        for j in np.arange(0, 2, 1):
            # for all spatial resolutions
            for r in np.arange(0, NumberOfElements.shape[0], 1):
                # for all BC types
                for bc in np.arange(0, 2, 1):
                    NumberOfTests  += 1
                    # current spatial resolution and number of nodes
                    nx  = NumberOfElements[r][0]
                    ny  = NumberOfElements[r][1]
                    nz  = NumberOfElements[r][2]
                    nn  = (2 * nx + 1) * (2 * ny + 1) * (2 * nz + 1)
                    nnP = (nx + 1) * (ny + 1) * (nz + 1)
                    # example folder: l2x1x1_n2x1x1_i2_s1_fd1_gm1_bc0
                    ####################################################################
                    # read reference data
                    foldername  = "current_run/l2x1x1_n"+str(nx)+"x"+str(ny)+"x"+str(nz)+"_i2_s"+str(s)+"_fd"+str(j)+"_gm"+str(m)+"_bc"+str(bc)+"/"
                    filename    = foldername + "Example.part0.exnode"
                    if not(os.path.isdir(foldername)):
                        status = "    Folder "+foldername+" does not exist."
                        print(status)
                        if (NumberOfFailedTests == 0):
                            failedtests_file.write("Failed tests:\n")
                        failedtests_file.write(status+"\n")
                        NumberOfFailedTests += 1
                        continue
                    linecount = 0
                    offset    = 0
                    errval    = 0.0
                    # load data
                    with open(filename) as f:
                        for line in f:
                            # skip header
                            if linecount < skiprows:
                                linecount += 1
                                continue
                            if ((linecount > skiprows+blocksize*nnP) and (linecount < skiprows*2+blocksize*nnP)):
                                linecount += 1
                                continue
                            if (linecount == skiprows*2+blocksize*nnP):
                                offset = skiprows - 1
                            linecount += 1
                            # node ID
                            if not(np.mod(linecount-skiprows-offset-1, blocksize)):
                                continue
                            # x-coordinate
                            elif not(np.mod(linecount-skiprows-offset-2, blocksize)):
                                ix      = float(line)
                            # y-coordinate
                            elif not(np.mod(linecount-skiprows-offset-3, blocksize)):
                                iy      = float(line)
                            # z-coordinate
                            elif not(np.mod(linecount-skiprows-offset-4, blocksize)):
                                iz      = float(line)
                            # x-displacement
                            elif not(np.mod(linecount-skiprows-offset-8, blocksize)):
                                # here we compare against the analytic solution
                                errval  += ((float(line) - lambd*ix) / 2.0)**2.0 / (3.0 * nn + nnP)
                            # y-displacement
                            elif not(np.mod(linecount-skiprows-offset-9, blocksize)):
                                # here we compare against the analytic solution
                                currentCoord_ocm    = float(line)
                                dispy_ocm           = currentCoord_ocm - iy
                                referenceCoordY_ocm = iy
                                referenceCoordY_ana = iy - 0.5
                                currentCoord_ana    = referenceCoordY_ana / np.sqrt(lambd)
                                dispy_ana           = currentCoord_ana - referenceCoordY_ana
                                errval             += (dispy_ocm - dispy_ana)**2.0 / (3.0 * nn + nnP)
                            # z-displacement
                            elif not(np.mod(linecount-skiprows-offset-10, blocksize)):
                                # here we compare against the analytic solution
                                currentCoord_ocm    = float(line)
                                dispz_ocm           = currentCoord_ocm - iz
                                referenceCoordZ_ocm = iz
                                referenceCoordZ_ana = iz - 0.5
                                currentCoord_ana    = referenceCoordZ_ana / np.sqrt(lambd)
                                dispz_ana           = currentCoord_ana - referenceCoordZ_ana
                                errval             += (dispz_ocm - dispz_ana)**2.0 / (3.0 * nn + nnP)
                            # pressure
                            elif not(np.mod(linecount-skiprows-offset-11, blocksize)):
                                # here we compare against the analytic solution
                                if (linecount > skiprows*2+blocksize*nnP):
                                    # this is a quadratic node, continue
                                    continue
                                else:
                                    errval  += ((float(line) - presval) / presval)**2.0 / (3.0 * nn + nnP)
                            else:
                                continue
                    errval  = np.sqrt(errval)
                    # check tolerance
                    if (errval > tolu):
                        status = filename+"       | analyticSolution - Iron |_2 = "+str(errval)
                        print(status)
                        failedtests_file.write(status+"\n")
                        NumberOfFailedTests += 1
if (NumberOfFailedTests == 0):
    failedtests_file.write("No failed tests.\n")
failedtests_file.close()
f       = open("results.summary", "w")
status  = "Passed tests: "+str(NumberOfTests-NumberOfFailedTests)+" / "+str(NumberOfTests)
print(status)
f.write(status+"\n")
f.close()
