#!/bin/python
#
# script to compare nodal field values between OpenCMISS-iron and CHeart
################################################################################
import numpy as np
import sys
################################################################################
print "Numpy version: "+np.version.version
################################################################################
NumberOfTests       = 0
NumberOfFailedTests = 0
################################################################################
# 3D
################################################################################
print "  Testing 3D"
################################################################################
# exnode header size
skiprows    = 20
# size of block containing node ID, coordinates, solution values in exnode file
blocksize   = 15
# tolerance used to match nodes between CHeart and iron
tolx         = 1.0e-6
# tolerance used to check numerical solution
tolu         = 1.0e-10

failedtests_file = open("failed.tests", "w")

# for all solver types
for s in np.arange(1, 2, 1):
    # for Jacobian type
    for i in np.arange(1, 2, 1):
        # for BC type
        for r in np.arange(0, 1, 1):
            NumberOfTests += 1
            ####################################################################
            # read reference data
            foldername  = "reference/cheart/"
            filename    = foldername + "SolidSpace-1.D"
            # load CHeart space variable, skip first line
            cxyz        = np.loadtxt(filename, dtype=float, skiprows=1)
            filename    = foldername + "Disp-1.D"
            filenameP   = foldername + "SolidPres-1.D"
            # load CHeart scalar variable, skip first line
            cvals       = np.loadtxt(filename,  dtype=float, skiprows=1)
            cvalsP      = np.loadtxt(filenameP, dtype=float, skiprows=1)
            # now start reading exnode data
            ivals0      = 0.0 * cvals
            ivals0P     = 0.0 * cvalsP
            foldername  = "reference/iron/"+"s"+str(s)+"_fd"+str(i)+"_bc"+str(r)+"/"
            filename    = foldername + "Example.part0.exnode"
            linecount = 0
            matched   = 0
            offset    = 0
            with open(filename) as f:
                for line in f:
                    # skip header
                    if linecount < skiprows:
                        linecount += 1
                        continue
                    if ((linecount > skiprows+blocksize*ivals0P.shape[0]) and (linecount < skiprows*2+blocksize*ivals0P.shape[0])):
                        linecount += 1
                        continue
                    if (linecount == skiprows*2+blocksize*ivals0P.shape[0]):
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
                    # scalar value
                    elif not(np.mod(linecount-skiprows-offset-8, blocksize)):
                        ival    = float(line)
                        # match nodes between CHeart and iron
                        # this can be expensive, but we only consider small spatial resolution here...
                        for node in np.arange(0, cxyz.shape[0], 1):
                            cx      = cxyz[node][0]
                            cy      = cxyz[node][1]
                            cz      = cxyz[node][2]
                            l2diff  = np.sqrt((cx-ix)**2.0+(cy-iy)**2.0+(cz-iz)**2.0)
                            if (l2diff < tolx):
                                ivals0[node] = ival - ix
                                matched     += 1
                                break
                    else:
                        continue
            if (matched != cxyz.shape[0]):
                print "wrong total number of nodes: ",matched," ",cxyz.shape[0]
                sys.exit()
            ####################################################################
            # read current_run
            ivals       = 0.0 * cvals
            foldername  = "current_run/"+"s"+str(s)+"_fd"+str(i)+"_bc"+str(r)+"/"
            filename    = foldername + "Example.part0.exnode"
            linecount = 0
            matched   = 0
            offset    = 0
            with open(filename) as f:
                for line in f:
                    # skip header
                    if linecount < skiprows:
                        linecount += 1
                        continue
                    if ((linecount > skiprows+blocksize*ivals0P.shape[0]) and (linecount < skiprows*2+blocksize*ivals0P.shape[0])):
                        linecount += 1
                        continue
                    if (linecount == skiprows*2+blocksize*ivals0P.shape[0]):
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
                    # scalar value
                    elif not(np.mod(linecount-skiprows-offset-8, blocksize)):
                        ival    = float(line)
                        # match nodes between CHeart and iron
                        # this can be expensive, but we only consider small spatial resolution here...
                        for node in np.arange(0, cxyz.shape[0], 1):
                            cx      = cxyz[node][0]
                            cy      = cxyz[node][1]
                            cz      = cxyz[node][2]
                            l2diff  = np.sqrt((cx-ix)**2.0+(cy-iy)**2.0+(cz-iz)**2.0)
                            if (l2diff < tolx):
                                ivals[node]  = ival - ix
                                matched     += 1
                                break
            if (matched != cxyz.shape[0]):
                print "wrong total number of nodes",matched," ",cxyz.shape[0]
                sys.exit()
            ####################################################################
            #l2diff_ci   = np.linalg.norm(cvals-ivals, 2)
            #l2diff_i0i  = np.linalg.norm(ivals0-ivals, 2)
            # compute difference; normalized RMSE
            l2diff_ci   = np.linalg.norm(cvals-ivals, 2) / np.sqrt(cvals.shape[0]) / 2.2
            l2diff_i0i  = np.linalg.norm(ivals0-ivals, 2) / np.sqrt(cvals.shape[0]) / 2.2
            if ((l2diff_ci > tolu) and (l2diff_i0i > tolu)):
                status = filename+"       | CHeart   - Iron |_2 = "+str(l2diff_ci) \
                    +"       | Iron_ref - Iron |_2 = "+str(l2diff_i0i)
                print status
                if (NumberOfFailedTests == 0):
                    failedtests_file.write("Failed tests:\n")
                failedtests_file.write(status+"\n")
                NumberOfFailedTests += 1
            elif (l2diff_ci > tolu):
                status = filename+"       | CHeart   - Iron |_2 = "+str(l2diff_ci)
                print status
                failedtests_file.write(status+"\n")
                NumberOfFailedTests += 1
            elif (l2diff_i0i > tolu):
                status = filename+"       | Iron_ref - Iron |_2 = "+str(l2diff_i0i)
                print status
                failedtests_file.write(status+"\n")
                NumberOfFailedTests += 1
if (NumberOfFailedTests == 0):
    failedtests_file.write("No failed tests.\n")
failedtests_file.close()
f       = open("results.summary", "w")
status  = "Passed tests: "+str(NumberOfTests-NumberOfFailedTests)+" / "+str(NumberOfTests)
print status
f.write(status+"\n")
f.close()
