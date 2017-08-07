#!/usr/bin/env python
#
# script to compare nodal field values between OpenCMISS-iron and Matlab

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
sys.path.append('../src/scripts/')    # add path of exnode_reader
import exnode_reader
import os

print "Numpy version: "+np.version.version
print ""

# tolerance used to check numerical solution
tol_matlab = 12345678    # matlab solution differs from iron solution
tol_iron = 0.05

# number of elements in each coordinate direction for different refinement levels
NumberOfElements = np.zeros((1,2), dtype=int)
NumberOfElements[0][0] = 24
NumberOfElements[0][1] = 24

NumberOfTests       = 0
NumberOfFailedTests = 0

failedtests_file = open("failed.tests", "w")

# for problem types
for [nx,i,s,p] in [[24,1,0,1], [24,1,1,1], [10,1,0,1], [24,1,0,2], [24,1,0,8], [2,1,0,2]]:
  
  
  #print "case ",[nx,i,s,p]
  ny = nx
  for t in [0.01, 0.1, 0.2, 1, 2, 3]:
    #print "t=",t
    
    ####################################################################
    # read reference data
    # read matlab reference data
    foldername  = "reference/matlab/l1x1_n"+str(nx)+"x"+str(ny)+"_i"+str(i)+"_s0/"
    matlab_filename = "vm_"+str(t)+".csv"
    
    status = "   matlab file: {}\n".format(foldername+matlab_filename)
    
    matlab_data = None
    if os.path.exists(foldername):
      matlab_data = np.loadtxt(foldername+matlab_filename)
    
    # read iron reference data
    foldername  = "reference/iron/l1x1_n"+str(nx)+"x"+str(ny)+"_i"+str(i)+"_s"+str(s)+"/"
    pde_time_step = 0.005
    filenumber = int(t / pde_time_step)
    
    iron_reference_filename = "Time_2_"+str(filenumber)+".part0.exnode"
    
    
    status += "   iron reference file: {}\n".format(foldername+iron_reference_filename)
    iron_reference_data = exnode_reader.parse_file(foldername+iron_reference_filename, [["Vm", 1]])   # extract field Vm, component 1
    
    # read iron data of current run
    foldername  = "current_run/l1x1_n"+str(nx)+"x"+str(ny)+"_i"+str(i)+"_s"+str(s)+"_p"+str(p)+"/"            
    pde_time_step = 0.005
    filenumber = int(t / pde_time_step)
    
    iron_filename = "Time_2_"+str(filenumber)+".part0.exnode"
    status += "   iron file: {}\n".format(foldername+iron_filename)
    iron_data = exnode_reader.parse_file(foldername+iron_filename, [["Vm", 1]])   # extract field Vm, component 1
    
    if iron_data is None:
      if False:   # do not warn, when the iron file is not present, silently skip test. This may occur in the 'fast' target    
        print "Warning! no current iron data available for:\n"+status
      continue
    
    # determine which iron file would produce the least l2 error to the current matlab file
    if False:
      files = exnode_reader.get_exnode_files(foldername)
      min_l2diff = -1
      min_filename = ""
      for filename in files:
        iron_data = exnode_reader.parse_file(foldername+filename, [["Vm", 1]])   # extract field Vm, component 1
        
        l2diff = np.linalg.norm(matlab_data-iron_data) / np.linalg.norm(matlab_data)
        print "file ",filename, ", l2 norm: ", l2diff
        
        if l2diff < min_l2diff or min_l2diff == -1:
          min_l2diff = l2diff
          min_filename = filename
          
      print "t="+str(t)+", matlab:",matlab_filename,", best fit for iron: ",min_filename,", l2: ",min_l2diff
    
    # plot values
    if False:
      X = range(nx+1)
      Y = range(ny+1)
      X, Y = np.meshgrid(X, Y)
      Z_matlab = X.copy()
      Z_iron = X.copy()
      Z_diff = X.copy()
      for x in range(nx+1):
        for y in range(ny+1):
          Z_matlab[x,y] = matlab_data[x*(ny+1)+y]
          Z_iron[x,y] = iron_data[x*(ny+1)+y]
          Z_diff[x,y] = Z_matlab[x,y] - Z_iron[x,y]
      
      print "max:",np.amax(Z_matlab)
      
      fig = plt.figure(figsize=(18,8))
      plt.title("t="+str(t))
      ax = fig.add_subplot(131, projection='3d')
      ax.plot_surface(X, Y, Z_matlab, cmap=cm.coolwarm, linewidth=1, antialiased=False, rstride=1, cstride=1)
      plt.title("matlab file "+str(matlab_filename)+" max:"+str(np.amax(Z_matlab)))
      
      ax = fig.add_subplot(132, projection='3d')
      ax.plot_surface(X, Y, Z_iron, cmap=cm.coolwarm, linewidth=1, antialiased=False, rstride=1, cstride=1)
      plt.title("iron file "+str(iron_filename)+" max:"+str(np.amax(Z_iron)))
      
      ax = fig.add_subplot(133, projection='3d')
      ax.plot_surface(X, Y, Z_diff, cmap=cm.coolwarm, linewidth=1, antialiased=False, rstride=1, cstride=1)
      plt.title("matlab - iron")
      
      plt.tight_layout()
      plt.show()
      
    ####################################################################
    # compute difference
    
    # matlab reference - iron
    if matlab_data is not None:
      NumberOfTests += 1
      l2diff_matlab_iron = np.linalg.norm(matlab_data-iron_data) / np.linalg.norm(matlab_data)
      
      status = foldername+iron_filename+"       | Matlab   - Iron |_2 = "+str(l2diff_matlab_iron)+"\n"
      if l2diff_matlab_iron > tol_matlab:
        print status
        if (NumberOfFailedTests == 0):
            failedtests_file.write("Failed tests:\n")
        failedtests_file.write(status)
        NumberOfFailedTests += 1
      else:
        print "SUCCESS: ",status
    
    # iron reference - iron
    if iron_reference_data is not None:
      NumberOfTests += 1
      
      l2diff_iron_iron = np.linalg.norm(iron_reference_data-iron_data) / np.linalg.norm(iron_data)
      
      status = foldername+iron_filename+"       | Iron   - Iron |_2 = "+str(l2diff_iron_iron)+"\n"
      if l2diff_iron_iron > tol_iron:
        print status
        if (NumberOfFailedTests == 0):
            failedtests_file.write("Failed tests:\n")
        failedtests_file.write(status)
        NumberOfFailedTests += 1
      else:
        print "SUCCESS:",status
    
if (NumberOfFailedTests == 0):
  failedtests_file.write("No failed tests.\n")
failedtests_file.close()
f       = open("results.summary", "w")
status  = "Passed tests: "+str(NumberOfTests-NumberOfFailedTests)+" / "+str(NumberOfTests)+"\n"
print status
f.write(status)
f.close()
