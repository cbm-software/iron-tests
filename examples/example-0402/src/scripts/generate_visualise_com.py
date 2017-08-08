#!/usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess

# This script adds the names of all *.exnode files in the current directory as cmgui load command to the beginning of visualise.com.
# It is intentionally not fully automated, e.g. you have to provide your cmgui executable in cmgui_command and check the output in visualise.com. It will overwrite that file, so multiple invocations of this script generate the input files more often.



import os
import sys
import numpy as np
import csv
import collections
import copy
import glob
import datetime
from sets import Set

show_file = "visualise.com"
cmgui_command = "/store/software/opencmiss/cmgui_2.7/cmgui-wx > /dev/null 2>&1"
 
"""
  compare function for sort
  compare OpenCMISS output files by their iteration number *_<no>_.part*
"""
def compare(x,y):
  if ".part" in x:
    posend = x.find(".part")
    pos = x.rfind("_", 0, posend)+1
    strx = x[pos:posend]

  if ".part" in y:
    posend = y.find(".part")
    pos = y.rfind("_", 0, posend)+1
    stry = y[pos:posend]
  
  try:  
    nx = int(strx)
  except:
    nx = -5
    
  try:
    ny = int(stry)
  except:
    ny = -5
    
  #print "x={}, y={}, nx={}, ny={}".format(x,y,nx,ny)
    
  if nx == ny:
    if "_M_" in x and not "_M_" in y:
      return 1
    elif "_M_" in y and not "_M_" in x:
      return -1
    elif "_FE" in x and not "_FE" in y:
      return -1
    elif "_FE" in y and not "_FE" in x:
      return 1
    
  return nx - ny
      

working_directory = subprocess.check_output("pwd", shell=True)
#print "pwd: [{}]".format(working_directory)

# get all input data in current directory, sorted by number of file name
ls = glob.glob('Time_2_*.exnode')
directory_content = sorted(ls, cmp=compare)

#print directory_content

if False:
  # replace NAN values in files by 0
  print "converting NANs in files . . ."
  i = 0
  n_nans_total = 0
  for filename in directory_content:
    contents = []
    
    sys.stdout.write("\b\b\b\b\b\b\b\b\b{}%".format(int(float(i)/len(directory_content)*100.0)))
    sys.stdout.flush()
    
    n_nans = 0
    with open(temporary_directory+"/"+filename, "rb") as file_in:
      for line in file_in:
        if "NAN" in line:
          contents += "  0.0000000000000000E+00\n"
          n_nans += 1
          n_nans_total += 1
        else:
          contents += line
        
    if n_nans > 0:
      print "   {} NANs in {}".format(n_nans, filename)
      with open(temporary_directory+"/"+filename, "wb") as file_out:
        for line in contents:
          file_out.write(line)
    i += 1
  if n_nans_total==0:
    print "\b\b\b\bno NANs"
  else:
    print "\b\b\b\b{} NANs total".format(n_nans_total)
  print "\b\b\b\bdone"

# read in show file
show_file_contents = []
with open(show_file, "rb") as f:
  for line in f:
    show_file_contents += line

# write out show file
with open(show_file, "w") as fout:

  fout.write("# automatically generated at {}\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
  fout.write("# nodes\n")

  # iterate over output files
  for filename in directory_content:
    print "  filename: {}".format(filename)
    if ".exnode" in filename and "Time" in filename and ".part" in filename:
      try:
        #print "pos: {} [{}]".format(filename.rfind("_")+1, filename[filename.rfind("_")+1:])
        time = int(filename[filename.rfind("_")+1:filename.find(".part")])
      except:
        continue
        
      #print "time: {}".format(time)  
      
      if "_M_" in filename:
        fout.write("gfx read node "+filename+" time "+str(time)+" node_offset 100000\n")
      else:
        fout.write("gfx read node "+filename+" time "+str(time)+"\n")
        
  fout.write("\n# elements\n")
        
  # write *.exelem files
  for filename in directory_content:
    if ".exelem" in filename:
      if "_M" in filename:
        fout.write("gfx read elem {} generate element_offset 100000 node_offset 100000 line_offset 100000 face_offset 100000\n".format(filename))
      else:
        fout.write("gfx read elem {} generate\n".format(filename))
        
  fout.write("\n# ----------------------\n")
        
  # write rest of file
  for line in show_file_contents:
    fout.write(line)
    
print "Found "+str(len(directory_content))+" files, added to "+show_file+"."
    
# call cmgui

#cmd = "cd {tmp} && {cmgui} show.com".format(cmgui=cmgui_command, tmp=temporary_directory)
#print cmd
#subprocess.check_output(cmd, shell=True)
