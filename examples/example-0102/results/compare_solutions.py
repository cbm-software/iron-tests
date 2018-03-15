#!/bin/python
import subprocess

# command to execute
command = "matlab -nodesktop -nosplash -r \"compare_solutions\""
print command

# call Matlab
subprocess.call(command, shell=True)

