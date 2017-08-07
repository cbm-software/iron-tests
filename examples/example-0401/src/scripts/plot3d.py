#!/usr/bin/env python
#
# script to plot Vm from all the exnode files as 2D surface for all timesteps in an animation

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import exnode_reader
from matplotlib import animation
import sys
import os

def animate(t, create_image=False):
  global foldername, image_name
  ax1.clear()
  
  # get all exnode files sorted in that directory
  files = exnode_reader.get_exnode_files(foldername)
  
  if t >= len(files):
    return
  iron_filename = files[t]
  
  iron_data = exnode_reader.parse_file(foldername+iron_filename, [["Vm", 1]])   # extract field Vm, component 1
  
  n = len(iron_data)
  nx = int(np.sqrt(n))-1
  ny = int(np.sqrt(n))-1
  print "t=",t,"/",len(files),", n=",n,", nx=",nx, ", ny=", ny

  X = range(nx+1)
  Y = range(ny+1)
  X, Y = np.meshgrid(X, Y)
  Z_iron = X.copy()
  for x in range(nx+1):
    for y in range(ny+1):
      try:
        Z_iron[x,y] = iron_data[x*(ny+1)+y]
      except:    # fails for nans
        Z_iron[x,y] = 0.0

  ax1.plot_surface(X, Y, Z_iron, cmap=cm.coolwarm, linewidth=1, antialiased=False, rstride=1, cstride=1)
  ax1.set_zlim(-80,30)
  plt.title("t = "+str(t*0.005))
  
  if t % 100 == 0 or create_image:
    filenamepng = image_name+"_t"+str(t)+".png"
    filenameeps = image_name+"_t"+str(t)+".eps"
    plt.savefig(filenamepng)
    print "Image saved to {}".format(filenamepng)
    os.system("convert "+filenamepng+" "+filenameeps)
      
# get folder name from command line argument
foldername = "."
if len(sys.argv) > 1:
  foldername = sys.argv[1]

# append slash if necessary
if foldername[-1] != "/":
  foldername = foldername + "/"

# get kind (fast or big)
kind = 'fast'
if len(sys.argv) > 2:
  kind = sys.argv[2]

# create file names
image_name = "{}".format(foldername).replace("/", "_")
animation_name = "{}.mp4".format(foldername).replace("/", "_")

# create figure
fig = plt.figure(1)
ax1 = fig.gca(projection='3d')
plt.tight_layout()

# create image of first time step
animate(0, create_image=True)

# create image of last time step
files = exnode_reader.get_exnode_files(foldername)
animate(len(files)-1, create_image=True)

# save to animation
if kind == 'fast':
  print "Do not save mp4 animation because of 'fast' target."
else:
  # call animator to create mp4 animation
  anim = animation.FuncAnimation(fig, animate,
                                 frames=1000, interval=20)
  anim.save(animation_name)
  print "Animation saved to {}".format(animation_name)

#plt.show()


