#!/usr/bin/python
# lamboot                        First before running an MPI session
# mpirun -np 4 python Main.py    To run MPI with 4 processors


import numpy as np
import RayTrace as rt

# initialize the environment
Inv = rt.RayTrace([0,0,0], [1,0,0], 0.0, 80, 80)

# add objects
Inv.addSphere([10,0,-5], 5)
Inv.addSphere([12,10,-5], 5)
Inv.addTriangle([-20, -20, -10], [-20,  20, -10], [ 20, -20, -10])
Inv.addTriangle([ 20,  20, -10], [-20,  20, -10], [ 20, -20, -10])
Inv.addTriangle([ 20,  20, -10], [ 20, -20, -10], [ 20,  20,  20])
Inv.addTriangle([ 20,  20,  20], [ 20, -20, -10], [ 20, -20,  20])
Inv.addTriangle([ 20,  20,  20], [ 20,  20, -10], [-20,  20, -10])
Inv.addTriangle([ 20, -20,  20], [ 20, -20, -10], [-20, -20, -10])
Inv.addTriangle([ 20, -15, -10], [ 20, -20, -5 ], [ 15, -20, -10])
#Inv.addCone([15,0,-8], [15,0,16], 2.0, 8.0)

# add light sources
Inv.addLightSource([-5,80,100])
Inv.addLightSource([-10,20,0])
#Inv.addLightSource([16,0,0])

# show results
Inv.process()
Inv.show()
#Inv.saveHDF5("The.h5", 7)
#Inv.saveImage("output.png")
Inv.finalize()
