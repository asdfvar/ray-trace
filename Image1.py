#!/usr/bin/python
# lamboot                        First before running an MPI session
# mpirun -np 4 python Main.py    To run MPI with 4 processors


import numpy as np
import RayTrace as rt

# initialize the environment
Inv = rt.RayTrace([0,0,0], [1,0,0], 0.0, 200, 200)

# add objects
Inv.addSphere([9,0,0], 5)

Inv.addTriangle([ 5, 0, -5], [ 10, 7, -5 ], [ 10, -7, -5])
Inv.addTriangle([ 5, 0, -5], [ 10, 7, -5 ], [ 9, 0, -30])
Inv.addTriangle([ 5, 0, -5], [ 10, -7, -5 ], [ 9, 0, -30])

# add light sources
Inv.addLightSource([-5,80,100])
Inv.addLightSource([-10,20,0])
Inv.addLightSource([-10,-20,0])

# show results
Inv.process()
Inv.show()
#Inv.saveHDF5("The.h5", 7)
Inv.saveImage("output.png")
Inv.finalize()
