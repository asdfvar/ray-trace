# added MPI parallelization. use mpirun -np (number of procs) python (program)

import numpy as np
import pylab as pl
import math
import Shape
import Ray
import copy
import pypar
import time

class RayTrace:

# + Objects
# + LightSources
# + WindowDistance
# + WindowHeight
# + WindowWidth
# + WindowRows
# + WindowCols
# + LookPos
# + LookDir
# + Yaw
# - Lat
# - Lon
# + Window
# + Visual
#----------------
# + addLightSource (self, Position)
# + addSphere (self, Position, Radius)
# + addTriangle (self, q1, q2, q3)
# + addCone(self, p1, p2, r1, r2)
# - RotVec (self, Vec, Lon, Lat, Yaw)
# - getEn (self, Energy, LightSource)
# + show (self)
# + saveImage (self, name)
# + saveHDF5 (self, name, fileNum)
# + clear (self)
# + finalize(self)

   Objects = []
   LightSources = []
   WindowDistance = .5
   WindowHeight = 2
   WindowWidth = 1.5
   WindowRows = 400 # must be divisable by the number of processes
   WindowCols = 400
   
   # must define the environment variables
   def __init__(self, LookPos, LookDir, LookYaw, WindowRows = 40, WindowCols = 40):
      self.LookPos = np.array(LookPos)
      self.LookDir = np.array(LookDir)
      self.Yaw = LookYaw
      self.WindowRows = WindowRows
      self.WindowCols = WindowCols
      rhop = np.linalg.norm(np.array([LookDir[0],LookDir[1]]))
      self.__Lon = math.atan2(LookDir[1], LookDir[0])
      self.__Lat = math.atan2(LookDir[2],rhop)
      self.start = time.time()
      
      # initialize the MPI
      self.numproc = pypar.size()
      self.myid =    pypar.rank()
      self.node =    pypar.get_processor_name()
      
      if self.myid != self.numproc - 1:
         self.Rows = self.WindowRows/self.numproc
         self.RowEnd = self.WindowRows/self.numproc * (self.myid+1) - 1
      else:
         self.Rows = self.WindowRows/self.numproc + self.WindowRows%self.numproc
         self.RowEnd = self.WindowRows

      self.RowStart = self.WindowRows/self.numproc * self.myid
      self.Window = np.zeros(shape = (self.Rows, self.WindowCols))


   # add a light source to the list of light sources
   def addLightSource(self, Position):
      self.LightSources.append(np.array(Position))
         
   # add sphere to the list of objects
   def addSphere(self, Position, Radius):
      Sphere = Shape.Sphere(Position, Radius)
      self.Objects.append(Sphere)
   
   # add a triangle to the list of objects
   def addTriangle(self, q1, q2, q3):
      Triangle = Shape.Triangle(q1, q2, q3)
      self.Objects.append(Triangle)
   
   # add a cone to the list of objects
   def addCone(self, p1, p2, r1, r2):
      Cone = Shape.Cone(p1, p2, r1, r2)
      self.Objects.append(Cone)
   
   def __RotVec(self, Vec, Lon, Lat, Yaw):
      Rotx = np.array([[1,0,0],[0,math.cos(Yaw),-math.sin(Yaw)],[0,math.sin(Yaw),math.cos(Yaw)]])
      Roty = np.array([[math.cos(Lat),0,math.sin(Lat)],[0,1,0],[-math.sin(Lat),0,math.cos(Lat)]])
      Rotz = np.array([[math.cos(Lon),-math.sin(Lon),0],[math.sin(Lon),math.cos(Lon),0],[0,0,1]])

      Vec = np.dot(Rotx,Vec)
      Vec = np.dot(Roty,Vec)
      Vec = np.dot(Rotz,Vec)

      return Vec
   
#########################################################################

   # check if the photon intersects any of the objects. if so, which one first   
   def __IntersectObj(self, Photon):
      MinDist = float("inf"); MinInd = float("nan")
      Intersects = False
      for ObjNum in range(len(self.Objects)):
         if self.Objects[ObjNum].Intersects(Photon):
            Distance = self.Objects[ObjNum].Dist(Photon)
            if Distance < MinDist:
               Intersects = True
               MinDist = Distance
               MinInd = int(ObjNum)
      return [Intersects, MinInd, MinDist]

#########################################################################

   # get the energy of the scattered ray after the Ray has hit a surface
   def __getEn(self, Ray, Light):
      # find out which object the photon hits first if any
      Obj = self.__IntersectObj(Ray)
      Intersects = Obj[0]; MinInd = Obj[1]
      
      # if the photon intersects an object, get the scattered and unscattered photons
      if Intersects:
         Cont = True
         Ret = self.Objects[MinInd].Reflection(Ray, Light) # [un-scattered ray, scattered ray]
         UnScatRay = Ret[0]
         ScatRay = Ret[1]
              
         # check if the scattered photon intersects any objects
         Obj = self.__IntersectObj(ScatRay)
         Intersects = Obj[0]
         DistObj = Obj[2]
         
         # get the distance to the light source from the scattered photon
         DistLight = np.linalg.norm(Light - ScatRay.Position)
         
         # if the photon hits an object before the light source, it is a shadow
         if Intersects and DistLight >= DistObj: ScatRay.Energy = 0.0
         
         # get the energy
         Energy = ScatRay.Energy
         
         return [Energy, UnScatRay, Cont]
      else:
         return [0, Ray, False]
      
      
#########################################################################
   
   # process the image with all the objects in the scene
   def process(self):
      # loop through all the rays that go through the window
      for row in xrange(self.Rows):
         for col in xrange(self.WindowCols):

            # loop through all the light sources
            for Light in self.LightSources:
               Energy = 0.0
                  
               # move the ray to the corresponding pixel on the window. dx and dy changes
               MoveUp = self.WindowHeight/2.0 - (row+self.RowStart - 0.5)*\
                        self.WindowHeight/self.WindowRows # starts from top most pixel
               MoveRight = -self.WindowWidth/2.0 + (col + 0.5)*\
                           self.WindowWidth/self.WindowCols # starts from left most pixel
                  
               # define the ray and rotate (in look direction)
               # then translate (window distance)
               Vec = np.array([self.WindowDistance, -MoveRight, MoveUp])
               Vec = self.__RotVec(Vec, self.__Lon, self.__Lat, self.Yaw)
               Photon = Ray.Ray(self.LookPos,Vec,1.0)
               
               Cont = True
               
		         # get energy after first contact
               Res = self.__getEn(Photon, Light)
               Energy += Res[0]
               UnScatRay = Res[1]
               Cont = Res[2]
               
               #pypar.barrier();
               
               for k in range(1):
                  # get energy after N contacts
                  if Cont:
                     Res = self.__getEn(UnScatRay, Light)
                     Energy += Res[0]
                     UnScatRay = Res[1]
                     Cont = Res[2]
                  
               self.Window[row, col] += Energy # set the energy to the corresponding window
  
      # stitch together the images from each process
      if self.myid != 0:
         pypar.send(self.Window, 0)
      else:
         self.Visual = np.zeros(shape = (self.WindowRows, self.WindowCols))
         RowEnd = self.WindowRows/self.numproc
         self.Visual[:RowEnd, :] = self.Window.copy()
         for Pr in range(1,self.numproc):            
            RowStart = self.WindowRows/self.numproc * Pr
            RowEnd = self.WindowRows/self.numproc * (Pr+1)
            self.Visual[RowStart:RowEnd,:] = pypar.receive(Pr)
      
         print "Elapsed time for process %d is %f" % (self.myid, time.time() - self.start)
      
#########################################################################
      
   # show the image with all the objects in the scene
   def show(self):
      if self.myid == 0:
         pl.gray()
         pl.imshow(self.Visual, interpolation='nearest')
         #pl.pcolor(self.Window)
         pl.show()
   
#########################################################################
   
   def saveImage(self, name):
      import Image
      if self.myid == 0:
         #Rescale to 0-255 and convert to uint8
         rescaled = (255.0 / self.Visual.max() * (self.Visual - self.Visual.min())).astype(np.uint8)
         im = Image.fromarray(rescaled)
         im.save(name)
         
#########################################################################

   def saveHDF5(self, name, fileNum):
      import h5py
      
      if self.myid == 0:
         # create file if it doesn't exist
         fid = h5py.File("Images.h5",'a')

         # create group
         grp = fid.create_group("frames")

         # create dataset
         dset = grp.create_dataset("Picture" + str(fileNum), (self.WindowRows, self.WindowCols), '=f8', maxshape=(None, None))
         dset.write_direct(self.Visual)
         
         # close hdf5
         fid.close()

#########################################################################

   # clear all shape objects and light sources
   def clear(self):
      self.Objects = []
      self.LightSources = []
      
#########################################################################

   def finalize(self):
      pypar.finalize()
      
