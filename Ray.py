import numpy as np
import math

class Ray:

# + Position      The initial position the ray has
# + Direction     The direction the ray is pointing
# + Energy        The energy the ray initially has
#-----------------------
# + Bounce        Bounce the photon given the normal to the surface of contact. Also determine its energy
# + Scatter       Get the scattered ray toward the light source. Also determine its energy

   ReflEn = 0.6

   def __init__(self, Position, Direction, Energy):
      self.Position = Position # array of the position
      self.Direction = Direction # array as the direction from the position
      self.Energy = Energy # energy of the ray
   
##########################################################
   
   # set to the unscattered ray after reflection
   def Bounce(self, Point, Normal):
      r = self.Direction
      
      # use the Point of contact and the surface normal to determine the new scattered photon
      Pu_n = np.dot(np.outer(Normal,Normal),r) # projection of photon onto normal
      newDir = r - 2*Pu_n # new direction of unscattered photon
      
      self.Direction = newDir
      self.Position = Point
      self.Energy *= self.ReflEn
   
###########################################################
   
   # set to the scattered ray toward the light source after reflection
   def Scatter(self, Point, Normal, LightSource, UnScatRay):
      newDir = LightSource - Point
      
      # get the difference in angle between the reflected unscattered photon
      # and the scattered photon towards the light source
      Theta = math.acos(sum(newDir*UnScatRay.Direction)/(np.linalg.norm(newDir)*np.linalg.norm(UnScatRay.Direction)))
      
      # determine the Energy for the scattered photon
      if sum(UnScatRay.Direction*Normal) >= 0.0 and Theta < math.pi/2.0:
         Energy = math.cos(Theta)
#        Energy = ((math.pi/2.0 - Theta)/(math.pi/2.0))**2
      else: Energy = 0
      Energy *= UnScatRay.Energy
      
      self.Position = Point
      self.Direction = newDir
      self.Energy = Energy
      
