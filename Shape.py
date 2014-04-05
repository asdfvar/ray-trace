import numpy as np
import Ray
import math
import copy
from numpy import linalg as la

##################################

class Shape: # parent class

# + Center
# + Radius
#----------------
# + Intersects(Ray.Ray)             determine if the ray intersects the object
# + Dist(Ray.Ray)                   determine the distance the ray travels until it intersects the object
# + Reflection(Ray.Ray, np.array)   get the scattered and unscattered rays
   
   Center = np.array([1,1,1])
   Radius = 1
   
####################################################################################
####################################################################################

class Sphere(Shape):

# + Center
# + Radius
#----------
# + Intersects(Ray.Ray)
# + Dist(Ray.Ray)
# + Reflection(Ray.Ray, numpy.array)
   
   def __init__(self, Center, Radius):
      print 'Sphere added'
      self.Center = np.array(Center)
      self.Radius = Radius

   ################

   # check if the Ray intersects the sphere
   def Intersects(self, Photon):
      if isinstance(Photon, Ray.Ray):
         # get vector from Photon position to sphere center (everything relative to photon position)
         r = Photon.Direction
         c = self.Center - Photon.Position
         
         # check to see if photon is even in the same direction as the sphere
         if sum(r*c) < 0.0: return False
         
         #project c onto r and get length
         Proj = sum(r*c)*r/np.linalg.norm(r)**2
         Length = Proj - c
         
         if np.linalg.norm(Length) <= self.Radius: return True
         else: return False
      else: print "Argument must be of type Ray"

   ################

   # get the distance to this object
   def Dist(self, Photon):
      if isinstance(Photon, Ray.Ray):
         if self.Intersects(Photon):
            # define the photon direction and sphere position (relative to the photon)
            r = Photon.Direction
            c = self.Center - Photon.Position
            
            # get the point of contact distance
            rNorm = r/np.linalg.norm(r)
            return sum(c*rNorm) - math.sqrt(sum(c*rNorm)**2 - sum(c*c) + self.Radius**2) # first point of contact
         else: return float("inf")
      else: print "Argument must be of type Ray"

   ###############

   # get the scattered and unscattered reflection of the photon (unscattered to the light source)
   def Reflection(self, Photon, LightSource):
      # find the point where the photon hits the sphere
      # determine where the photon will project from here without scatter
      # determine the direction the scattered photon will reach the light source
      
      if isinstance(Photon, Ray.Ray) and self.Intersects(Photon):
         # define the photon direction and sphere position (relative to the photon)
         r = Photon.Direction
         c = self.Center - Photon.Position
         
         # get the point of contact
         rNorm = r/np.linalg.norm(r)
         Length = sum(c*rNorm) - math.sqrt(sum(c*rNorm)**2 - sum(c*c) + self.Radius**2) # first point of contact
         Point = Photon.Position + Length*rNorm # the point of contact
         
         # get direction where the photon will project without scatter
         n = (Point - c)/np.linalg.norm(Point - c) # normal to surface

         UnScatRay = copy.copy(Photon)
         UnScatRay.Bounce(Point, n)
         
         ScatRay = copy.copy(Photon)
         ScatRay.Scatter(Point, n, LightSource, UnScatRay)
         
         return [UnScatRay, ScatRay] # unscattered ray, scattered ray (toward light source)
      elif not(isinstance(Photon, Ray.Ray)): print "Photon is not of the Ray type"; return -1
      elif not(self.Intersects(Photon)): print "Photon does not intersect sphere"; return -1
      else: return -1

#############################################################################
#############################################################################

class Triangle(Shape):
   # + Center
   # + Radius
   #----------
   # + Intersects(Ray.Ray)
   # + Dist(Ray.Ray)
   # + Reflection(Ray.Ray, numpy.array)
   
   def __init__(self, q1, q2, q3):
      print 'Triangle added'
      self.q1 = np.array(q1)
      self.q2 = np.array(q2)
      self.q3 = np.array(q3)
   
######################################################

   def Intersects(self, Photon):
      l = Photon.Position
      m = Photon.Direction
   
      # find vectors u and v relative to q1
      u = self.q2 - self.q1
      v = self.q3 - self.q1

      n = np.cross(u,v)
      n = n/np.linalg.norm(n)

      # get p from p = n1x + n2y + n3z
      p = sum(n*self.q1)

      # find the time (t)
      if sum(n*m)!=0:
         t = (p - sum(n*l))/sum(n*m)
         
         # check to see if photon is even in the same direction as the triangle
         if t < 1.0e-10: return False
         
         # get the intersection point
         Point = l + m*t
         
         # test if the intersection point falls within the 3 points.
         # check if the cross product direction of the triangle corners is the same
         # compared to the point in question relative to all the corners of the triangle (in the same order)
         Direction = np.cross(self.q2 - self.q1, self.q3 - self.q2) # get direction of cross product
         q12P = np.cross(self.q2 - self.q1, Point - self.q1)
         q23P = np.cross(self.q3 - self.q2, Point - self.q2)
         q31P = np.cross(self.q1 - self.q3, Point - self.q3)
         if sum(q12P*Direction) > 0 and sum(q23P*Direction) > 0 and sum(q31P*Direction) > 0:
            return True
         else:
            return False
      else:
         return False
         
###############################################

   def Reflection(self, Photon, LightSource):
      r = Photon.Direction
      
      u = self.q2 - self.q1
      v = self.q3 - self.q1
      
      # get the direction of the unscattered photon
      n = np.cross(u,v)
      if sum(r*n) > 0: n = -n # make sure the normal is facing the photon
      
      # get the point of contact
      rNorm = r/np.linalg.norm(r)
      Length = self.Dist(Photon)
      Point = Photon.Position + Length*rNorm # the point of contact

      UnScatRay = copy.copy(Photon)
      UnScatRay.Bounce(Point, n)
         
      ScatRay = copy.copy(Photon)
      ScatRay.Scatter(Point, n, LightSource, UnScatRay)

      return [UnScatRay, ScatRay] # unscattered ray, scattered ray (toward light source)

####################################

   def Dist(self, Photon):
      l = Photon.Position
      m = Photon.Direction
      
      # find vectors u and v relative to q1
      u = self.q2 - self.q1
      v = self.q3 - self.q1

      n = np.cross(u,v)
      n = n/np.linalg.norm(n)

      # get p from p = n1x + n2y + n3z
      p = sum(n*self.q1)

      # find the time (t)
      if sum(n*m)!=0:
         t = (p - sum(n*l))/sum(n*m)
         
         # get the intersection point
         self.Point = l + m*t
         
         # get the distance
         return np.linalg.norm(self.Point - Photon.Position)
      else:
         return -1

####################################################################################
####################################################################################

class Cone(Shape):
# Finite cone (the general concept of the cylindar). This cuts the infinite (meaning the shape used in studies of conics) in length in all directions, cone.
# + p1
# + p2
# + r1
# + r2
# - length        finite section length of the infinite cone
# - vertex        the vertex of the cone
# - direction     direction from vertex to p1 and p2
# - angle         angle the cone makes from the direction of the cone
#----------
# + Intersects(Ray.Ray)
# + Dist(Ray.Ray)
# + Reflection(Ray.Ray, numpy.array)

   import math

   def __init__(self, p1, p2, r1, r2):
      print 'Cone added'
      from numpy import linalg as la
      
      p1 = np.array(p1)
      p2 = np.array(p2)
      self.p1 = p1
      self.p2 = p2
      self.r1 = r1
      self.r2 = r2
      
      # get the vertex, the unit direction from the vertex, and the angle the cone spans
      self.__length = la.norm(p2 - p1)
      r = min(r1,r2)
      R = max(r1,r2)
      if r == R:
         self.__vertex = float('inf')
         self.__direction = p2 - p1
         self.__angle = 0.0
      else:
         L = self.__length*R/(R - r)
         PointList = [p1,p2]
         minInd = np.argmin(np.array([r1,r2]))
         maxInd = np.argmax(np.array([r1,r2]))
         a = PointList[maxInd] - PointList[minInd]
         a = a/la.norm(a)
         self.__vertex = PointList[maxInd] - a*L
         self.__direction = a
         self.__angle = math.atan(R/L)
   
####################################

   def Intersects(self, Ray): # Ray has Position, Direction
      from numpy import linalg as la
      
      # http://geomalgorithms.com/a07-_distance.html
      
      # the closest distance between two lines is through the common perpendicular.
      # Here we find the point on each line that is closest to the other line.
      # we make use of wc*p = 0 = wc*q where wc = P(tc[0]) - Q(tc[1])
      p0 = Ray.Position; p = Ray.Direction
      q0 = self.__vertex; q = self.__direction
      
      A = np.array([[sum(p*p), -sum(q*p)],[sum(p*q), -sum(q*q)]])
      C = np.array([[sum(q0*p - p0*p)],[sum(q0*q - p0*q)]])
      if la.det(A) != 0: tc = np.dot(la.inv(A), C) # Ray[tc[0]], Cone[tc[1]]
      
#      print 'tc = %f, %f' % (tc[0], tc[1])
      
      # now we get (p1,q1) = min(P,Q)
      p1 = p0 + p*tc[0]
      q1 = q0 + q*tc[1]
      
#      print 'cone vertex = %f, %f, %f' % (self.__vertex[0], self.__vertex[1], self.__vertex[2])
#      print 'cone direction = %f, %f, %f' % (self.__direction[0], self.__direction[1], self.__direction[2])
#      print 'at angle = %f' % (self.__angle)
#      print 'angle at minimum dist = %f' % (math.atan2(la.norm(p1-q1), la.norm(q1-q0)))
#      print 'p0 = %f, %f, %f   p = %f, %f, %f\n' % (p0[0],p0[1],p0[2], p[0],p[1],p[2])
      
      if math.atan2(la.norm(p1-q1), la.norm(q1-q0)) <= self.__angle and \
            la.norm(q1-q0) > la.norm(self.p1-q0) and \
            la.norm(q1-q0) < la.norm(self.p2-q0) and \
            tc[1] >= 0:
         return True
      else:
         return False
      
      
####################################

   def Dist(self, Ray):
      p0 = Ray.Position; p = Ray.Direction
      q0 = self.__vertex; q = self.__direction
      theta = self.__angle
      
      # find the time to intersection
      #t = [(sum(p*q)*(math.cos(theta)**2*sum(q*q) - 2*sum(p0*q) + 2*sum(q0*q0)) - math.sqrt(math.cos(theta)**4*sum(p*q)**2*sum(q*q)**2))/(2*sum(p*q)**2), (sum(p*q)*(math.cos(theta)**2*sum(q*q) - 2*sum(p0*q) + 2*sum(q0*q0)) + math.sqrt(math.cos(theta)**4*sum(p*q)**2*sum(q*q)**2))/(2*sum(p*q)**2)]
      try:
         t = [(sum(p*q)*(sum(p0*q) - sum(q0*q)) + math.sqrt(math.cos(theta)**2*sum(q*q)*(-math.cos(theta)**2*sum(p0*p0)*sum(p*p)*sum(q*q) + 2*math.cos(theta)**2*sum(p0*q0)*sum(p*p)*sum(q*q) - math.cos(theta)**2*sum(p*p)*sum(q0*q0)*sum(q*q) + sum(p0*p0)*sum(p*q)**2 + sum(p0*q)**2*sum(p*p) - 2*sum(p0*q)*sum(p*p)*sum(q0*q) - 2*sum(p0*q0)*sum(p*q)**2 + sum(p*p)*sum(q0*q)**2 + sum(p*q)**2*sum(q0*q0))))/(math.cos(theta)**2*sum(p*p)*sum(q*q) - sum(p*q)**2), (sum(p0*q)*sum(p*q) - sum(p*q)*sum(q0*q) - math.sqrt(math.cos(theta)**2*sum(q*q)*(-math.cos(theta)**2*sum(p0*p0)*sum(p*p)*sum(q*q) + 2*math.cos(theta)**2*sum(p0*q0)*sum(p*p)*sum(q*q) - math.cos(theta)**2*sum(p*p)*sum(q0*q0)*sum(q*q) + sum(p0*p0)*sum(p*q)**2 + sum(p0*q)**2*sum(p*p) - 2*sum(p0*q)*sum(p*p)*sum(q0*q) - 2*sum(p0*q0)*sum(p*q)**2 + sum(p*p)*sum(q0*q)**2 + sum(p*q)**2*sum(q0*q0))))/(math.cos(theta)**2*sum(p*p)*sum(q*q) - sum(p*q)**2)]
      except: t = [0.0, 0.0]
      
      hitPnt = p0 + p*min(t) # need to consider negative distances later
      return la.norm(hitPnt - p0)
            
   def Dist2(self, Ray):
      # find the time to intersection first with Newtons method
      t = 0.0
      for i in range(8):
         if self.__fp(t, Ray) != 0:
            t -= self.__f(t, Ray)/self.__fp(t, Ray)
      
      return la.norm(Ray.Position + Ray.Direction*t)
      
####################################

####################################
   def Reflection(self, Ray, LightSource):
      p0 = Ray.Position; p = Ray.Direction
      q0 = self.__vertex; q = self.__direction
      theta = self.__angle
      
      # find the time to intersection
      #t = [(sum(p*q)*(math.cos(theta)**2*sum(q*q) - 2*sum(p0*q) + 2*sum(q0*q0)) - math.sqrt(math.cos(theta)**4*sum(p*q)**2*sum(q*q)**2))/(2*sum(p*q)**2), (sum(p*q)*(math.cos(theta)**2*sum(q*q) - 2*sum(p0*q) + 2*sum(q0*q0)) + math.sqrt(math.cos(theta)**4*sum(p*q)**2*sum(q*q)**2))/(2*sum(p*q)**2)]
      try:
         t = [(sum(p*q)*(sum(p0*q) - sum(q0*q)) + math.sqrt(math.cos(theta)**2*sum(q*q)*(-math.cos(theta)**2*sum(p0*p0)*sum(p*p)*sum(q*q) + 2*math.cos(theta)**2*sum(p0*q0)*sum(p*p)*sum(q*q) - math.cos(theta)**2*sum(p*p)*sum(q0*q0)*sum(q*q) + sum(p0*p0)*sum(p*q)**2 + sum(p0*q)**2*sum(p*p) - 2*sum(p0*q)*sum(p*p)*sum(q0*q) - 2*sum(p0*q0)*sum(p*q)**2 + sum(p*p)*sum(q0*q)**2 + sum(p*q)**2*sum(q0*q0))))/(math.cos(theta)**2*sum(p*p)*sum(q*q) - sum(p*q)**2), (sum(p0*q)*sum(p*q) - sum(p*q)*sum(q0*q) - math.sqrt(math.cos(theta)**2*sum(q*q)*(-math.cos(theta)**2*sum(p0*p0)*sum(p*p)*sum(q*q) + 2*math.cos(theta)**2*sum(p0*q0)*sum(p*p)*sum(q*q) - math.cos(theta)**2*sum(p*p)*sum(q0*q0)*sum(q*q) + sum(p0*p0)*sum(p*q)**2 + sum(p0*q)**2*sum(p*p) - 2*sum(p0*q)*sum(p*p)*sum(q0*q) - 2*sum(p0*q0)*sum(p*q)**2 + sum(p*p)*sum(q0*q)**2 + sum(p*q)**2*sum(q0*q0))))/(math.cos(theta)**2*sum(p*p)*sum(q*q) - sum(p*q)**2)]
      except: t = [0.0, 0.0]
      t = min(t)
      
      # get the vector from the vertex to the intersecting point of the cone
      Point = p0 + p*t
      u = Point - q0
      
      # get the vector from the vertex of the cone to the point along the cone direction for finding the normal
      w = q/la.norm(q)*la.norm(u)/math.cos(theta)
      
      # get the normal
      normal = Point - (q0 + w)
      normal = normal/la.norm(normal)
      
      UnScatRay = copy.copy(Ray)
      UnScatRay.Bounce(Point, normal)
      
      ScatRay = copy.copy(Ray)
      ScatRay.Scatter(Point, normal, LightSource, UnScatRay)
      
      ScatRay.Energy = 1.0; UnScatRay.Energy = 1.0; #
      return [UnScatRay, ScatRay] # unscattered ray, scattered ray (toward light source)
      

   def Reflection2(self, Ray, LightSource):
      # find the time to intersection first with Newtons method
      t = 0.0
      for i in range(8):
         if self.__fp(t, Ray) != 0:
            t -= self.__f(t, Ray)/self.__fp(t, Ray)
      
      # get the vector from the vertex to the intersecting point of the cone
      Point = Ray.Position + Ray.Direction*t
      u = Point - self.__vertex
      
      # project u onto the direction of the cone
      q = sum(u*self.__direction)*self.__direction
      
      # get the normal vector to the cone
      w = Point-q
      qp = (la.norm(q) + la.norm(w)*math.tan(self.__angle))*self.__direction
      normal = (Point-self.__vertex) - qp
      normal = normal/la.norm(normal)
      
      UnScatRay = copy.copy(Ray)
      UnScatRay.Bounce(Point, normal)
      
      ScatRay = copy.copy(Ray)
      ScatRay.Scatter(Point, normal, LightSource, UnScatRay)
      
 #     ScatRay.Energy = 1.0; UnScatRay.Energy = 1.0; #
      return [UnScatRay, ScatRay] # unscattered ray, scattered ray (toward light source)
      
####################################

   def __f(self, t, Ray):
      RayV = Ray.Position - self.__vertex + Ray.Direction * t
      return sum(RayV*self.__direction)/la.norm(RayV)/la.norm(self.__direction)
      
####################################

   def __fp(self, t, Ray):
      RayV = Ray.Position - self.__vertex + Ray.Direction * t
      return sum(Ray.Direction*self.__direction)/la.norm(RayV)/la.norm(self.__direction)
   
      
