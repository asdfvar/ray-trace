import Ray
import Shape
import numpy as np
import copy

a = Ray.Ray(7,7,7)

print a.Position

b = copy.copy(a)

b.Position = 8

print b.Position

print a.Position
