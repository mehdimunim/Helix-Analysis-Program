import math

def arrdistsd(a1, a2):
      """
      
      dist2 = sum(ai - bi)**2

      """
      x = a1[0] - a2[0]
      y = a1[1] - a2[1]
      z = a1[2] - a2[2]
      dist2 = x*x + y*y + z*z
      dist = math.sqrt(dist2)
      return dist