import math
import numpy as np

def dihangl(a, ix, jx, kx, lx, noneg):
  """
  
  Returns the dihedral angle of atoms of
  atoms a[ix], a[jx], a[kx], a[lx]

  ***
  Parameters:
  a: input coordinates
  ix, jx, kx, lx : index

  Returns:
  dihangle: dihedral angle

  """

  PI = math.pi 
  dcoszero = 0.99999
  print("DIHANGL ix = {}, jx = {}, kx = {}, lx = {}".format(ix,jx,kx,lx))
  for k in range(3):
    ai[k] = a[k][ix]
    aj[k] = a[k][jx]
    ak[k] = a[k][kx]
    al[k] = a[k][lx]
  rjk = arrdistsd(ak,aj)
  rij = arrdistsd(ai,aj)  
  rkl = arrdistsd(ak,al)
  dij = ai - aj
  dkj = ak - aj
  dlk = al - ak
  cosijk1 = np.dot(dij, dkj)
  cosjkl1 = np.dot(dij, dkj)
  cosijk = cosijkl/(rij*rjk)
  cosjkl = cosjkl1/(rkl*rjk)
  if (cosijk < -dcoszero or cosijk > dcoszero or cosjkl < - dcoszero or cosjkl > dcoszero):
    print("DIHANGL: colinear atoms, torision set to 0 cosijk =, jkl",cosijk,cosjkl)
    return 0
  costa = 0
  for k in range(3):
    dii[k] = dij[k] - dkj[k]*cosijk1/rjk2
    djj[k] = dlk[k] + dkj[k]/cosjkl/rjk2
    costa += dii[k]*djj[k]
  costa = costa/(rij*rkl*math.sqrt((1.0-cosijk*cosijk)*(1.0-cosjkl*cosjkl)))
  try:
    angle = dacoscheck(costa)
  except
    print("Problem in dacoscheck costa = ",costa)
  # Decide sign of the angle 
  cp = dvprd(dii, dkj)
  if ((cp[0]*djj[0]+cp[1]*djj[1] + cp[2]*djj[2]) > 0.0):
    angle = - angle
  if (angle < 0 and noneg == 1):
    angle += 2*PI 
  dihangl = angle
  return dihangl  
  