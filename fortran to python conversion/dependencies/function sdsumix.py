import numpy as np

def sdsumix(n, c1,c2,atw,numsel,index,devmax,maxat):
  """
  Calculate the RMS over the atomlist in index, weighted with aw
  """
  sd=0.0
  awsum=0.0
  devmax=0.0
  if (numsel == 0):
    #Use all atoms
    for i in range(n):
      sd0 = (c1[0][i] - c2[0][i])**2 + (c1[1][i] - c2[1][i])**2 + (c1[2][i] - c2[2][i])**2
      sd+=sd0*atw[i]
      awsum+=atw[i]
      if (devmax < sd0):
        devmax = sd0
  else:
    for ii in range(numsel):
      i = index[ii]
      sd0 = (c1[0][i] - c2[0][i])**2 + (c1[1][i] - c2[1][i])**2 + (c1[2][i] - c2[2][i])**2
      sd+=sd0*atw[i]
      awsum+=atw[i]
      if (devmax < sd0):
        devmax=sd0
  devmax = math.sqrt(devmax)
  sdsumix = sd/awsum
  return sdsumix