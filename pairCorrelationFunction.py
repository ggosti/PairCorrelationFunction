import numpy as np
import matplotlib.pyplot as plt

from skimage import io
from skimage.measure import regionprops
from scipy.spatial import distance



def find_centers(labels):
  h,w = labels.shape
  regions = regionprops(labels)
  xs = []
  ys = []
  for props in regions:
    y0, x0 = props.centroid
    xs.append(x0)
    ys.append(y0)
  xs = np.array(xs)
  ys = np.array(ys)
  return(xs,ys)

def makeBox(xs,ys, bs = 1000):
  #rmax = bs #1000
  centers = np.stack((centersX, centersY), axis = 0).transpose()
  insideBox = (bs<=centers[:,1])*((h-bs)>=centers[:,1])*(bs<=centers[:,0])*((w-bs)>=centers[:,0])
  xsBox,ysBox = centersX[insideBox2], centersY[insideBox2]
  return insideBox,xsBox,ysBox

def distances(xs,ys,xsBox,ysBox):
  #n = len(xs)
  #nBox = len(xsBox)
  #rs = []
  #for i,xi,yi in zip(np.arange(n),xsBox,ysBox):
  #  if i%500==0: print(i)
  #  Dx = xs - xi  
  #  Dy = ys - yi
  #  if (Dx ** 2 < 0).any(): print('Dx min 0', Dx)
  #  if (Dy ** 2 < 0).any(): print('Dx min 0', Dy)
  #  if ((Dx ** 2 + Dy ** 2) < 0).any(): print('SumofSq min 0', (Dx ** 2 + Dy ** 2))
  #  distances = np.sqrt(Dx ** 2 + Dy ** 2)
  #  rs = rs + list(distances)#[distances < 2*bs]) #bs*Conv]) 
  centersA = np.stack((xs, ys), axis = 0).transpose()
  centersB = np.stack((xsBox,ysBox), axis = 0).transpose()
  R = distance.cdist(centersA, centersB, 'euclidean')
  rs = R.flatten()
  return rs[rs < 2*bs]  

def gdr(xs,xsBox, rs, bs = 1000): 
  n = len(xs)
  nBox = len(xsBox)  
  dr = 3 #(pixel)
  hist, bin_edges = np.histogram(rs, bins= np.arange(1,bs, dr),normed=False)
  rhoBox = float(nBox) / float((w - 2 * bs) * (h - 2 * bs)  ) #(micron)
  rho = float(n) / float(w * h  )  #(micron)
  r = bin_edges[:-1] + 0.5 * dr
  deltaV = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)
  g_tot = hist / ( nBox * rho * deltaV )
  return r, g_tot
