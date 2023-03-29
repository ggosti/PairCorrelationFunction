import numpy as np
import matplotlib.pyplot as plt

from skimage import io
from skimage.measure import regionprops
from scipy.spatial import distance



def find_centers(labels):
  """Given masked image 'labels', returns the centers of the masks."""
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

def makeBox(xs,ys, h, w, bs = 1000):
  """Identidies positions of particles inside the image. Inside is defined ignoring a 
  border of 'bs' pixels from sides, top and bottom of image. 'xs' and 'ys' are the x and y positions of the centers, 
  respectively. 'h' and 'w' are height and width of the image."""
  #rmax = bs #1000
  centers = np.stack((xs, ys), axis = 0).transpose()
  insideBox = (bs<=centers[:,1])*((h-bs)>=centers[:,1])*(bs<=centers[:,0])*((w-bs)>=centers[:,0])
  xsBox,ysBox = xs[insideBox], ys[insideBox]
  return insideBox,xsBox,ysBox

def distances(xs,ys,xsBox,ysBox)#, bs = 1000):
  """Calculates distances between centers inside the box and all other centers. 
  'xs' and 'ys' are the x and y positions of all the centers, 
  'xsBox' and 'ysBox' are the x and y positions of the centers inside the central box."""
  centersA = np.stack((xs, ys), axis = 0).transpose()
  centersB = np.stack((xsBox,ysBox), axis = 0).transpose()
  R = distance.cdist(centersA, centersB, 'euclidean')
  rs = R.flatten()
  return rs #[rs < 2*bs]  

def gdr(xs,xsBox, w, h, rs, bs = 1000, dr = 3):
  """Function to calculate the gdr. 'xs' and 'xsBox' are the x coordinates of centers in the image and in the inside Box, respectively. 'w' and 'h' are width and height of the image, 'bs' is the size of the border. 'dr' is the number of pixel in a bin of the histogram, i.e. the resolution of the gdr."""
  n = len(xs)
  nBox = len(xsBox)  
  hist, bin_edges = np.histogram(rs, bins= np.arange(1,bs, dr),normed=False)
  rhoBox = float(nBox) / float((w - 2 * bs) * (h - 2 * bs)  ) #(micron)
  rho = float(n) / float(w * h  )  #(micron)
  r = bin_edges[:-1] + 0.5 * dr
  deltaV = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)
  g_tot = hist / ( nBox * rho * deltaV )
  return r, g_tot

def calc_gdr(labels, b = 1000):
  """Function to calculate the gdr. Takes masked image 'label' and boarder size 'b' as arguments."""
  h,w = labels.shape
  xs, ys = find_centers(labels)
  temp, xsBox, ysBox = makeBox(xs,ys, h, w, bs = b)
  rs = distances(xs,ys,xsBox,ysBox)
  return gdr(xs,xsBox, w, h, rs, bs = b)
