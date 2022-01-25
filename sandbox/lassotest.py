#!/usr/bin/python3


from scipy.fftpack import dct, idct
from sklearn.linear_model import Lasso
import scipy.fftpack as spfft
import numpy as np
import matplotlib.pyplot as plt
import math


def gauss(cx,cy,wx,wy,szx,szy):
  x = np.arange(szx)
  y = np.arange(szy)
  zx = np.exp(-((x-cx)/wx)**2)
  zy = np.exp(-((y-cy)/wy)**2)
  return np.outer(zx,zy)

def genImage(nc=5,nx=32,ny=32):
  rng = np.random.default_rng(5)
  xcenters = rng.uniform(low=0.,high=float(nx),size=nc)
  ycenters = rng.uniform(low=0.,high=float(ny),size=nc)
  xwidths = rng.poisson(2.,size=nc)+0.5
  ywidths = rng.poisson(3.,size=nc)+0.5
  X = np.zeros((nx,ny),dtype=float)
  #print(xcenters,ycenters,xwidths,ywidths)
  for i in range(nc):
    X = X + gauss(xcenters[i],ycenters[i],xwidths[i],ywidths[i],nx,ny) + 1.
  return X

def idct2(x):
  return spfft.idct(spfft.idct(x.T, norm='ortho', axis=0).T, norm='ortho', axis=0)

def idct2_rnc(x):
  return spfft.idct(spfft.idct(x, axis=0),axis=1)


def main():
    alphalist = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]
    nx=ny=16
    X = genImage(nc=4,nx=nx,ny=ny)
    k = round(nx * ny * 0.5) # 50% sample
    ri = np.random.choice(nx * ny, k, replace=False) # random sample of indices
    print(X)
    #b = X.T.flat[ri]
    #b = np.expand_dims(b, axis=1)
    #blist = [b for i in range(len(alphalist))]
    blist = [np.expand_dims(X.T.flat[ri],axis=1) for i in range(len(alphalist))]
    Alist = [ np.kron(
        spfft.idct(np.identity(nx), norm='ortho', axis=0),
        spfft.idct(np.identity(ny), norm='ortho', axis=0)
        ) for i in range(len(alphalist))]
    for i in range(len(alphalist)):
        Alist[i] = Alist[i][ri,:] # same as phi times kron

    lassolist = []
    lassolist += [Lasso(a) for a in alphalist]
    _ = [lassolist[i].fit(Alist[i],blist[i]) for i in range(len(lassolist))]

    Xatlist = [np.array(l.coef_).reshape(nx, ny).T for l in lassolist] # stack columns
    #Xat[:, ny//2:ny] = 0
    # Get the reconstructed image
    print(np.max(Xatlist[0]))
    Xalist = [idct2_rnc(Xat) for Xat in Xatlist]
    fig = plt.figure()
    imlist = []
    Xatlist = []
    Xalist = []
    for i in range(len(Xalist)):
        imlist += ['im'%i]
        Xatlist += ['Xat%i'%(i)]
        Xalist += ['Xa%i'%(i)]
        axstr = '%i%i%i'%(len(Xalist),3,i+1)
        imlist += fig.add_subplot(axstr)
        axstr = '%i%i%i'%(len(Xalist),3,i+2)
        Xatlist += fig.add_subplot(axstr)
        axstr = '%i%i%i'%(len(Xalist),3,i+3)
        Xalist += fig.add_subplot(axstr)
        imlist[-1].imshow(X)
        Xatlist[-1].imshow(Xatlist[i])
        Xalist[-1].imshow(Xalist[i])
    return

if __name__ == '__main__':
    main()
