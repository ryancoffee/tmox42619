#!/usr/bin/python3


from scipy.fftpack import dct
from sklearn.linear_model import Lasso
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
  X = np.ones((nx,ny),dtype=int)
  #print(xcenters,ycenters,xwidths,ywidths)
  for i in range(nc):
    X = X + (2**4*gauss(xcenters[i],ycenters[i],xwidths[i],ywidths[i],nx,ny)).astype(int)
  return X

def dct_2d(x):
    return(dct(dct(x,axis=0,type=2),axis=1,type=2))

def idct_2d(x):
  return dct(dct(x,axis=0,type=3), axis=1,type=3)/4/((x.shape[0]*x.shape[1]))

def main():
    alphalist = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]
    nx=ny=32
    X = genImage(nc=8,nx=nx,ny=ny)
    k = round(nx * ny * 0.5) # 50% sample
    np.savetxt('testImg.dat',X,fmt='%i')
    Y = dct_2d(X)
    Y00 = Y[0][0]
    Y[0,0] = 0
    np.savetxt('testImg_dct.dat',Y,fmt='%i')
    np.savetxt('testImg_idct.dat',idct_2d(Y),fmt='%i')
    return
'''
    ri = np.random.choice(nx * ny, k, replace=False) # random sample of indices
    print(X)
    #b = X.T.flat[ri]
    #b = np.expand_dims(b, axis=1)
    #blist = [b for i in range(len(alphalist))]
    b = np.expand_dims(X.T.flat[ri],axis=1)
    A = np.kron(
        spfft.idct(np.identity(nx), norm='ortho', axis=0),
        spfft.idct(np.identity(ny), norm='ortho', axis=0)
        )
    Asparse = A[ri,:] # same as phi times kron

    lassolist = [Lasso(a) for a in alphalist]
    lassolist = [lasso.fit(Asparse,b) for lasso in lassolist]

    Xatlist = [np.array(l.coef_).reshape(nx, ny).T for l in lassolist] # stack columns
    #Xat[:, ny//2:ny] = 0
    # Get the reconstructed image
    Xalist = [idct2_rnc(Xat) for Xat in Xatlist]
    fig = plt.figure()
    ax =[]
    ax += [fig.add_subplot(331)]
    ax += [fig.add_subplot(332)]
    ax += [fig.add_subplot(333)]
    ax += [fig.add_subplot(334)]
    ax += [fig.add_subplot(335)]
    ax += [fig.add_subplot(336)]
    ax += [fig.add_subplot(337)]
    ax += [fig.add_subplot(338)]
    ax += [fig.add_subplot(339)]
  
    for i in range(3):#len(Xalist)):
        print(i)
        ax[i*3].imshow(X)
        ax[i*3+1].imshow(Xatlist[i])
        ax[i*3+2].imshow(Xalist[i])
    fig.show()
    return
'''

if __name__ == '__main__':
    main()
