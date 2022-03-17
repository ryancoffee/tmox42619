#!/usr/bin/python3

import numpy as np
from scipy.fftpack import dct,dst

def show(x):
    offset = 5+np.abs(np.min(x))
    print(' '*(offset-1)+'|')
    _ = [print(' '*int(offset+v-1)+'*') for v in x]
    print(' '*(offset-1)+'|')
    return

def dctmat(sz):
    return 2*dct( np.identity(2*sz),axis=1 , type=2)[::2].T
def dstmat(sz):
    return 2*dst( np.identity(2*sz),axis=1 , type=2)[1::2].T
def idctmat(sz):
    return 1./2*dct( np.identity(2*sz),axis=0 , type=3)[::2].T
def idstmat(sz):
    return 1./2*dst( np.identity(2*sz),axis=0 , type=3)[1::2].T

def mydct(cmat,x):
    xin = np.append(x[::2],x[-1::-2])
    return np.inner(cmat,xin)[::2]

def mydst(smat,x):
    xin = np.append(x[1::2],-1*x[-2::-2])
    return np.inner(smat,xin)[1::2]

# then add the frequency scaled version
# then add the zeros append for oversampling


def pkey(p):
    return 'port_%i'%p

def mypoly(x,order=4):
	result = np.ones((x.shape[0],order+1),dtype=float)
	result[:,1] = x.copy()
	if order < 2:
		return result
	for p in range(2,order+1):
		result[:,p] = np.power(result[:,1],int(p))
	return result

def fitpoly(x,y,order=4):
    assert len(x)==len(y)
    assert len(x)>order
    x0 = np.mean(np.array(x))
    theta = np.linalg.pinv( mypoly(np.array(x-x0).astype(float),order=order) ).dot(np.array(y).astype(float)) # fit a polynomial (order 3) to the points
    return x0,theta

def getcentroid(inds,spec):
    x = inds
    y = spec 
    num = np.sum(x*y)
    denom = np.sum(y)
    if (denom>0):
        return int(num/denom)
    return 0
