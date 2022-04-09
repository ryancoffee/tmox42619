#!/usr/bin/python3

import numpy as np
import h5py

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

def makehist(fname,bins=(np.arange(2**10+1)-2**9)*2**10):
    with h5py.File(fname,'r') as f:
        for p in f.keys():
            data = []
            for k in f[p]['waves'].keys():
                data += list(f['port_0']['waves'][k][()])
            h,b = np.histogram(data,bins)
            np.savetxt('histsig.%s.dat'%p,np.column_stack((b[:-1],h)),fmt='%i')
    return

def samples2ascii(ipath,fname,opath):
    logics={}
    waves={}
    with h5py.File('%s/%s'%(ipath,fname),'r') as f:
        for p in f.keys(): 
            logics[p] = np.stack([f[p]['logics'][k][()] for k in f[p]['logics'].keys()],axis=-1)
        for p in f.keys(): 
            waves[p] = np.stack([f[p]['waves'][k][()] for k in f[p]['waves'].keys()],axis=-1)
        _=[np.savetxt('%s/%s.logics.%s.dat'%(opath,fname,p),logics[p],fmt='%i') for p in logics.keys()]
        _=[np.savetxt('%s/%s.waves.%s.dat'%(opath,fname,p),waves[p],fmt='%i') for p in waves.keys()]
    return

