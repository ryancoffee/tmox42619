#!/usr/bin/python3

import numpy as np
import h5py
import math
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


def randomround(x:float,rng):
    return (np.int64(x) + np.int64(x%1>rng.random()))


def checkdet(runlist,detname):
    for r in runlist:
        if not detname in r.detnames:
            return False 
    return True

def PWRspectrum(wv):
    return np.power(abs(np.fft.fft(wv).real),int(2))

def rollon(vec,n):
    vec[:int(n)] = vec[:int(n)]*np.arange(int(n),dtype=float)/float(n)
    return vec

def tanhFloat(x,thresh=1):
    y = 0.5*(1.0 + np.tanh(x.astype(float)/thresh))
    return y.astype(type(x[0]))

def tanhInt(x,bits=8):
    y = (1<<(bits-1))*(1+np.tanh(x.astype(float)))
    return y.astype(type(x[0]))

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

