#!/usr/bin/python3

import numpy as np

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
