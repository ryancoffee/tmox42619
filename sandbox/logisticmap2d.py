#!/usr/bin/python3

import numpy as np
import sys
import h5py

def getlogistic(x0,r):
    tol = np.finfo(float).resolution
    x = np.copy(x0)
    for i in range(2**8):
        xprev = np.copy(x)
        x = r*xprev*(1-xprev)
        if np.abs(x-xprev)<tol:
            return x
    return x

def getlogistic2d(x0,y0,r,c):
    x = np.copy(x0)
    y = np.copy(y0)
    tol = np.finfo(float).resolution
    for i in range(2**6):
        xprev = np.copy(x)
        yprev = np.copy(y)
        xp=getlogistic(xprev,r)
        yp=getlogistic(yprev,r)
        x = (1-c)*xp + (c)*yp
        y = (1-c)*yp + (c)*xp
        if max(np.abs(x-xprev),np.abs(y-yprev))<tol:
            return x,y
    return x,y

def val2key(val):
    strfmt = np.log10(np.finfo(type(val)).resolution)
    decprec='%i'%( int( abs( np.log10(np.finfo(type(val)).resolution) ) ) )
    fmtstring = '%.' + decprec + 'f'
    return fmtstring%val

def getNpoints(n,r,c):
    xx = []
    yy = []
    for i in range(n):
        x,y = getlogistic2d(np.random.random(),np.random.random(),r,c)
        xx += [x]
        yy += [y]
    return xx,yy

def main():
    if len(sys.argv)<5:
        print('Give me the follwoing:' )
        print('%s <npoints> <rval> <cval> <outfile.h5>'%sys.argv[0] )
        return
    npoints = int(sys.argv[1])
    rval = np.float32(sys.argv[2])
    cval = np.float32(sys.argv[3])
    ofname = sys.argv[4]
    with h5py.File(ofname,'a') as f:
        print(rval,cval)
        x=[]
        y=[]
        rgrp = f.create_group(val2key(rval))
        cgrp = rgrp.create_group(val2key(cval))
        x,y = getNpoints(npoints,rval,cval)
        cgrp.create_dataset('x',data=x)
        cgrp.create_dataset('y',data=y)
    return

if __name__ == '__main__':
    main()
