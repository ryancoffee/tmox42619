#!/usr/bin/python3

import numpy as np
import sys
import h5py

def getlogistic(x0,r):
    x = np.copy(x0)
    tol=np.finfo(x0).resolution
    for i in range(2**8):
        xprev = np.copy(x)
        x = r*xprev*(1-xprev)
        if np.abs(x-xprev)<tol:
            return x
    return x

def rval2key(rval):
    strfmt = np.log10(np.finfo(type(rval)).resolution)
    decprec='%i'%( int( abs( np.log10(np.finfo(type(rval)).resolution) ) ) )
    fmtstring = '%.' + decprec + 'f'
    return fmtstring%rval

def getNxvals(n,r):
    x = []
    for i in range(n):
        x += [getlogistic(np.float32(np.random.random()),r)]
    return x


def main():
    if len(sys.argv)<6:
        print('Give me the number of xvals, xbins, rvalue and outfilename' )
        print('%s <nxvals> <nxbins> <rval> <cval> <outfilename>'%sys.argv[0] )
        return
    nxvals = int(sys.argv[1])
    nxbins = int(sys.argv[2])
    rval = np.float32(sys.argv[3])
    cval = np.float32(sys.argv[4])
    ofname = '%s'%(sys.argv[5])
    with h5py.File(ofname,'a') as f:
        rkey = rval2key(rval)
        rgrp = None
        if rkey in f.keys():
            rgrp = f[rkey]
            if 'xlist' not in rgrp.keys():
                xlist = getNxvals(nxvals,rval)
                rgrp.create_dataset('xlist',data=xlist,dtype=np.float64,maxshape=(None,))
            else:
                rgrp['xlist'][()] += getNxvals(nxvals,rval)
            h,b = np.histogram(rgrp['xlist'][()],nxbins)
            if 'hist' not in rgrp.keys():
                rgrp.create_dataset('hist',data=h,dtype=int,maxshape=(None,))
                rgrp['hist'].attrs.create('min',b[0])
                rgrp['hist'].attrs.create('max',b[-1])
                rgrp['hist'].attrs.create('ncounts',len(rgrp['xlist'][()]))
            else:
                for i in range(len(h)):
                    rgrp['hist'][i] = h[i]
                rgrp['hist'].attrs['min'] = b[0]
                rgrp['hist'].attrs['max'] = b[-1]
                rgrp['hist'].attrs['ncounts'] = len(rgrp['xlist'][()])
        else:
            rgrp = f.create_group(rval2key(rval))
            xlist = getNxvals(nxvals,rval)
            rgrp.create_dataset('xlist',data=xlist,dtype=np.float64,maxshape=(None,))
            h,b = np.histogram(xlist,nxbins)
            rgrp.create_dataset('hist',data=h,dtype=int,maxshape=(None,))
            rgrp['hist'].attrs.create('min',b[0])
            rgrp['hist'].attrs.create('max',b[-1])
            rgrp['hist'].attrs.create('ncounts',len(rgrp['xlist'][()]))
    return

if __name__ == '__main__':
    main()
