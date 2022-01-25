#!/usr/bin/python3

import numpy as np
import sys

def getlogistic(x0,r,tol=0.0000001):
    x = np.copy(x0)
    for i in range(2**8):
        xprev = np.copy(x)
        x = r*xprev*(1-xprev)
        if np.abs(x-xprev)<tol:
            return x
    return x

def getlogistic2d(x0,y0,r,c):
    tol = min(np.finfo(x0).resolution,np.finfo(y0).resolution)
    x = np.copy(x0)
    y = np.copy(y0)
    for i in range(2**8):
        xprev = np.copy(x)
        yprev = np.copy(y)
        x = (1-c)*xprev + (c)*yprev
        y = (1-c)*yprev + (c)*xprev
        if np.abs(x-xprev)<tol && np.abs(y-yprev)<tol:
            return x,y
    return x,y

def main():
    if len(sys.argv)<5:
        print('Give me the number of bins and then the range you want low' )
        print('%s <xbins> <rbins> <lowr> <highr>'%sys.argv[0] )
        return
    nxbins = int(sys.argv[1])
    nrbins = int(sys.argv[2])
    lowr = float(sys.argv[3])
    highr = float(sys.argv[4])
    for r in lowr+(highr-lowr)*np.random.random(nrbins):
        f = open('coords_rx_%.3f_%.3f.out'%(lowr,highr),'a')
        print(r)
        for i in range(nxbins):
            f.write('%.6f\t%.6f\n'%(r,getlogistic(np.random.random(),r)))
        f.close()
    return

if __name__ == '__main__':
    main()
