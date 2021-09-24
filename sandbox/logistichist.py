#!/usr/bin/python3

import numpy as np
import sys
import h5py

def getlogistic(x0,r,tol=0.0000001):
    x = np.copy(x0)
    for i in range(2**8):
        xprev = np.copy(x)
        x = r*xprev*(1-xprev)
        if np.abs(x-xprev)<tol:
            return x
    return x

def main():
    if len(sys.argv)<6:
        print('Give me the number of bins and then the range you want low' )
        print('%s <xbins> <rbins> <lowr> <highr> <outfilename>'%sys.argv[0] )
        return
    nxbins = int(sys.argv[1])
    nrbins = int(sys.argv[2])
    lowr = float(sys.argv[3])
    highr = float(sys.argv[4])
    hbins = np.arange(2**10,dtype=float)
    d = {}
    ofname = '%s.%.3f_%.3f.h5'%(sys.argv[5],lowr,highr)
    with h5py.File(ofname,'w') as f:
        f.create_dataset('xbins',data=hbins,dtype=np.float16)
        for r in lowr+(highr-lowr)*np.random.random(nrbins):
            print(r)
            xlist = []
            for i in range(nxbins):
                xlist += [getlogistic(np.random.random(),r)]
            xmin = np.min(xlist)
            xmax = np.max(xlist)
            xlist -= xmin
            xlist *= np.power(2.,int(10))/(xmax-xmin)
            f.create_dataset('%.6f'%r,data=np.histogram(xlist,hbins)[0])
            d.update( {r:np.histogram(xlist,hbins)[0]} )
    
    return

if __name__ == '__main__':
    main()
