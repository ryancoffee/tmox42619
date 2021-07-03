#!/usr/bin/python3

import h5py
import numpy as np
import sys
from scipy.interpolate import interp1d 
import re

def fill2dhist(csum):
    x = csum
    y = np.arange(x.shape[1],dtype=np.uint8) 
    h = np.zeros(x.shape,dtype=np.int8)
    for i in range(x.shape[0]):
        cs = {}
        for j in range(x.shape[1]):
            cs.update({x[i,j]:y[j]})
        xx = [v for v in cs.keys()]
        yy = [v for v in cs.values()]
        f = interp1d(xx,yy)
        xnew = np.random.random((100,))*(xx[-2]-xx[1]) + xx[1]
        ynew = f(xnew).astype(np.uint8)
        for k in ynew:
            h[i,k] += 1
    return h

def main():
    if len(sys.argv)<2:
        print('syntax: make2Dhist.py <fnames>')
        return
    fnames = list(sys.argv[1:])
    data = {}
    fullname = sys.argv[1]
    datadir = './'
    filefront = 'hits.h5'

    for fname in fnames:
        #with h5py.File('./data_fs/h5files/%s'%(fname),'r') as f:
        m = re.search('^(.+)/h5files/(.+).h5',fname)
        if m:
            fullname = m.group(0)
            datadir = m.group(1)
            filefront = m.group(2)
            print('processing:\t%s'%fullname)
        else:
            return
        with h5py.File('%s'%(fullname),'r') as f:
            print(list(f.keys()))
            for key in f.keys():
                data[key] = np.squeeze(f[key][()])
        # here now the data is stored in system memory and we close the h5 file
        vlsmean = np.mean(data['vls'],axis=0)
        vlscsum = np.cumsum(data['vls'],axis=1)
        h = fill2dhist(vlscsum)
        tof = data['port_12_tofs']
        tofaddresses = data['port_12_addresses']
        tofnedges= data['port_12_nedges']
        print(tof.shape)
        print('\t'.join((tof[:20]-25000).astype(str)))
        print('\t'.join(tofaddresses[:20].astype(str)))
        print('\t'.join(tofnedges[:20].astype(str)))

        tofhist = np.zeros((2**16,vlsmean.shape[0]),dtype=np.int8)
        print(tofhist.shape)
        return

        x = vlscsum
        y = np.arange(x.shape[1],dtype=np.uint8) 


        np.savetxt('%s/ascii/%s.origvls.dat'%(datadir,filefront),data['vls'][:10,:].T,fmt='%i')
        np.savetxt('%s/ascii/%s.tmpvls.dat'%(datadir,filefront),vlscsum[:10,:].T,fmt='%i')
        np.savetxt('%s/ascii/%s.hist.dat'%(datadir,filefront),h.T,fmt='%i')
    return


if __name__ == '__main__':
    main()
