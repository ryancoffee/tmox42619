#!/usr/bin/python3

import h5py
import numpy as np
import sys
from scipy.interpolate import interp1d 

def fill2dhist(csum):
    x = csum
    y = np.arange(x.shape[1],dtype=np.uint8) # inflating by facotr of 4
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

    for fname in fnames:
        #with h5py.File('./data_fs/h5files/%s'%(fname),'r') as f:
        with h5py.File('%s'%(fname),'r') as f:
            print(list(f.keys()))
            for key in f.keys():
                data[key] = np.squeeze(f[key][()])
        # here now the data is stored in system memory and we close the h5 file
        vlscsum = np.cumsum(data['vls'],axis=1)
        h = fill2dhist(vlscsum)

        np.savetxt('./data_fs/ascii/origvls.dat',data['vls'][:10,:].T,fmt='%i')
        np.savetxt('./data_fs/ascii/tmpvls.dat',vlscsum[:10,:].T,fmt='%i')
        np.savetxt('./data_fs/ascii/hist.dat',h.T,fmt='%i')
    return


if __name__ == '__main__':
    main()
