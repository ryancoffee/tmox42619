#!/usr/bin/python3

import h5py
import numpy as np
import sys

from utils import mypoly,fitpoly

def main():
    if len(sys.argv)<2:
        print('syntax:%s <calib file> <portnum>'%(sys.argv[0]))
        return
    indices = {}
    energies = {}
    for p in (0,1,2,4,5,12,13,14,15):
        indices.update( {'port_%i'%p:[]} )
        energies.update( {'port_%i'%p:[]} )
    with h5py.File(sys.argv[1],'r+') as f:
        for key in f.keys():
            for port in indices.keys():
                indices[port] += [np.log2(float(v)) for v in f[key][port]['indices'][()]]
                energies[port] += [np.log2(float(v)) for v in f[key][port]['energies'][()]]
        for key in f.keys():
            for port in indices.keys():
                x0,theta=fitpoly(indices[port],energies[port],order=2)
                f[key][port].attrs['x0'] = x0
                f[key][port].attrs['theta'] = theta
                if key == 'vret_50':
                    xout = np.arange(-2.,2.,.1)
                    X = mypoly(xout,order=2)
                    Y = X.dot(f[key][port].attrs['theta'])
                    np.savetxt('tmp.fit.%s.dat'%port,np.column_stack((xout+f[key][port].attrs['x0'],Y)),fmt='%.3f')
                    np.savetxt('tmp.%s.dat'%port,np.column_stack((indices[port],energies[port])),fmt='%.3f')
    print('finished')
    return

if __name__ == '__main__':
    main()

