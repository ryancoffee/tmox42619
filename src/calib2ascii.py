#!/usr/bin/python3

import h5py
import numpy as np
import sys

from utils import mypoly,fitpoly

def main():
    fitorder = 2
    if len(sys.argv)<2:
        print('syntax:%s <calib file> <fitorder> '%(sys.argv[0]))
        return
    if len(sys.argv)>2:
        fitorder = int(sys.argv[2])
    ports = [0,1,4,5,12,14,15]
    '''
    USING retardation results independently
    for p in ports:
        log2tofs.update( {'port_%i'%p:[]} )
        energies.update( {'port_%i'%p:[]} )
    '''
    with h5py.File(sys.argv[1],'r+') as f:
        for ret in f.keys():
            print(ret)
            for p in ports:
                #log2tofs['port_%i'%p] += [float(v) for v in f[ret]['port_%i'%p]['log2tofs'][()]]
                #energies['port_%i'%p] += [float(v) for v in f[ret]['port_%i'%p]['energies'][()]]
                log2tofs = [float(v) for v in f[ret]['port_%i'%p]['log2tofs'][()]]
                log2energies = [np.log2(float(v)) for v in f[ret]['port_%i'%p]['energies'][()]]
                x0,theta=fitpoly(log2tofs,log2energies,order=fitorder)
                f[ret]['port_%i'%p].attrs['x0'] = x0
                f[ret]['port_%i'%p].attrs['theta'] = theta
                xout = np.arange(-1,1,.05)
                X = mypoly(xout,order=fitorder)
                Y = X.dot(f[ret]['port_%i'%p].attrs['theta'])
                np.savetxt('tmp.fit.%s.port_%i.dat'%(ret,p),np.column_stack((xout+f[ret]['port_%i'%p].attrs['x0'],Y)),fmt='%.6f')
                np.savetxt('tmp.vals.%s.port_%i.dat'%(ret,p),np.column_stack((log2tofs,log2energies)),fmt='%.6f')
    print('finished')
    return

if __name__ == '__main__':
    main()

