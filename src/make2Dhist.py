#!/usr/bin/python3

import h5py
import numpy as np
import sys


def main():
    if len(sys.argv)<2:
        print('syntax: make2Dhist.py <fnames>')
        return
    fnames = list(sys.argv[1:])
    data = {}

    for fname in fnames:
        with h5py.File('./data_fs/h5files/%s'%(fname),'r') as f:
            #print(list(f.keys()))
            for key in f.keys():
                data[key] = np.squeeze(f[key][()])
        # here now the data is stored in system memory and we close the h5 file
        vlscsum = np.cumsum(data['vls'],axis=1)
        np.savetxt('./data_fs/ascii/tmpvls.dat',vlscsum[:10,:].T,fmt='%i')
    return


if __name__ == '__main__':
    main()
