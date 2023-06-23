#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import re
from Calibration import CalibData

def main():
    calib = CalibData()
    calib.fit()
    with h5py.File(sys.argv[-1]) as f:
        for w in calib.wins.keys():
            for k in calib.portkeys:
                if w in calib.theta[k].keys():
                    q=np.arange(len(f[k]['qbins'][()]) -1)
                    s=1./(f[k]['qbins'][()][1:] - f[k]['qbins'][()][:-1])
                    dim = len(calib.theta[k][w])
                    sz = len(s)
                    X = np.ones((dim,sz),dtype=float)
                    for i in range(1,dim):
                        X[i,:] = np.power(q,i)
                    Y = np.dot(calib.theta[k][w],X)
                    plt.plot(Y,s,label='%s, %s'%(k,w))
                    plt.show()
    return

if __name__ == '__main__':
    if len(sys.argv) > 1:
        main()
    else:
        print('syntax: main <quantized.h5 file path and name>')
