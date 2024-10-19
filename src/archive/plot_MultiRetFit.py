#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import re
from Calibration import CalibData

def main():
    wins = {'n1s':[180,200]
        ,'nam':[360,375]
        ,'oam':[480,510]
        ,'IV':[553,567]
        }
    ywins = {'n1s':[0,0.75]
        ,'nam':[0,0.9]
        ,'oam':[0,1.2]
        ,'IV':[0,0.6]
        }
    scale = {'n1s':[2,1,1]
        ,'nam':[1.5,2,1,2]
        ,'oam':[.5,.5,0.75,1.5,1,0.5]
        ,'IV':[1,1.5,1,1,1,1]
        }
    ports = {'n1s':['port_12','port_13','port_5']
        ,'nam':['port_12','port_13','port_4','port_5']
        ,'oam':['port_0','port_1','port_12','port_14','port_15','port_4']
        ,'IV':['port_1','port_12','port_14','port_15','port_4','port_5']
        }
    calib = CalibData()
    calib.fit()
    with h5py.File(sys.argv[-1]) as f:
        for w in wins.keys():
        #for w in calib.wins.keys():
            legendvec = []
            plt.figure(figsize=(4,8))
            for j,k in enumerate(ports[w]):
                if w in calib.theta[k].keys():
                    q=np.arange(len(f[k]['qbins'][()]) -1)
                    s=1./(f[k]['qbins'][()][1:] - f[k]['qbins'][()][:-1])
                    dim = len(calib.theta[k][w])
                    sz = len(s)
                    X = np.ones((dim,sz),dtype=float)
                    for i in range(1,dim):
                        X[i,:] = np.power(q,i)
                    Y = np.dot(calib.theta[k][w],X)
                    plt.plot(Y,scale[w][j]*s)
                    legendvec += ['%s'%(k)]
            plt.xlim(wins[w][0],wins[w][1])
            plt.ylim(ywins[w][0],ywins[w][1])
            plt.xlabel('electron energy [eV]')
            plt.ylabel('signal [arb.]')
            plt.legend(legendvec)
            #plt.title(w)
            plt.savefig('./figures/runs132-134_%s.png'%(w))
    return

if __name__ == '__main__':
    if len(sys.argv) > 1:
        main()
    else:
        print('syntax: main <quantized.h5 file path and name>')
