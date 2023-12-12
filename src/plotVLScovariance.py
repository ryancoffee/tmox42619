#!/sdf/group/lcls/ds/ana/sw/conda2/inst/envs/h5-1.0.1/bin/python3

import h5py
import numpy as np
import sys
from math import log2
import matplotlib.pyplot as plt
from Vls import Vls,shouldSplit,getCentroid,getSplit

def main():
    docovmat:bool = False
    for fname in sys.argv[1:]:
        vlsdata = {}
        issplit:bool = False
        with h5py.File(fname,'r') as f:
            print(f.keys())
            vlsdata.update({'full':f['vls']['data'][()]})
            print(vlsdata['full'].shape)
            meanspec = np.sum(vlsdata['full'],axis=0)
            rightshift = int(log2(np.max(meanspec))) - int(log2(np.max(vlsdata['full']))) -1
            meanspec >>= rightshift
            plt.plot(meanspec)
            plt.show()
            if shouldSplit(meanspec):
                vlsdata.update({'order3':vlsdata['full'][:,:getSplit(meanspec)]})
                vlsdata.update({'order2':vlsdata['full'][:,getSplit(meanspec):]})
                issplit = True
            if issplit:
                splitkeys = [k for k in vlsdata.keys() if k != 'full']
                if docovmat:
                    _ = [plt.plot(np.sum(vlsdata[k],axis=0),label=k) for k in splitkeys]
                    plt.show()
                    covmat = np.cov(vlsdata['full'].T)[getSplit(meanspec):,:getSplit(meanspec)]
                    plt.imshow(covmat)
                    plt.show()
                _=[plt.plot(100*float(i)+vlsdata['order2'][i,:]) for i in range(0,100,10)]
                _=[plt.plot(np.arange(vlsdata['order3'].shape[1])+37,100*float(i)+1.5*vlsdata['order3'][i,:]) for i in range(0,100,10)]
                plt.xlim(100,200)
                plt.show()
                '''
                plt.imshow(covmat[:vlsdata[splitkeys[0]].shape[0],vlsdata[splitkeys[0]].shape[0]:])
                plt.show()
                '''
    '''
    vls2=vlsdata[:,((1<<10)+(1<<9)):((1<<11)-(1<<8))]
    covmat=np.cov(vls2.T,vls3.T)
    ridge = [np.argmax(covmat[i,:]) for i in range(covmat.shape[0])]
    plt.plot(ridge,'.');plt.show()
    '''
    return

if __name__ == '__main__':
    if len(sys.argv)>1:
        main()
    else:
        print('I need an .h5 file to process.')
