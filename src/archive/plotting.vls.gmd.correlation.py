#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3

import sys
import h5py
import matplotlib.pyplot as plt
import numpy as np

from Quantizers import Quantizer
from utils import inlims


def main(nbins,fname):
    gmdquant = Quantizer(style='nonuniform',nbins=nbins)
    vlsmeans = np.zeros((nbins,2048),dtype=float)
    nspectra = np.zeros((nbins),dtype=int)
    vlsdata = np.zeros((1,),dtype=int)
    binIDs = np.zeros((1,),dtype=int)
    with h5py.File(fname,'r') as f:
        gmdquant.setbins(data=f['gmd']['gmdenergy'][()])
        binIDs = gmdquant.getbins(f['gmd']['gmdenergy'][()])
        vlsdata = np.copy(f['vls']['data'][()])

    for s,b in enumerate(binIDs):
        vlsmeans[b,:] += vlsdata[s,:]
        nspectra[b] += 1

    for i,s in enumerate(vlsmeans):
        if i>0:
            plt.plot(s/float(nspectra[i]),label='%i uJ'%(gmdquant.bincenters()[i]))
    plt.ylabel('intensity')
    plt.xlabel('vls index')
    plt.xlim((100,200))
    plt.ylim((0,2500))
    plt.legend()
    plt.savefig('./figs/VlsVsGmd.png')
    plt.show()
    return

if __name__ == '__main__':
    if (len(sys.argv)<3):
        print('syntax: plotting.vls.gmd.correlation.py <nbins_gmd> <h5filename>')
    else:
        main(int(sys.argv[1]),sys.argv[2])

