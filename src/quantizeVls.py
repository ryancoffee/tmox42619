#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3
import numpy as np
import sys
import h5py
import re
import matplotlib.pyplot as plt
from Quantizers import Quantizer

def main():
    if len(sys.argv)<3:
        print('syntax: quantizeHits.py <ntofbins> <nvlsbins> <fnames>')
        return

    fnames = sys.argv[3:]
    ntofbins = np.uint32(sys.argv[1])
    nvlsbins = np.uint32(sys.argv[2])
    vlsoffset = 1000
    tofs = {} 
    addresses = {} 
    nedges = {} 
    quants = {}
    hist = {}
    vlscenters = []
    portkeys = []
    for fname in fnames:
        with h5py.File(fname,'r') as f:
            portkeys = [k for k in f.keys() if (re.search('port',k) and not re.search('_16',k) and not re.search('_2',k))]
            if len(quants.keys())==0:
                for k in portkeys:
                    quants[k] = Quantizer(style='nonuniform',nbins=ntofbins)
                    tofs[k] = list(f[k]['tofs'][()])
                    addresses[k] = list(f[k]['addresses'][()])
                    nedges[k] = list(f[k]['nedges'][()])
                    hist[k] = np.zeros((nvlsbins,ntofbins),dtype=int)
                vlscenters = list(f['vls']['centroids'][()])
            else:
                for k in portkeys:
                    offsetTofs = len(tofs[k])
                    addresses[k] += [int(offsetTofs+a) for a in f[k]['addresses'][()]]
                    nedges[k] += list(f[k]['nedges'][()])
                    tofs[k] += list(f[k]['tofs'][()])
                vlscenters += list(f['vls']['centroids'][()])
    _=[print('%s\t%i\t%i'%(k,len(tofs[k]),addresses[k][-1]+nedges[k][-1])) for k in portkeys]
    for k in portkeys:
        quants[k].setbins(data=tofs[k])

    vlsbins = [np.uint32(max(0,min(int(v-vlsoffset),nvlsbins-1))) for v in vlscenters]
    for shot,vlsbin in enumerate(vlsbins):
        for k in portkeys:
            a = addresses[k][shot]
            n = nedges[k][shot]
            try:
                hist[k][vlsbin,:] += quants[k].histogram(tofs[k][a:a+n])
            except:
                print(vlsbin,a,n)

    for k in portkeys:
        plt.imshow(hist[k],origin='lower')
        plt.show()
        #outname = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_vernier/ascii/test_%s_hist.dat'%(k)

    return

if __name__=='__main__':
    main()

