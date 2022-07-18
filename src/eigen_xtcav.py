#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.2.5/bin/python3

import numpy as np
import sys
import re
import h5py
from scipy.fftpack import dct

'''
data = np.array([images[i,:,:].flatten() for i in range(images.shape[0])])
avgim = np.mean(data,axis=0)
residual = np.array([data[i,:]-avgim for i in range(data.shape[0])])
COVMAT = np.cov(residual.T)
vals,vecs = np.linalg.eig(COVMAT) # this peice takes a long time.

Then remove mean, then remove a handfull of np.dot(eigvecs,data) and plot the >0 pixels
'''

def main():
    if len(sys.argv) < 3:
        print('Syntax: eigen_xtcav.py <h5 full path and name> <outputpath>')
        return
    fullinfile = sys.argv[1]
    outpath = sys.argv[2]
    m = re.search('^(.+)/(.+)\.h5$',fullinfile)
    if not m:
        print(fullinfile)
        return
    print('%s\t%s'%(m.group(1),m.group(2)))
    inpath = m.group(1)
    infile = m.group(2)
    with h5py.File('%s/%s'%(inpath,infile)) as inh5:
        datablock = []
        for i in range(10):
            datablock += [inh5['xtcav']['images'][()][i::10,:,:] ]
    return

if __name__ == '__main__':
    main()
