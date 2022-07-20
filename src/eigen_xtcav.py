#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.2.5/bin/python3

import numpy as np
import sys
import re
import h5py
from scipy.fftpack import dct
import multiprocessing as mp

'''
data = np.array([images[i,:,:].flatten() for i in range(images.shape[0])])
avgim = np.mean(data,axis=0)
residual = np.array([data[i,:]-avgim for i in range(data.shape[0])])
COVMAT = np.cov(residual.T)
vals,vecs = np.linalg.eig(COVMAT) # this peice takes a long time.

Then remove mean, then remove a handfull of np.dot(eigvecs,data) and plot the >0 pixels
'''

def processBlock(params):
    with h5py.File('%s/%s.h5'%(params['inpath'],params['infile']),'r') as inh5:
        (nims,xdim,ydim) = inh5['xtcav']['images'][()].shape 
        nims //= params['nblocks'] 
        print('newshape = (%i,%i,%i)'%(nims,xdim,ydim))
        datablock = np.column_stack([inh5['xtcav']['images'][()][params['tid']+params['nblocks']*j,:,:].flatten() for j in range(nims) ] )
        print('datablock %i shape\t%s'%(params['tid'],datablock.shape))
    outfile = '%s.eigenxtcav.tid%i.h5'%(params['infile'],params['tid'])
    params['outfile'] = outfile
    covMat = np.cov(datablock)
    print('covMat.shape',covMat.shape)
    eigMat,eigVals = np.linalg.eig(covMat)
    print(eigMat.shape)
    return params

def main():
    if len(sys.argv) < 3:
        print('Syntax: eigen_xtcav.py <h5 full path and name> <outputpath>')
        return
    fullinfile = sys.argv[1]
    m = re.search('^(.+)/(.+)\.h5$',fullinfile)
    if not m:
        print(fullinfile)
        return
    paramslist = [{} for i in range(10)]
    for i,p in enumerate(paramslist):
        p['inpath'] = str(m.group(1))
        p['infile'] = str(m.group(2))
        p['outpath'] = str(sys.argv[2])
        p['nblocks'] = np.uint8(100)
        p['tid'] = np.uint8(i)

    pool = mp.Pool(len(paramslist))
    returnlist = pool.map(processBlock,paramslist)
    return

if __name__ == '__main__':
    main()
