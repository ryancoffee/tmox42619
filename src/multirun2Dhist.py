#!/usr/bin/python3

import h5py
import numpy as np
import sys
from scipy.interpolate import interp1d 
from scipy.sparse import coo_matrix,csr_matrix
import re

def samplePDF(csum,nsamples = 10):
    b = []
    x = csum
    y = np.arange(x.shape[1],dtype=np.uint8) 
    for i in range(x.shape[0]):
        cs = {}
        for j in range(x.shape[1]):
            cs.update({x[i,j]:y[j]})
        xx = [v for v in cs.keys()]
        yy = [v for v in cs.values()]
        f = interp1d(xx,yy)
        xnew = np.random.random((nsamples,))*(xx[-2]-xx[1]) + xx[1]
        b.append(list(f(xnew).astype(np.uint16)))
    return b

def fill2dhist(csum):
    x = csum
    y = np.arange(x.shape[1],dtype=np.uint8) 
    h = np.zeros(x.shape,dtype=np.int8)
    for i in range(x.shape[0]):
        cs = {}
        for j in range(x.shape[1]):
            cs.update({x[i,j]:y[j]})
        xx = [v for v in cs.keys()]
        yy = [v for v in cs.values()]
        f = interp1d(xx,yy)
        xnew = np.random.random((100,))*(xx[-2]-xx[1]) + xx[1]
        ynew = f(xnew).astype(np.uint16)
        for k in ynew:
            h[i,k] += 1
    return h

def processmultitofs(fnames,toft,tofv,tofd,portstring):
    for fname in fnames:
        m = re.search('^(.+)/h5files/(.+).h5',fname)
        if m:
            fullname = m.group(0)
            datadir = m.group(1)
            filefront = m.group(2)
            print('processing:\t%s'%fullname)
        else:
            return
        with h5py.File('%s'%(fullname),'r') as f:
            print(list(f.keys()))
            for key in f.keys():
                data[key] = np.squeeze(f[key][()])
        # here now the data is stored in system memory and we close the h5 file
        vlsmean = np.mean(data['vls'],axis=0)
        vlscsum = np.cumsum(data['vls'],axis=1)

        spinds = samplePDF(vlscsum,256) # repeate this for every tof individually

        tof = data['%s_tofs'%portstring]
        tofaddresses = data['%s_addresses'%portstring]
        tofnedges= data['%s_nedges'%portstring]

        shotlist = []
        toflist = []

        for i in range(len(tofnedges)):
            if tofnedges[i] > 0:
                toflist = tof[ int(tofaddresses[i]):int(tofaddresses[i]+tofnedges[i]) ].astype(np.uint32)
                shotlist.append(np.uint16(i))
                if i%100==0: print(shotlist[-1],toflist)
                for j in spinds[shotlist[-1]]: ### change this to choose a fresh random set of 16 for each hit
                    for t in toflist:
                        if t<(2**16-1):
                            toft += [t]
                            tofv += [j]
                            tofd += [1]
    return toft,tofv,tofd

def main():
    if len(sys.argv)<2:
        print('syntax: make2Dhist.py <fnames>')
        return
    fnames = list(sys.argv[1:])
    data = {}
    fullname = sys.argv[1]
    datadir = './'
    filefront = 'hits.h5'

    m = re.search('^(.+)/h5files/(.+).h5',fname)
    if m:
        datadir = m.group(1)
    else:
        print('syntax: main h5file list')
        return

    portstring = 'port_14'
    tofv = [] # for vls spectrometer index
    toft = [] # for time-of-flight index
    tofd = [] # for counting

    toft,tofv,tofd = processmultitofs(fnames,toft,tofv,tofd,portstring)
    tofhist_dok = coo_matrix((tofd,(toft,tofv)),shape=((2**16,vlsmean.shape[0])),dtype=np.uint16).tocsr().todok()

    histname = '%s/ascii/%s.hist.dat'%(datadir,portstring)
    f = open(histname,'w')
    for key in tofhist_dok.keys():
        print('%i\t%i\t%i'%(key[0],key[1],tofhist_dok[key] ),file=f)
    f.close()
    return


if __name__ == '__main__':
    main()
