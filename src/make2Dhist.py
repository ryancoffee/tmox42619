#!/usr/bin/python3

import h5py
import numpy as np
import sys
from scipy.interpolate import interp1d 
from scipy.sparse import coo_matrix,csr_matrix
import re

def samplePDFmat(csum,nsamples = 10):
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

def main():
    if len(sys.argv)<2:
        print('syntax: make2Dhist.py <fnames>')
        return
    fnames = list(sys.argv[1:])
    data = {}
    fullname = sys.argv[1]
    datadir = './'
    filefront = 'hits.h5'

    for fname in fnames:
        #with h5py.File('./data_fs/h5files/%s'%(fname),'r') as f:
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
        h = fill2dhist(vlscsum)

        spinds = samplePDFmat(vlscsum,32) # repeate this for every tof individually
        portstring = 'port_14'
        tof = data['%s_tofs'%portstring]
        tofaddresses = data['%s_addresses'%portstring]
        tofnedges= data['%s_nedges'%portstring]
        print(tof.shape)
        print('tofs:\t\t'+'\t\t'.join((tof[:10]).astype(str)))
        print('adds:\t\t'+'\t\t'.join(tofaddresses[:10].astype(str)))
        print('edges:\t\t'+'\t\t'.join(tofnedges[:10].astype(str)))

        shotlist = []
        toflist = []

        tofv = [] # for vls spectrometer index
        toft = [] # for time-of-flight index
        tofd = [] # for counting
        for i in range(len(tofnedges)):
            if tofnedges[i] > 0:
                toflist = tof[ int(tofaddresses[i]):int(tofaddresses[i]+tofnedges[i]) ].astype(np.uint32)
                shotlist.append(np.uint16(i))
                if i%10==0: print(shotlist[-1],toflist)
                for j in spinds[shotlist[-1]]: ### change this to choose a fresh random set of 16 for each hit
                    #tofhist[toflist[i],j] += 1
                    for t in toflist:
                        if t<2**16:
                            toft += [t]
                            tofv += [j]
                            tofd += [1]

        #tofhist = np.zeros((2**16,vlsmean.shape[0]),dtype=np.int8)
        print(np.max(toft),np.max(tofv))

        tofhist_dok = coo_matrix((tofd,(toft,tofv)),shape=((2**16,vlsmean.shape[0])),dtype=np.uint16).tocsr().todok()




        np.savetxt('%s/ascii/%s.origvls.dat'%(datadir,filefront),data['vls'][:10,:].T,fmt='%i')
        np.savetxt('%s/ascii/%s.tmpvls.dat'%(datadir,filefront),vlscsum[:10,:].T,fmt='%i')
        #np.savetxt('%s/ascii/%s.hist.dat'%(datadir,filefront),tofhist_csr,fmt='%i')
        histname = '%s/ascii/%s.%s.hist.dat'%(datadir,filefront,portstring)
        f = open(histname,'w')
        for key in tofhist_dok.keys():
            print('%i\t%i\t%i'%(key[0],key[1],tofhist_dok[key] ),file=f)
        f.close()
    return


if __name__ == '__main__':
    main()
