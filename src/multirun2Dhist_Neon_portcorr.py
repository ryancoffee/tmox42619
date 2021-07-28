#!/usr/bin/python3

import h5py
import numpy as np
import sys
from scipy.interpolate import interp1d 
from scipy.sparse import coo_matrix,csr_matrix
import re

def getcentroid(inds,spec):
    x = inds
    y = spec 
    num = np.sum(x*y)
    denom = np.sum(y)
    if (denom>0):
        return [int(num/denom)]
    return [0]

def samplePDF(inds,csum,nsamples = 16):
    x = csum
    y = inds
    cs = {}
    for j in range(len(x)):
        cs.update({x[j]:y[j]})
    if len(cs.keys())<50:
        return []
    xx = [v for v in cs.keys()]
    yy = [v for v in cs.values()]
    f = interp1d(xx,yy)
    xnew = np.random.random((nsamples,))*(xx[-2]-xx[1]) + xx[1]
    return list(f(xnew).astype(np.uint16))

def loadh5(fname):
    data = {}
    attrs = {}
    m = re.search('^(.+)/h5files/(.+).h5',fname)
    if m:
        fullname = m.group(0)
        datadir = m.group(1)
        filefront = m.group(2)
        print('processing:\t%s'%fullname)
    else:
        return
    with h5py.File('%s'%(fullname),'r') as f:
        for key in f.keys():
            if re.search('port',key):
                data[key] = {}
                attrs[key]= {}
                for k in f[key].keys():
                    data[key].update({k:np.squeeze(f[key][k][()])})
                for k in f[key].attrs.keys():
                    attrs[key].update({k:np.squeeze(f[key].attrs[k])})
            if re.search('vls',key):
                data[key] = {}
                for k in f[key].keys():
                    data[key].update({k: np.squeeze(f[key][k][()])})
    return data,attrs

def processmultitofs_log(fnames,portnums):
    t0 = {'port_0':26900,'port_1':24625,'port_2':25000,'port_4':23800,'port_5':24200,'port_12':24250,'port_13':24300,'port_14':24600,'port_15':26200,'port_16':25000}
    vlsmean = np.array((None,),dtype=float)
    toflogt = [] # for time-of-flight index
    tofd = [] # for counting
    ncor = 0
    tof = {}
    tofslope = {}
    tofaddresses = {}
    tofnedges = {}
    toflogt = {}
    for fname in fnames:
        data = {}
        attrs = {}
        data,attrs = loadh5(fname)
        # here now the data is stored in system memory and we close the h5 file
        print(data.keys())
        for key in data.keys():
            print(key,data[key].keys())


        for p in portnums:
            tof[p] = data['port_%i'%p]['tofs']
            tofslope[p] = data['port_%i'%p]['slopes']
            tofaddresses[p] = data['port_%i'%p]['addresses']
            tofnedges[p] = data['port_%i'%p]['nedges']

        cor = coo_matrix(([],([],[])),shape = (2**14,2**14),dtype=np.int16)
        for i in range(len(tofnedges[portnums[0]])):
            for p in portnums:
                toflogt[p] = []
                if tofnedges[p][i] > 0:
                    tmpt = tof[p][ int(tofaddresses[p][i]):int(tofaddresses[p][i]+tofnedges[p][i]) ] - t0['port_%i'%p]
                    tmpslope = tofslope[p][ int(tofaddresses[p][i]):int(tofaddresses[p][i]+tofnedges[p][i]) ]
                    if len(tmpt)==len(tmpslope):
                        toflist = []
                        slopelist = []
                        for j in range(len(tmpt)):
                            if tmpslope[j]>500: ## HERE HERE HERE HERE setting by hand
                                toflist += [np.int16(tmpt[j])]
                                slopelist += [np.int32(tmpslope[j])]
                                for t in toflist:
                                    if t>1 and t<(2**14-1):
                                        toflogt[p] += [t]

            t1s = []
            t2s = []
            for t1 in toflogt[portnums[0]]:
                for t2 in toflogt[portnums[1]]:
                    t1s += [t1]
                    t2s += [t2]
            ds = [1]*len(t1s)
            cor += coo_matrix((ds,(t1s,t2s)),shape=(2**14,2**14),dtype=np.int16)
            ncor += 1

    return cor.tocsr(),ncor


def main():
    if len(sys.argv)<2:
        print('syntax: make2Dhist.py <fnames>')
        return
    fnames = list(sys.argv[1:])
    datadir = './'
    filefront = 'hits.h5'

    m = re.search('^(.+)/h5files/(.+).h5',fnames[0])
    if m:
        datadir = m.group(1)
    else:
        print('syntax: main h5file list')
        return


    portnums = [12,15]
    cor,ncor = processmultitofs_log(fnames,portnums)

    histname = '%s/ascii/ports_%i_%i.hist.dat'%(datadir,portnums[0],portnums[1])
    outername = '%s/ascii/ports_%i_%i.outer.dat'%(datadir,portnums[0],portnums[1])
    diffname = '%s/ascii/ports_%i_%i.diff.dat'%(datadir,portnums[0],portnums[1])

    np.savetxt(histname,cor.toarray(),fmt='%i')
    h = cor.toarray().astype(float)
    a = np.sum(h,axis=0)
    b = np.sum(h,axis=1) 
    c = np.outer(b,a)/ncor
    np.savetxt(outername,c,fmt='%.6f')
    np.savetxt(diffname,h-c,fmt='%.4f')
    print(np.max(h-c))
    return


if __name__ == '__main__':
    main()
