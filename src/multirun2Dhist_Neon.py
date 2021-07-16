#!/usr/bin/python3

import h5py
import numpy as np
import sys
from scipy.interpolate import interp1d 
from scipy.sparse import coo_matrix,csr_matrix
import re

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

'''
def samplePDFmat(inds,csum,nsamples = 16):
    b = []
    x = csum
    y = inds
    for i in range(x.shape[0]):
        cs = {}
        for j in range(x.shape[1]):
            cs.update({x[i,j]:y[j]})
        if len(cs.keys())<50:
            b.append(np.zeros((nsamples,),dtype=np.uint16))
        f = interp1d(xx,yy)
        xnew = np.random.random((nsamples,))*(xx[-2]-xx[1]) + xx[1]
        b.append(f(xnew).astype(np.uint16))
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
'''

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

def processmultitofs_log(fnames,portnum):
    t0 = {'port_0':26900,'port_1':24625,'port_2':25000,'port_4':23800,'port_5':24200,'port_12':24250,'port_13':24300,'port_14':24600,'port_15':26200,'port_16':25000}
    vlsmean = np.array((None,),dtype=float)
    tofv = [] # for vls spectrometer index
    toflogt = [] # for time-of-flight index
    tofd = [] # for counting
    for fname in fnames:
        data = {}
        attrs = {}
        data,attrs = loadh5(fname)
        # here now the data is stored in system memory and we close the h5 file
        print(data.keys())
        for key in data.keys():
            print(key,data[key].keys())

        if vlsmean.shape[0]<data['vls']['data'].shape[1]:
            vlsmean = np.mean(data['vls']['data'],axis=0)
        else:
            vlsmean += np.mean(data['vls']['data'],axis=0)
        vlswindow = (-100,100)
        peaks = np.argmax(data['vls']['data'][:,1024:],axis=1)
        vlscsum = []
        vlsinds = []
        for i,p in enumerate(peaks):
            '''
            if len(vlsinds)==0:
                vlsinds = [i for i in range( max(0,p-vlswindow[0]) , min(1900,p+vlswindow[1]) )]
                vlscsum = np.cumsum( (data['vls'][i,vlsinds[-1]] - np.max(data['vls'][i,:20])) )
            else:
            '''
            vlsinds.append([i for i in range( max(0,1024+p+vlswindow[0]) , min(1900,1024+p+vlswindow[1]) )])
            tmp = data['vls']['data'][i,vlsinds[-1]]   - np.max(data['vls']['data'][i,1024:1044]) ## restricting to high indices for sake of capturing only 3rd order (4th order near pixel 100)
            tmp *= (tmp>0)
            vlscsum.append( np.cumsum( tmp ) )

        tof = data['port_%i'%portnum]['tofs']
        #t0 = attrs['port_%i'%portnum]['t0']
        tofaddresses = data['port_%i'%portnum]['addresses']
        tofnedges= data['port_%i'%portnum]['nedges']


        for i in range(len(tofnedges)):
            if tofnedges[i] > 0:
                tmp = tof[ int(tofaddresses[i]):int(tofaddresses[i]+tofnedges[i]) ] - t0['port_%i'%portnum]
                toflist = ((np.array([np.log2(t) for t in tmp if t>5])-5)*2**8).astype(np.int16)
                if i%1000==0: print(i,t0['port_%i'%portnum],toflist)
                for t in toflist:
                    if len(vlsinds[i])>150:
                        for j in samplePDF(vlsinds[i],vlscsum[i],16): ### choose a fresh random set of 16 for each hit
                            if t>1 and t<(2**16-1):
                                toflogt += [t]
                                tofv += [j]
                                tofd += [1]
    return toflogt,tofv,tofd,vlsmean


def processmultitofs(fnames,portnum):
    vlsmean = np.array((None,),dtype=float)
    tofv = [] # for vls spectrometer index
    toft = [] # for time-of-flight index
    tofd = [] # for counting
    for fname in fnames:
        data = {}
        data = loadh5(fname)
        # here now the data is stored in system memory and we close the h5 file
        print(data.keys())
        for key in data.keys():
            print(key,data[key].keys())

        if vlsmean.shape[0]<data['vls']['data'].shape[1]:
            vlsmean = np.mean(data['vls']['data'],axis=0)
        else:
            vlsmean += np.mean(data['vls']['data'],axis=0)
        vlswindow = (-100,100)
        peaks = np.argmax(data['vls']['data'][:,1024:],axis=1)
        vlscsum = []
        vlsinds = []
        for i,p in enumerate(peaks):
            vlsinds.append([i for i in range( max(0,1024+p+vlswindow[0]) , min(1900,1024+p+vlswindow[1]) )])
            tmp = data['vls']['data'][i,vlsinds[-1]]   - np.max(data['vls']['data'][i,1024:1044])
            tmp *= (tmp>0)
            vlscsum.append( np.cumsum( tmp ) )

        tof = data['port_%i'%portnum]['tofs']
        tofaddresses = data['port_%i'%portnum]['addresses']
        tofnedges= data['port_%i'%portnum]['nedges']


        for i in range(len(tofnedges)):
            if tofnedges[i] > 0:
                toflist = tof[ int(tofaddresses[i]):int(tofaddresses[i]+tofnedges[i]) ].astype(np.uint32)
                if i%1000==0: print(i,toflist)
                for t in toflist:
                    if len(vlsinds[i])>150:
                        for j in samplePDF(vlsinds[i],vlscsum[i],16): ### choose a fresh random set of 16 for each hit
                            if t<(2**16-1):
                                toft += [t]
                                tofv += [j]
                                tofd += [1]
    return toft,tofv,tofd,vlsmean

def main():
    if len(sys.argv)<3:
        print('syntax: make2Dhist.py <portnum> <fnames>')
        return
    portnum = np.int8(sys.argv[1])
    fnames = list(sys.argv[2:])
    fullname = sys.argv[2]
    datadir = './'
    filefront = 'hits.h5'

    m = re.search('^(.+)/h5files/(.+).h5',fnames[0])
    if m:
        datadir = m.group(1)
    else:
        print('syntax: main h5file list')
        return


    toft,tofv,tofd,vlsmean = processmultitofs_log(fnames,portnum)
    #tofhist_dok = coo_matrix((tofd,(toft,tofv)),shape=((2**12,vlsmean.shape[0])),dtype=np.uint16).tocsr().todok()
    tofhist_csr = coo_matrix((tofd,(toft,tofv)),shape=((2**14,vlsmean.shape[0])),dtype=np.uint16).tocsr()

    histname = '%s/ascii/port_%i.hist.dat'%(datadir,portnum)
    dokname = '%s/ascii/port_%i.dok.dat'%(datadir,portnum)
    vlsname = '%s/ascii/port_%i.vls.dat'%(datadir,portnum)

    np.savetxt(vlsname,vlsmean,fmt='%.2f')
    np.savetxt(histname,tofhist_csr.toarray(),fmt='%i')
    f = open(dokname,'w')
    tofhist_dok = tofhist_csr.todok()
    for key in tofhist_dok.keys():
        print('%i\t%i\t%i'%(key[0],key[1],tofhist_dok[key] ),file=f)
    f.close()
    return


if __name__ == '__main__':
    main()
