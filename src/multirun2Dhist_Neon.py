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
'''
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
'''

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
    t0 = {'port_0':27460,'port_1':25114,'port_2':24981,'port_4':24295,'port_5':24768,'port_12':24645,'port_13':24669,'port_14':25087,'port_15':26742,'port_16':24507}
    slopethresh = {'port_0':500,'port_1':500,'port_2':300,'port_4':150,'port_5':500,'port_12':500,'port_13':500,'port_14':500,'port_15':500,'port_16':300}
    vlsoffsetdict = {82:141,83:141,84:141,85:141,86:141,87:0,88:0,89:0,90:0,93:0,94:0,95:0}
    vlsmean = np.array((None,),dtype=float)
    tofv = [] # for vls spectrometer index
    toflogt = [] # for time-of-flight index
    tofd = [] # for counting
    ncor = 0
    for fname in fnames:
        data = {}
        attrs = {}
        data,attrs = loadh5(fname)
        m = re.search('run(\d+).h5$',fname)
        vlsoffset = 0
        if m:
            vlsoffset = vlsoffsetdict[int(m.group(1))]

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
        vlsspec = []
        vlsinds = []
        for i,p in enumerate(peaks):
            vlsinds.append([i for i in range( max(0,1024+p+vlswindow[0]) , min(1900,1024+p+vlswindow[1]) )])
            tmp = data['vls']['data'][i,vlsinds[-1]]   - np.max(data['vls']['data'][i,1024:1044]) ## restricting to high indices for sake of capturing only 3rd order (4th order near pixel 100)
            tmp *= (tmp>0)
            vlscsum.append( np.cumsum( tmp ) )
            vlsspec.append( tmp )

        tof = data['port_%i'%portnum]['tofs']
        tofslope = data['port_%i'%portnum]['slopes']
        #t0 = attrs['port_%i'%portnum]['t0']
        tofaddresses = data['port_%i'%portnum]['addresses']
        tofnedges= data['port_%i'%portnum]['nedges']

        for i in range(len(tofnedges)):
            if tofnedges[i] > 0:
                tmpt = tof[ int(tofaddresses[i]):int(tofaddresses[i]+tofnedges[i]) ] - t0['port_%i'%portnum]
                tmpslope = tofslope[ int(tofaddresses[i]):int(tofaddresses[i]+tofnedges[i]) ]
                if len(tmpt)==len(tmpslope):
                    toflist = []
                    slopelist = []
                    for j in range(len(tmpt)):
                        if tmpslope[j]>slopethresh['port_%i'%portnum]: ## HERE HERE HERE HERE setting by hand
                            toflist += [np.int16(tmpt[j])]
                            slopelist += [np.int32(tmpslope[j])]
                    if i%1000==0: 
                        print(i,t0['port_%i'%portnum],toflist,slopelist)
                    if len(vlsinds[i])>150:
                        #for j in samplePDF(vlsinds[i],vlscsum[i],len(toflist)): ### choose a fresh random set of 16 for each hit
                        for j in getcentroid(vlsinds[i],vlsspec[i]): 
                            for t in toflist:
                                if t>1 and t<(2**15-1):
                                    toflogt += [np.log2(t)/16.0*2**12] ### SCALING HERE ###
                                    #toflogt += [t]
                                    tofv += [j+vlsoffset]
                                    tofd += [1]
                                    ncor += 1
    return toflogt,tofv,tofd,vlsmean,ncor


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


    toft,tofv,tofd,vlsmean,ncor = processmultitofs_log(fnames,portnum)
    tofhist_csr = coo_matrix((tofd,(toft,tofv)),shape=((2**14,vlsmean.shape[0])),dtype=np.uint16).tocsr()

    histname = '%s/ascii/port_%i.hist.dat'%(datadir,portnum)
    tofname = '%s/ascii/port_%i.tof.dat'%(datadir,portnum)
    #outername = '%s/ascii/port_%i.outer.dat'%(datadir,portnum)
    #diffname = '%s/ascii/port_%i.diff.dat'%(datadir,portnum)
    dokname = '%s/ascii/port_%i.dok.dat'%(datadir,portnum)
    vlsname = '%s/ascii/port_%i.vls.dat'%(datadir,portnum)

    np.savetxt(vlsname,vlsmean,fmt='%.2f')
    h= tofhist_csr.toarray()[2000:4000,:] # this index selection is a function of the scaling at ### SCALING HERE ###
    np.savetxt(histname,h,fmt='%i')
    np.savetxt(tofname,np.sum(h,axis=1),fmt='%i')
    f = open(dokname,'w')
    tofhist_dok = tofhist_csr.todok()
    for key in tofhist_dok.keys():
        print('%i\t%i\t%i'%(key[0],key[1],tofhist_dok[key] ),file=f)
    f.close()
    #h = tofhist_csr.toarray().astype(float)
    #a = np.sum(h,axis=0)
    #b = np.sum(h,axis=1) 
    #c = np.outer(b,a)/ncor
    #np.savetxt(outername,c,fmt='%.6f')
    #np.savetxt(diffname,h-c,fmt='%.4f')
    #print(np.max(h-c))
    return


if __name__ == '__main__':
    main()
