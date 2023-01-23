#!/usr/bin/python3

import h5py
import numpy as np
import sys
from scipy.interpolate import interp1d 
from scipy.sparse import coo_matrix,csr_matrix
import re

def logt2e(lt,coeffs,x0):
    if lt<1.:
        return 0.
    x = np.array([(np.log(float(lt))-x0)**i for i in range(4)])
    return np.exp(x.dot(coeffs))

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

def loadh5calib(fname):
    thetas = {}
    x0s = {}
    with h5py.File(fname,'r') as f:
        for key in f.keys():
            for p in f[key].keys():
                if re.search('port_\d+',p):
                    thetas.update( {p:f[key][p].attrs['theta']} )
                    x0s.update( {p:f[key][p].attrs['x0']} )
    return thetas,x0s

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

def processmultitofs_log(fnames,calibfname,portnum):
    ### eventually, use t0['Neon'] = {...} and then have a <sample> parameter that is either 'Neon', 'Argon', or 'NNO'
    # from Neon runs 80s-90s 
    # t0 = {'port_0':27460,'port_1':25114,'port_2':24981,'port_4':24295,'port_5':24768,'port_12':24645,'port_13':24669,'port_14':25087,'port_15':26742,'port_16':24507}
    ### Argon runs 7-17 t0s
    ### now we need to account for hte Newton-Raphson expansion as well
    ### t0 = {'port_0':27913,'port_1':25570,'port_2':24900,'port_4':24752,'port_5':25225,'port_12':25100,'port_13':25120,'port_14':25540,'port_15':27198,'port_16':24000}
    #t0 = {'port_0':27913,'port_1':25570,'port_2':24900,'port_4':24752,'port_5':25225,'port_12':25100,'port_13':25120,'port_14':25540,'port_15':27198,'port_16':24000}
    #t0 = {0:109830,1:100451,2:99810,4:97180,5:99071,12:98561,13:98657,14:100331,15:106956,16:97330}
    slopethresh = {'port_0':500,'port_1':500,'port_2':300,'port_4':150,'port_5':500,'port_12':500,'port_13':500,'port_14':500,'port_15':500,'port_16':300}
    vlsoffsetdict = {82:141,83:141,84:141,85:141,86:141,87:0,88:0,89:0,90:0,93:0,94:0,95:0}
    vlsoffsetdict.update({7:0,8:50,9:100,10:150,11:175,12:225,13:250,14:275,15:300,16:325,17:125})
    vlsoffsetdict.update({21:0,22:0,23:25,26:0,27:0,28:0,29:0,30:0,31:0,36:0,39:0,40:0})
    vlsoffsetdict.update({59:0,61:50,62:100,63:150,66:200,67:250,68:300,69:325,70:350,72:375,73:400,74:425})
    vlsmean = np.array((None,),dtype=float)
    tofv = [] # for vls spectrometer index
    toflogt = [] # for time-of-flight index
    en = []
    tofd = [] # for counting
    ncor = 0

    ########## setting up for energy calibration ##############
    lt2e_coeffs = [5. ,-4., -0.5, 1.5] 
    x0 = 6.7
    with h5py.File(calibfname,'r') as lt2ecalib:
        lt2e_coeffs = lt2ecalib['port_%i'%portnum]['lt2e_coeffs'][()]
        x0 = lt2ecalib['x0'][()]
    ##########################################################

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
        #peaks = np.argmax(data['vls']['data'][:,1024:],axis=1) ## for 700eV photons
        peaks = np.argmax(data['vls']['data'][:,256:],axis=1)
        vlscsum = []
        vlsspec = []
        vlsinds = []
        for i,p in enumerate(peaks):
            vlsinds.append([i for i in range( max(0,256+p+vlswindow[0]) , min(1900,256+p+vlswindow[1]) )]) ## for below 600eV or 700eV photons when the 2nd and 3rd order don't fit in window
            #vlsinds.append([i for i in range( max(0,1024+p+vlswindow[0]) , min(1900,1024+p+vlswindow[1]) )]) ## for 700eV photons
            tmp = data['vls']['data'][i,vlsinds[-1]] - np.max(data['vls']['data'][i,128:256]) ## restricting to high indices for sake of capturing only 3rd order (4th order near pixel 100)
            #tmp = data['vls']['data'][i,vlsinds[-1]] - np.max(data['vls']['data'][i,1024:1044]) ## only for 700eV photons.. restricting to high indices for sake of capturing only 3rd order (4th order near pixel 100)
            tmp *= (tmp>0)
            vlscsum.append( np.cumsum( tmp ) )
            vlsspec.append( tmp )

        tof = data['port_%i'%portnum]['tofs']
        tofslope = data['port_%i'%portnum]['slopes']
        t0 = attrs['port_%i'%portnum]['t0']
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
                        if tmpslope[j]>slopethresh['port_%i'%portnum]: 
                            toflist += [np.int16(tmpt[j])]
                            slopelist += [np.int32(tmpslope[j])]
                    if i%1000==0: 
                        print(i,t0['port_%i'%portnum],toflist,slopelist)
                    if len(vlsinds[i])>150:
                        #for j in samplePDF(vlsinds[i],vlscsum[i],len(toflist)): ### choose a fresh random set of 16 for each hit
                        for j in getcentroid(vlsinds[i],vlsspec[i]): 
                            for t in toflist:
                                if t>1 and t<(2**15-1):
                                    tmp = np.log2(t)/16.0*2**14 ### SCALING HERE ###
                                    toflogt += [tmp] 
                                    if tmp>2001:
                                        tmpen = logt2e(tmp-2000.,lt2e_coeffs,x0) ### SCALING HERE ###
                                        if tmpen<2**14:
                                            en += [int(tmpen)]
                                        else:
                                            en += [int(0)]
                                    else:
                                        en += [int(0)]
                                    #toflogt += [t]
                                    tofv += [j+vlsoffset]
                                    tofd += [1]
                                    ncor += 1
    return toflogt,tofv,tofd,en,vlsmean,ncor


def main():
    if len(sys.argv)<3:
        print('syntax: main <portnum> <h5file list>')
        return
    portnum = np.int8(sys.argv[1])
    fnames = list(sys.argv[2:])
    fullname = sys.argv[2]
    datadir = './'
    calibfname = 'calib/Neon_lt2e.h5'
    filefront = 'hits.h5'

    m = re.search('^(.+)/h5files/(.+).h5',fnames[0])
    if m:
        datadir = m.group(1)
        calibfname = '%s/calib/Neon_lt2e.h5'%(m.group(1))
    else:
        print('syntax: main <portnum> <h5file list>')
        return


    toft,tofv,tofd,en,vlsmean,ncor = processmultitofs_log(fnames,calibfname,portnum)
    '''
    ## HERE HERE HERE the vls shape is getting larger than it should be because I'm adding the vls offset for retardation... I really shouldn't do that.
    print(vlsmean.shape[0])
    print(np.max(tofv))
    print(np.max(toft))
    '''
    tofhist_csr = coo_matrix((tofd,(toft,tofv)),shape=((2**14,vlsmean.shape[0]+512)),dtype=np.uint16).tocsr()
    enhist_csr = coo_matrix((tofd,(en,tofv)),shape=((2**14,vlsmean.shape[0]+512)),dtype=np.uint16).tocsr()

    histname = '%s/ascii/Argon_port_%i.hist.dat'%(datadir,portnum)
    enhistname = '%s/ascii/Argon_port_%i.enhist_high.dat'%(datadir,portnum)
    tofname = '%s/ascii/Argon_port_%i.tof.dat'%(datadir,portnum)
    dokname = '%s/ascii/Argon_port_%i.dok.dat'%(datadir,portnum)
    vlsname = '%s/ascii/Argon_port_%i.vls.dat'%(datadir,portnum)

    np.savetxt(vlsname,vlsmean,fmt='%.2f')
    h= tofhist_csr.toarray()[8000:,256:] # this index selection is a function of the scaling at ### SCALING HERE ###
    #h= tofhist_csr.toarray()[8000:,1024:] # this index selection is a function of the scaling at ### SCALING HERE ###
    np.savetxt(histname,h,fmt='%i')
    np.savetxt(tofname,np.sum(h,axis=1),fmt='%i')
    e= enhist_csr.toarray()[:4000,1024:]
    np.savetxt(enhistname,e,fmt='%i')
    f = open(dokname,'w')
    tofhist_dok = tofhist_csr.todok()
    for key in tofhist_dok.keys():
        print('%i\t%i\t%i'%(key[0],key[1],tofhist_dok[key] ),file=f)
    f.close()
    return


if __name__ == '__main__':
    main()
