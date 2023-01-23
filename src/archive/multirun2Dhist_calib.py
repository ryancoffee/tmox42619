#!/usr/bin/python3

import h5py
import numpy as np
import sys
from scipy.interpolate import interp1d 
from scipy.sparse import coo_matrix,csr_matrix
import re
from utils import mypoly,fitpoly,getcentroid

def logt2e(lt,coeffs,x0):
    if lt<1.:
        return 0.
    x = np.array([(np.log(float(lt))-x0)**i for i in range(4)])
    return np.exp(x.dot(coeffs))


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

def loadh5calib(fname,portnum):
    with h5py.File(fname,'r') as f:
        for n,key in enumerate(f.keys()):
            if n>0:
                continue
            theta = f[key]['port_%i'%portnum].attrs['theta']
            x0 = f[key]['port_%i'%portnum].attrs['x0']
    return theta,x0

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
    t0 = {0:109830,1:100451,2:99810,4:97180,5:99071,12:98561,13:98657,14:100331,15:106956,16:97330}
    #slopethresh = {'port_0':50,'port_1':50,'port_2':30,'port_4':15,'port_5':50,'port_12':50,'port_13':50,'port_14':50,'port_15':50,'port_16':30}
    slopethresh = {'port_0':90,'port_1':90,'port_2':90,'port_4':90,'port_5':90,'port_12':90,'port_13':90,'port_14':90,'port_15':90,'port_16':90}
    vlsoffsetdict = {82:141,83:141,84:141,85:141,86:141,87:0,88:0,89:0,90:0,93:0,94:0,95:0}
    vlsoffsetdict.update({7:0,8:50,9:100,10:150,11:175,12:225,13:250,14:275,15:300,16:325,17:125})
    vlsoffsetdict.update({21:0,22:0,23:25,26:0,27:0,28:0,29:0,30:0,31:0,36:0,39:0,40:0})
    vlsoffsetdict.update({59:0,61:50,62:100,63:150,66:200,67:250,68:300,69:325,70:350,72:375,73:400,74:425})
    vlsmean = np.array((None,),dtype=float)
    tofv = [] # for vls spectrometer index
    toft = [] # for time-of-flight index
    en = [] # for energy index
    tofd = [] # for counting
    ncor = 0

    ########## setting up for energy calibration ##############
    theta,x0 = loadh5calib(calibfname,portnum)

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
        print('port_%i\tt0 = %i'%(portnum,t0))
        
        tofaddresses = data['port_%i'%portnum]['addresses']
        tofnedges= data['port_%i'%portnum]['nedges']

        for i in range(len(tofnedges)):
            if tofnedges[i] > 0:
                tmpt = (tof[ int(tofaddresses[i]):int(tofaddresses[i]+tofnedges[i]) ])/2. -t0 # OK, this is getting silling with all the scaling
                tmpslope = tofslope[ int(tofaddresses[i]):int(tofaddresses[i]+tofnedges[i]) ]
                if len(tmpt)==len(tmpslope):
                    toflist = []
                    slopelist = []
                    for j in range(len(tmpt)):
                        if tmpslope[j]>slopethresh['port_%i'%portnum]: 
                            if tmpt[j] < 2**18:
                                toflist += [tmpt[j]]
                                slopelist += [np.int32(tmpslope[j])]
                    if i%1000==0: 
                        print(i,toflist,slopelist)
                    if len(vlsinds[i])>150:
                        j = getcentroid(vlsinds[i],vlsspec[i])
                        x = []
                        for t in toflist:
                            t *= 2.
                            if t>0 and t<2**16:
                                toft += [t]
                                x += [np.log2(t)]
                        #x = np.log2(toft) - x0
                                tofv += [j+vlsoffset]
                                tofd += [1]
                                ncor += 1
                        X = mypoly(x-x0,order=2)
                        Y = X.dot(np.array(theta))
                        ens = np.power(2.,3+Y)#.astype(np.uint16)
                        en += list(ens)

    return toft,tofv,tofd,en,vlsmean,ncor


def main():
    if len(sys.argv)<3:
        print('syntax: main <calibfilename> <h5file list>')
        return
    calibfname = sys.argv[1]
    fnames = list(sys.argv[2:])
    fullname = sys.argv[2]
    datadir = './'
    filefront = 'hits.h5'

    m = re.search('^(.+)/h5files/(.+).h5',fnames[0])
    if m:
        datadir = m.group(1)
    else:
        print('syntax: main <calibfilename> <h5file list>')
        return


    for portnum in (0,12):#,4,5,12,13,14,15):
        toft,tofv,tofd,en,vlsmean,ncor = processmultitofs_log(fnames,calibfname,portnum)
        print(len(toft),len(tofv),len(tofd))
        tofhist_csr = coo_matrix((tofd,(toft,tofv)),shape=((2**18,vlsmean.shape[0]+512)),dtype=np.int16).tocsr()

        #enhist_csr = coo_matrix((tofd,(en,tofv)),shape=((2**10,vlsmean.shape[0]+512)),dtype=np.uint16).tocsr()

        tofhistname = '%s/ascii/Argon_port_%i.tofhist.dat'%(datadir,portnum)
        bins = np.arange(2**16)
        h,b = np.histogram(toft,bins)
        np.savetxt(tofhistname,np.column_stack((b[:-1],h)),fmt='%i')

        #histname = '%s/ascii/Argon_port_%i.hist.dat'%(datadir,portnum)
        enhistname = '%s/ascii/Argon_port_%i.enhist.dat'%(datadir,portnum)
        h,b = np.histogram(en,np.arange(2**12))
        np.savetxt(enhistname,np.column_stack((b[:-1],h)),fmt='%i')
        #tofname = '%s/ascii/Argon_port_%i.tof.dat'%(datadir,portnum)
        #dokname = '%s/ascii/Argon_port_%i.dok.dat'%(datadir,portnum)
        #vlsname = '%s/ascii/Argon_port_%i.vls.dat'%(datadir,portnum)

        #np.savetxt(vlsname,vlsmean,fmt='%.2f')
        #h= tofhist_csr.toarray()
        #np.savetxt(histname,h,fmt='%i')
        #np.savetxt(tofname,np.sum(h,axis=1),fmt='%i')
        #e= enhist_csr.toarray()
        #np.savetxt(enhistname,e,fmt='%i')
        #f = open(dokname,'w')
        #tofhist_dok = tofhist_csr.todok()
        #for key in tofhist_dok.keys():
        #    print('%i\t%i\t%i'%(key[0],key[1],tofhist_dok[key] ),file=f)
        #f.close()
        print('processed port %i'%portnum)
    return


if __name__ == '__main__':
    main()
