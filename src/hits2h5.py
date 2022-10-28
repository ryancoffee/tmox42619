#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.2.5/bin/python3


import psana
import numpy as np
import sys
import re
import h5py
from scipy.fftpack import dct,dst

from Ports import *
from Ebeam import *
from Vls import *
from utils import *


def fillconfigs(cfgname):
    params = {'chans':{},'t0s':{},'logicthresh':{}}
    with h5py.File(cfgname,'r') as f:
        params['inflate'] = f.attrs['inflate']
        params['expand'] = f.attrs['expand']
        for p in f.keys():
            m = re.search('^\w+_(\d+)$',p)
            if m:
                k = int(m.group(1))
                params['chans'][k] = f[p].attrs['hsd']
                params['t0s'][k] = f[p].attrs['t0']
                params['logicthresh'][k] = f[p].attrs['logicthresh']
    return params

def xtcav_crop(inimg,win=(256,256)):
    # hard coded factor of 2 scale down
    xprof = np.mean(inimg,axis=0)
    yprof = np.mean(inimg,axis=1)
    y0 = np.argmax(xprof)
    x0 = np.argmax(yprof)
    resimg = (np.roll(inimg,(-x0+win[0]//2,-y0+win[1]//2),axis=(0,1)))[:win[0],:win[1]]
    tmp= np.column_stack((resimg,np.flip(resimg,axis=1)))
    outimg=np.row_stack((tmp,np.flip(tmp,axis=0)))
    W = dct(dct(outimg,axis=1,type=2),axis=0,type=2)
    xenv = np.zeros(W.shape[0])
    yenv = np.zeros(W.shape[1])
    xenv[:win[0]//2] = 0.5*(1+np.cos(np.arange(win[0]//2)*np.pi/(win[0]/2)))
    yenv[:win[1]//2] = 0.5*(1+np.cos(np.arange(win[1]//2)*np.pi/(win[1]/2)))
    for i in range(W.shape[1]//2):
        W[:,i] *= xenv
    for i in range(W.shape[0]//2):
        W[i,:] *= yenv
    W *= 4.0/np.product(W.shape)
    out = dct( dct(W[:win[0]//2,:win[1]//2],type=3,axis=0),type=3,axis=1)[:win[0]//4,:win[1]//4]
    return out,x0//2,y0//2
    #return dct(dct(W,axis=2,type=3),axis=1,type=3),x0,y0
    #print(x0,y0)
    #return inimg[:win[0],:win[1]],x0,y0
    

def main():
        ############################################
        ###### Change this to your output dir ######
        ############################################
    #scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output/h5files'
    #scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_2022/h5files'
    scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_multicolorhack/h5files'

    if len(sys.argv)<3:
        print('syntax: ./hits2h5.py <nshots> <expname> <list of run numbers>')
    expname = sys.argv[1]
    nshots = int(sys.argv[2])
    runnums = [int(run) for run in sys.argv[3:]]

    print('starting analysis exp %s for run %i'%(expname,int(runnum)))
    cfgname = '%s/%s.hsdconfig.h5'%(scratchdir,expname)
    params = fillconfigs(cfgname)
    chans = params['chans']
    t0s = params['t0s']
    logicthresh = params['logicthresh']

    nr_expand = params['expand']

    spect = [Vls() for r in runnums]
    ebunch = [Ebeam() for r in runnums]
    port = [{} for r in runnums] 
    scale = int(1) # to better fill 16 bit int
    inflate = params['inflate']
    for key in logicthresh.keys():
        logicthresh[key] *= scale # inflating by factor of 4 since we are also scaling the waveforms by 4 in vertical to fill bit depth.

    for r in range(len(runnums)):
        for key in chans.keys():
            port[r][key] = Port(key,chans[key],t0=t0s[key],logicthresh=logicthresh[key],inflate=inflate,expand=nr_expand,scale=scale,nrolloff=2**6)

    ds = [psana.DataSource(exp=expname,run=r) for r in runnums]

    #for run in ds.runs():
    runs = [next(ds[r].runs()) for r in range(len(runnums))]
        #np.savetxt('%s/waveforms.%s.%i.%i.dat'%(scratchdir,expname,runnum,key),wv[key],fmt='%i',header=headstring)
    for r in range(len(runnums)):
        print(runs[r].detnames)
    eventnum = 0
    runhsd=True
    runvls=True
    runebeam=True
    runxtcav=False
    hsds = []
    vlss = []
    ebeams = []
    xtcavs = []
    for r in range(len(runnums)):
        if runhsd and 'hsd' in runs[r].detnames:
            hsds += [runs[r].Detector('hsd')]
        if runvls and 'andor' in runs[r].detnames:
            vlss += [runs[r].Detector('andor')]
        if runebeam and 'ebeam' in run.detnames:
            ebeams += [run[r].Detector('ebeam')]
        if runxtcav and 'xtcav' in run.detnames:
            xtcavs += [runs[r].Detector('xtcav')]

####### HERE HERE HERE HERE ###########
####### finish working with lists of runs ###

        wv = {}
        wv_logic = {}
        v = [] # vls data matrix
        vc = [] # vls centroids vector
        vs = [] # vls sum is I think not used, maybe for normalization or used to be for integration and PDF sampling
        l3 = [] # e-beam l3 (linac 3) in GeV.
        xtcavImages = []
        xtcavX0s = []
        xtcavY0s = []
        xtcavEvents = []

        vlsEvents = []
        hsdEvents = []


        init = True 
        vsize = 0

        print('chans: ',chans)
        for evt in run.events():
            if eventnum > nshots:
                break

            if runxtcav:
                ## HERE HERE HERE HERE ##
                ## change this to xtcav.process() style... build xtcav object like the others
                try:
                    if type(xtcav.raw.value(evt)) == None:
                        print(eventnum,'skip per problem with XTCAV')
                        continue
                    img = np.copy(xtcav.raw.value(evt)).astype(np.int16)
                    mf = np.argmax(np.histogram(img,np.arange(2**8))[0])
                    img -= mf
                    imgcrop,x0,y0 = xtcav_crop(img,win=(512,256))
                    xtcavImages += [imgcrop]
                    xtcavX0s += [x0]
                    xtcavY0s += [y0]
                    xtcavEvents += [eventnum]
                except:
                    print(eventnum,'skipping xtcav, skip per failed try:')
                    continue



            if runvls:
                ''' VLS specific section, do this first to slice only good shots '''
                try:
                    if type(vls) == None:
                        print(eventnum,'skip per problem with VLS')
                        continue
                    vlswv = np.squeeze(vls.raw.value(evt))
                    vlswv = vlswv-int(np.mean(vlswv[1900:])) # this subtracts baseline
                    if np.max(vlswv)<1:  # too little amount of xrays
                        print(eventnum,'skip per negative vls')
                        #eventnum += 1
                        continue
                    spect.process(vlswv)
                    vlsEvents += [eventnum]
                    #spect.print_v()

                except:
                    print(eventnum,'skip per vls')
                    continue

            if runebeam:
                ''' Ebeam specific section '''
                try:
                    thisl3 = ebeam.raw.ebeamL3Energy(evt)
                    thisl3 += 0.5
                    ebunch.process(thisl3)
                except:
                    print(eventnum,'skipping ebeam, skip per l3')
                    continue

            if runhsd:
    
                ''' HSD-Abaco section '''
                for key in chans.keys(): # here key means 'port number'
                    #try:
                    s = np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=np.int16) 
                    port[key].process(s)

                    if init:
                        init = False
                        ebunch.set_initState(False)
                        spect.set_initState(False)
                        for key in chans.keys():
                            port[key].set_initState(False)
                    #except:
                     #   print(eventnum, 'failed hsd for some reason')
                      #  continue

                hsdEvents += [eventnum]

                if eventnum<10:
                        print('ports = %s'%([k for k in chans.keys()]))
                if eventnum<100:
                    if eventnum%10<2: 
                        print('working event %i,\tnedges = %s'%(eventnum,[port[k].getnedges() for k in chans.keys()] ))
                elif eventnum<1000:
                    if eventnum%100<2: 
                        print('working event %i,\tnedges = %s'%(eventnum,[port[k].getnedges() for k in chans.keys()] ))
                else:
                    if eventnum%1000<2: 
                        print('working event %i,\tnedges = %s'%(eventnum,[port[k].getnedges() for k in chans.keys()] ))
                eventnum += 1

        f = h5py.File('%s/hits.%s.run%i.h5'%(scratchdir,expname,runnum),'w') 
                # use f.create_group('port_%i'%i,portnum)
        #_ = [print(key,chans[key]) for key in chans.keys()]
        if runhsd:
            for key in chans.keys(): # remember key == port number
                g = f.create_group('port_%i'%(key))
                g.create_dataset('tofs',data=port[key].tofs,dtype=np.int32) 
                g.create_dataset('slopes',data=port[key].slopes,dtype=np.int32) 
                g.create_dataset('addresses',data=port[key].addresses,dtype=np.uint64)
                g.create_dataset('nedges',data=port[key].nedges,dtype=np.uint32)
                wvgrp = g.create_group('waves')
                lggrp = g.create_group('logics')
                for k in port[key].waves.keys():
                    wvgrp.create_dataset(k,data=port[key].waves[k],dtype=np.int16)
                    lggrp.create_dataset(k,data=port[key].logics[k],dtype=np.int32)
                g.attrs.create('inflate',data=port[key].inflate,dtype=np.uint8)
                g.attrs.create('expand',data=port[key].expand,dtype=np.uint8)
                g.attrs.create('t0',data=port[key].t0,dtype=float)
                g.attrs.create('logicthresh',data=port[key].logicthresh,dtype=np.int32)
                g.attrs.create('hsd',data=port[key].hsd,dtype=np.uint8)
                g.attrs.create('size',data=port[key].sz*port[key].inflate,dtype=int) ### need to also multiply by expand #### HERE HERE HERE HERE
                g.create_dataset('events',data=hsdEvents)
        if runxtcav:
            grpxtcav = f.create_group('xtcav')
            grpxtcav.create_dataset('images',data=xtcavImages,dtype=np.int16)
            grpxtcav.create_dataset('x0s',data=xtcavX0s,dtype=np.float16)
            grpxtcav.create_dataset('y0s',data=xtcavY0s,dtype=np.float16)
            grpxtcav.create_dataset('xtcavEvents',data=xtcavEvents,dtype=np.float16)
            grpxtcav.create_dataset('events',data=xtcavEvents)

        if runvls:
            grpvls = f.create_group('vls')
            grpvls.create_dataset('data',data=spect.v,dtype=np.int16)
            grpvls.create_dataset('centroids',data=spect.vc,dtype=np.int16)
            grpvls.create_dataset('sum',data=spect.vs,dtype=np.uint64)
            grpvls.attrs.create('size',data=spect.vsize,dtype=np.int32)
            grpvls.create_dataset('events',data=vlsEvents)
        if runebeam:
            grpebeam = f.create_group('ebeam')
            grpebeam.create_dataset('l3energy',data=ebunch.l3,dtype=np.uint16)
        f.close()

    print("Hello, I'm done now!")
    return

if __name__ == '__main__':
    main()
