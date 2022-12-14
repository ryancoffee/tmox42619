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
        params['vlsthresh'] = f.attrs['vlsthresh']
        params['l3offset'] = f.attrs['l3offset']
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
    #scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_vernier/h5files'
    scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_debug/h5files'

    if len(sys.argv)<3:
        print('syntax: ./hits2h5.py <nshots> <expname> <list of run numbers>')
    expname = sys.argv[2]
    nshots = int(sys.argv[1])
    runnums = [int(run) for run in sys.argv[3:]]
    outnames = ['%s/hits.%s.run_%03i.h5'%(scratchdir,expname,r) for r in runnums]

    _=[print('starting analysis exp %s for run %i'%(expname,int(r))) for r in runnums]
    cfgname = '%s/%s.configs.h5'%(scratchdir,expname)
    params = fillconfigs(cfgname)
    chans = params['chans']
    t0s = params['t0s']
    logicthresh = params['logicthresh']

    nr_expand = params['expand']


    spect = [Vls(params['vlsthresh']) for r in runnums]
    #_ = [s.setthresh(params['vlsthresh']) for s in spect]
    ebunch = [Ebeam() for r in runnums]
    _ = [e.setoffset(params['l3offset']) for e in ebunch]
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
    runhsd=True
    runvls=True
    runebeam=True
    runxtcav=False
    rungmd=False
    hsds = []
    vlss = []
    ebeams = []
    xtcavs = []
    gmds = []
    for r in range(len(runnums)):
        eventnum = 0
        print('processing run %i'%runnums[r])
        if runhsd and 'hsd' in runs[r].detnames:
            hsds += [runs[r].Detector('hsd')]
        if runvls and 'andor' in runs[r].detnames:
            vlss += [runs[r].Detector('andor')]
        if runebeam and 'ebeam' in runs[r].detnames:
            ebeams += [runs[r].Detector('ebeam')]
        if runxtcav and 'xtcav' in runs[r].detnames:
            xtcavs += [runs[r].Detector('xtcav')]
        if rungmd and 'gms' in runs[r].detnames:
            gmds += [runs[r].Detector('gmd')]

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
        ebeamEvents = []

        init = True 
        vsize = 0

        print('chans: ',chans)
        for evt in runs[r].events():
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
                if type(vlss[r]) == None:
                    print(eventnum,'skip per problem with VLS')
                    continue
                vlswv = np.squeeze(vlss[r].raw.value(evt))
                if spect[r].process(vlswv):
                    vlsEvents += [eventnum]
                else:
                    #print(eventnum,'skip per low vls')
                    continue

            if runebeam:
                ''' Ebeam specific section '''
                thisl3 = ebeams[r].raw.ebeamL3Energy(evt)
                if type(thisl3)==None:
                    print(eventnum,'l3 is None')
                    continue
                if ebunch[r].process(thisl3):
                    ebeamEvents += [eventnum]
                else:
                    print(eventnum,'skipping for l3')
                    continue

            if runhsd:
    
                ''' HSD-Abaco section '''
                goodevents = 0
                for key in chans.keys(): # here key means 'port number'
                    try:
                        s = np.array(hsds[r].raw.waveforms(evt)[ chans[key] ][0] , dtype=np.int16) 
                        if port[r][key].process(s):
                            goodevents += 1
                        else:
                            print(eventnum, 'hsd process == False for %s'%key)
                            continue
                    except:
                        print(eventnum, 'failed hsd for some reason')
                        continue
                hsdEvents += [eventnum]

                if init:
                    init = False
                    ebunch[r].set_initState(False)
                    spect[r].set_initState(False)
                    for key in chans.keys():
                        port[r][key].set_initState(False)


                if eventnum<2:
                        print('ports = %s'%([k for k in chans.keys()]))
                if eventnum<100:
                    if eventnum%10==0: 
                        print('working event %i,\tnedges = %s'%(eventnum,[port[r][k].getnedges() for k in chans.keys()] ))
                elif eventnum<1000:
                    if eventnum%100==0: 
                        print('working event %i,\tnedges = %s'%(eventnum,[port[r][k].getnedges() for k in chans.keys()] ))
                else:
                    if eventnum%1000==0: 
                        print('working event %i,\tnedges = %s'%(eventnum,[port[r][k].getnedges() for k in chans.keys()] ))


            if eventnum > 1 and eventnum <1000 and eventnum%100==0:
                with h5py.File(outnames[r],'w') as f:
                    print('writing to %s'%outnames[r])
                    if runhsd:
                        Port.update_h5(f,port[r],hsdEvents,chans)
                    if runvls:
                        Vls.update_h5(f,spect[r],vlsEvents)
                    if runebeam:
                        Ebeam.update_h5(f,ebunch[r],ebeamEvents)
                    if rungmd:
                        Gmd.update_h5(f,gmd[r],gmdEvents)

            elif eventnum>900 and eventnum%1000==0:
                with h5py.File(outnames[r],'w') as f:
                    print('writing to %s'%outnames[r])
                    if runhsd:
                        Port.update_h5(f,port[r],hsdEvents,chans)
                    if runvls:
                        Vls.update_h5(f,spect[r],vlsEvents)
                    if runebeam:
                        Ebeam.update_h5(f,ebunch[r],ebeamEvents)
                    if rungmd:
                        Gmd.update_h5(f,gmd[r],gmdEvents)

            eventnum += 1
 

        with h5py.File(outnames[r],'w') as f:
            print('writing to %s'%outnames[r])
            if runhsd:
                Port.update_h5(f,port[r],hsdEvents,chans)
            if runvls:
                Vls.update_h5(f,spect[r],vlsEvents)
            if runebeam:
                Ebeam.update_h5(f,ebunch[r],ebeamEvents)
            if rungmd:
                Gmd.update_h5(f,gmd[r],gmdEvents)

        print('Finished with run %i'%runnums[r])
    print("Hello, I'm done now!")
    return

if __name__ == '__main__':
    main()
