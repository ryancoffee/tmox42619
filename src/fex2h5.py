#!/sdf/group/lcls/ds/ana/sw/conda2/manage/bin/psconda.sh

import psana
import numpy as np
import sys
import re
import h5py
import os
import socket

from Ports import *
from Ebeam import *
from Vls import *
from Gmd import *
from Config import Config
from utils import *




def main(nshots:int,expname:str,runnums:list,scratchdir:str):
  
    runnums = [int(run) for run in sys.argv[3:]]
    outnames = ['%s/hits.%s.run_%03i.h5'%(scratchdir,expname,r) for r in runnums]

    _=[print('starting analysis exp %s for run %i'%(expname,int(r))) for r in runnums]
    cfgname = '%s/%s.%s.configs.h5'%(scratchdir,expname,os.environ.get('USER'))
    configs = Config()
    params = configs.writeconfigs(cfgname).getparams()
    chans = params['chans']
    t0s = params['t0s']
    logicthresh = params['logicthresh']
    offsets = params['offsets']

    nr_expand = params['expand']


    '''
    spect = [Vls(params['vlsthresh']) for r in runnums]
    _ = [s.setwin(params['vlswin'][0],params['vlswin'][1]) for s in spect]
    #_ = [s.setthresh(params['vlsthresh']) for s in spect]
    ebunch = [Ebeam() for r in runnums]
    gmd = [Gmd() for r in runnums]
    _ = [e.setoffset(params['l3offset']) for e in ebunch]
    '''

    port = [{} for r in runnums] 
    inflate = params['inflate']

    for r in range(len(runnums)):
        for key in chans.keys():
            port[r][key] = Port(key,chans[key],t0=t0s[key],logicthresh=logicthresh[key],inflate=inflate,expand=nr_expand,nrolloff=2**6)

    ds = [psana.DataSource(exp=expname,run=r) for r in runnums]

    #for run in ds.runs():
    runs = [next(ds[r].runs()) for r in range(len(runnums))]
        #np.savetxt('%s/waveforms.%s.%i.%i.dat'%(scratchdir,expname,runnum,key),wv[key],fmt='%i',header=headstring)
    for r in range(len(runnums)):
        print(runs[r].detnames)
    return
"""
    runhsd=True
    runvls=True
    runebeam=True
    runxtcav=False
    rungmd=True
    hsds = []
    vlss = []
    ebeams = []
    xtcavs = []
    xgmds = []
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
        if rungmd and 'xgmd' in runs[r].detnames:
            xgmds += [runs[r].Detector('xgmd')]

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
        gmdEvents =[]

        init = True 
        vsize = 0

        print('chans: ',chans)
        for evt in runs[r].events():
            if eventnum > nshots:
                break

            completeEvent = [True]

            if runxtcav and bool(np.prod(completeEvent)):
                try:
                    if type(xtcav.raw.value(evt)) == None:
                        print(eventnum,'skip per problem with XTCAV')
                        completeEvent += [False]
                        continue
                    else:
                        img = np.copy(xtcav.raw.value(evt)).astype(np.int16)
                        mf = np.argmax(np.histogram(img,np.arange(2**8))[0])
                        img -= mf
                        imgcrop,x0,y0 = xtcav_crop(img,win=(512,256))
                        xtcavImages += [imgcrop]
                        xtcavX0s += [x0]
                        xtcavY0s += [y0]
                        completeEvent += [True]
                except:
                    print(eventnum,'skipping xtcav, skip per failed try:')
                    completeEvent += [False]
                    continue


## test vlswv
            vlswv = None
            if runvls and bool(np.prod(completeEvent)):
                ''' VLS specific section, do this first to slice only good shots '''
                if type(vlss[r]) == None:
                    print(eventnum,'skip per problem with VLS')
                    completeEvent += [False]
                    continue
                vlswv = np.squeeze(vlss[r].raw.value(evt))
                completeEvent += [spect[r].test(vlswv)]


## test thisl3
            thisl3 = None
            if runebeam and bool(np.prod(completeEvent)):
                ''' Ebeam specific section '''
                if type(ebeams[r]) == None:
                    print(eventnum,'ebeam is None')
                    completeEvent += [False]
                    continue
                thisl3 = ebeams[r].raw.ebeamL3Energy(evt)
                completeEvent += [ebunch[r].test(thisl3)]


## test thisgmde
            thisgmde = None
            if rungmd and bool(np.prod(completeEvent)):
                if type(xgmds[r]) == None:
                    print(eventnum,'gmd is None')
                    completeEvent += [False]
                    continue
                thisgmde = xgmds[r].raw.energy(evt)
                completeEvent += [gmd[r].test(thisgmde)]


## test hsds
            if runhsd and bool(np.prod(completeEvent)):
                if type(hsds[r]) == None:
                    print(eventnum,'hsds is None')
                    completeEvent += [False]
                for key in chans.keys(): # here key means 'port number'
                    completeEvent += [port[r][key].test(hsds[r].raw.waveforms(evt)[ chans[key] ][0])]


## process VLS
            if runvls and bool(np.prod(completeEvent)):
                spect[r].process(vlswv)

## process ebeam
            if runebeam and bool(np.prod(completeEvent)):
                ebunch[r].process(thisl3)

## process gmd
            if rungmd and bool(np.prod(completeEvent)):
                gmd[r].process(thisgmde)


## process hsds
            if runhsd and bool(np.prod(completeEvent)):
                ''' HSD-Abaco section '''
                for key in chans.keys(): # here key means 'port number'
                    s = np.array(hsds[r].raw.waveforms(evt)[ chans[key] ][0] , dtype=np.int16) # presumably 12 bits unsigned input, cast as int16_t since will immediately in-place subtract baseline
                    port[r][key].process(s)

## redundant events vec
            if bool(np.prod(completeEvent)):
                if runebeam:
                    ebeamEvents += [eventnum]
                if runvls:
                    vlsEvents += [eventnum]
                if runxtcav:
                    xtcavEvents += [eventnum]
                if rungmd:
                    gmdEvents += [eventnum]
                if runhsd:
                    hsdEvents += [eventnum]



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


            if init:
                init = False
                ebunch[r].set_initState(False)
                spect[r].set_initState(False)
                gmd[r].set_initState(False)
                for key in chans.keys():
                    port[r][key].set_initState(False)

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
    print("Hello, I'm done now.  Have a most exceelent day!")
    return

"""




if __name__ == '__main__':
    if len(sys.argv)>3:
        nshots = int(sys.argv[1])
        expname = sys.argv[2]
        runnums = [int(r) for r in list(sys.argv[3:])]
        scratchdir = '/sdf/data/lcls/ds/tmo/%s/scratch/%s/h5files/%s/'%(expname,os.environ.get('USER'),socket.gethostname())
        if not os.path.exists(scratchdir):
            os.makedirs(scratchdir)
        main(nshots,expname,runnums,scratchdir)
    else:
        print('Please give me a number of shots (-1 for all), experiment, a list of run numbers, and output directory is to expt scratch)')
