#!/sdf/group/lcls/ds/ana/sw/conda2/manage/bin/psconda.sh

import psana
import numpy as np
import sys
import re
import h5py
import os
import socket

from typing import Type,List

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
    inflate = params['inflate']

    _= [print(k) for k in chans.keys()]


    '''
    spect = [Vls(params['vlsthresh']) for r in runnums]
    _ = [s.setwin(params['vlswin'][0],params['vlswin'][1]) for s in spect]
    #_ = [s.setthresh(params['vlsthresh']) for s in spect]
    ebunch = [Ebeam() for r in runnums]
    gmd = [Gmd() for r in runnums]
    _ = [e.setoffset(params['l3offset']) for e in ebunch]
    '''
    runhsd=True
    runtiming=True
    runfzp=False

    runvls=False
    runebeam=False
    runxtcav=False
    rungmd=False
    hsds = []
    chankeys = {}
    timings = []
    fzps = []
    vlss = []
    ebeams = []
    xtcavs = []
    xgmds = []

    '''
    Setting up the datasource
    '''
    dsdict = {}
    runs = {}
    port = {}
    chankeys = {}
    hsds = {}
    for r in runnums:
        dsdict.update({r:psana.DataSource(exp=expname,run=r)})
        runs.update({r:next(dsdict[r].runs())})
        port.update({r:{}})
        chankeys.update({r,[]})
        detslist.update({r:runs[r].detnames})
        if runhsd and 'hsd_mrco' in detlist[r]:
            hsds.update({r:runs[r].Detector('hsd_mrco')})
            for k in list(hsds[r].raw._seg_configs().keys()):
                chankeys[r] += [k]
                port[r].update({k:Port(k,chankeys[k],t0=t0s[k],logicthresh=logicthresh[k],inflate=inflate,expand=nr_expand)})
                port[r][k].setRollOn((3*int(hsds[r].raw._seg_configs()[k].config.user.fex.xpre))>>2)
                port[r][k].setRollOff((3*int(hsds[r].raw._seg_configs()[k].config.user.fex.xpost))>>2)
        else:
            runhsd = False

    _ = [print(d) for d in dslist]


    THIS IS GOING TO BE SUPER BROKEN... NOT LISTS ANYMORE!! DICTIONARIES




    for r in range(len(runnums)):
        chankeys.update({r,[]})
        eventnum:int = 0
        print('processing run %i'%runnums[r])
        if runfzp and 'tmo_fzppiranha' in runs[r].detnames:
            fzps += [runs[r].Detector('tmo_fzppiranha')]
        else:
            runfzp = False

        if runtiming and '' in runs[r].detnames:
            timings += [runs[r].Detector('timing')]
        else:
            runtiming = False

        if runvls and 'andor' in runs[r].detnames:
            vlss += [runs[r].Detector('andor')]
        else:
            runvls = False

        if runebeam and 'ebeam' in runs[r].detnames:
            ebeams += [runs[r].Detector('ebeam')]
        else:
            runebeam = False

        if runxtcav and 'xtcav' in runs[r].detnames:
            xtcavs += [runs[r].Detector('xtcav')]
        else:
            runxtcav = False

        if rungmd and 'xgmd' in runs[r].detnames:
            rungmd = True
            xgmds += [runs[r].Detector('xgmd')]
        else:
            rungmd = False

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

        completeEvent:List[bool] = [True]

        for evt in runs[r].events():
            if eventnum > nshots:
                break
            print(eventnum)
            ## test hsds
            if runhsd and bool(np.prod(completeEvent)):
                if type(hsds[r]) == None:
                    print(eventnum,'hsds is None')
                    completeEvent += [False]

                for key in chans.keys(): # here key means 'port number'
                    completeEvent += [port[r][key].test(hsds[r].raw.peaks(evt)[ chans[key] ][0])]
                    #completeEvent += [port[r][key].test(hsds[r].raw.waveforms(evt)[ chans[key] ][0])]


            ## finish testing all detectors to measure ##
            ## before processing ##

            ## process hsds
            if runhsd and bool(np.prod(completeEvent)):
                ''' HSD-Abaco section '''
                for key in chans.keys(): # here key means 'port number'
                    nwins:int = len(hsd.raw.peaks(evt)[ chans[key] ][0][1]
                    for i in range(nwins):
                        x = hsds[r].raw.peaks(evt)[ chans[key] ][0][0][i] 
                        s = np.array(hsds[r].raw.peaks(evt)[ chans[key] ][0][1][i],)
                        port[r][key].process(s,x)
                    #s = np.array(hsds[r].raw.waveforms(evt)[ chans[key] ][0] , dtype=np.int16) # presumably 12 bits unsigned input, cast as int16_t since will immediately in-place subtract baseline
                    #port[r][key].process(s,x=0).advance_event()
            eventnum += 1

        print('Finished with run %i'%runnums[r])
        print('returning')
        return

"""
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




## process VLS
            if runvls and bool(np.prod(completeEvent)):
                spect[r].process(vlswv)

## process ebeam
            if runebeam and bool(np.prod(completeEvent)):
                ebunch[r].process(thisl3)

## process gmd
            if rungmd and bool(np.prod(completeEvent)):
                gmd[r].process(thisgmde)





            eventnum += 1

    print("Hello, I'm done now.  Have a most excellent day!")
    return




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
