#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.2.5/bin/python3


import psana
import numpy as np
import sys
import re
import h5py
from scipy.fftpack import dct,dst
import time

from Ports import *
from Ebeam import *
from Vls import *
from Gmd import *
from utils import *
from config import read_config


def main():
        ############################################
        ###### Change this to your output dir ######
        ############################################
    #scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output/h5files'
    #scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_2022/h5files'
    scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_multicolorhack/h5files'

    rng = np.random.default_rng(seed = int(time.time()%1*1e6))

    if len(sys.argv)<3:
        print('syntax: ./hits2h5_twocolor.py <nshots> <expname> <list of run numbers>')
    expname = sys.argv[1]
    nshots = int(sys.argv[2])
    runnums = [int(r) for r in sys.argv[3:]]

    print('starting analysis exp %s for runs '%(expname))
    print(' '.join([str(r) for r in runnums]) )
    cfgname = '%s/%s.hsdconfig.h5'%(scratchdir,expname)
    params = read_config(cfgname)
    chans = params.chans
    t0s = params.t0s
    logicthresh = params.logicthresh

    nr_expand = params.expand

    spect = Vls()
    ebunch = Ebeam()
    port = {} 
    scale = int(1) # to better fill 16 bit int
    inflate = params.inflate
    for key in logicthresh.keys():
        logicthresh[key] *= scale # inflating by factor of 4 since we are also scaling the waveforms by 4 in vertical to fill bit depth.

    for key in chans.keys():
        port[key] = Port(key,chans[key],t0=t0s[key],logicthresh=logicthresh[key],inflate=inflate,expand=nr_expand,scale=scale,nrolloff=2**6)

    dss = [psana.DataSource(exp=expname,run=r) for r in runnums]

    #for run in ds.runs():
    runs = [next(ds.runs()) for ds in dss]
        #np.savetxt('%s/waveforms.%s.%i.%i.dat'%(scratchdir,expname,runnum,key),wv[key],fmt='%i',header=headstring)
    for r in runs:
        print(r.detnames)
    eventnum = 0
    runhsd=True
    runvls=True
    runebeam=True
    rungmd=False
    runxtcav=False
    hsd = None
    vls = None 
    ebeam = None
    gmd = None
    if runhsd and checkdet(runs,'hsd'):
        hsd = runs[0].Detector('hsd')
    if runvls and checkdet(runs,'andor'):
        vls = runs[0].Detector('andor')
    if runebeam and checkdet(runs,'ebeam'):
        ebeam = runs[0].Detector('ebeam')
    if rungmd and checkdet(runs,'gmd'):
        gmd = runs[0].Detector('gmd')

    wv = {}
    wv_logic = {}
    v = [] # vls data matrix
    vc = [] # vls centroids vector
    vs = [] # vls sum is I think not used, maybe for normalization or used to be for integration and PDF sampling
    l3 = [] # e-beam l3 (linac 3) in GeV.
    gm = [] # gmd energy in unknown

    init = True 
    vsize = 0

    vlsEvents = []
    hsdEvents = []
    ebeamEvents = []
    gmdEvents = []
    runstrings = ['%03i'%i for i in runnums]
    outname = '%s/hits.%s.runs_'%(scratchdir,expname) + '-'.join(runstrings) + '.h5'
    print('chans: ',chans)
    for eventnum in range(nshots): # careful, going to pull shots as if from same event... so not processing evt by evt anymore
        evts = [next(r.events()) for r in runs]

        if rungmd:
            ''' GMD specific section '''
            try:
                goodchoice = True
                ens = [gmd.raw.energy(evt) for evt in evts]
                for e in ens:
                    if e==None or str(e)=='nan':
                        goodchoice = False

                if goodchoice:
                    gmd.process_list(ens,max_len=len(runs))
                    gmdEvents += [eventnum]
                else:
                    print(eventnum,"skipping since badchoice on evts")
                    continue
            except:
                print(eventnum,"skipping, GMD failed as None or 'nan' for unkiwn opaque reason")
                continue

        if runvls:
            ''' VLS specific section, do this first to slice only good shots '''
            try:
                if type(vls) == None:
                    print(eventnum,'skip per problem with VLS')
                    continue
                #if np.max(vlswv)<1:  # too little amount of xrays
                #    print(eventnum,'skip per negative vls')
                #    #eventnum += 1
                #    continue
                spect.process_list([np.squeeze(vls.raw.value(evt)) for evt in evts],max_len=len(runs))
                vlsEvents += [eventnum]
                #spect.print_v()
            except:
                print(eventnum,'skip per vls')
                continue

        if runebeam:
                ''' Ebeam specific section '''
                try:
                    ebunch.process_list([ebeam.raw.ebeamL3Energy(evt) for evt in evts],max_len=len(runs))
                    ebeamEvents += [eventnum]
                except:
                    print(eventnum,'skipping ebeam, skip per l3')
                    continue



        if runhsd:
            ''' HSD-Abaco section '''
            for key in chans.keys(): # here key means 'port number'
                #try:
                ss = [np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=np.int16) for evt in evts]
                port[key].process_list(ss,max_len=len(runs))

            hsdEvents += [eventnum]

            if eventnum<10:
                print('ports = %s'%([k for k in chans.keys()]))
            if eventnum<100:
                if eventnum%10<1: 
                    print('working event %i,\tnedges = %s'%(eventnum,[port[k].getnedges() for k in chans.keys()] ))
            elif eventnum<1000:
                if eventnum%100<1: 
                    print('working event %i,\tnedges = %s'%(eventnum,[port[k].getnedges() for k in chans.keys()] ))
            else:
                if eventnum%1000<1: 
                    print('working event %i,\tnedges = %s'%(eventnum,[port[k].getnedges() for k in chans.keys()] ))



        if eventnum > 1 and eventnum <1000 and eventnum%100==0:
            with h5py.File(outname,'w') as f:
                print('writing to %s'%outname)
                if runhsd:
                    Port.update_h5(f,port,hsdEvents,chans)
                if runvls:
                    Vls.update_h5(f,spect,vlsEvents)
                if runebeam:
                    Ebeam.update_h5(f,ebunch,ebeamEvents)
                if rungmd:
                    Gmd.update_h5(f,gmd,gmdEvents)

        elif eventnum>900 and eventnum%1000==0:
            with h5py.File(outname,'w') as f:
                print('writing to %s'%outname)
                if runhsd:
                    Port.update_h5(f,port,hsdEvents,chans)
                if runvls:
                    Vls.update_h5(f,spect,vlsEvents)
                if runebeam:
                    Ebeam.update_h5(f,ebunch,ebeamEvents)
                if rungmd:
                    Gmd.update_h5(f,gmd,gmdEvents)

        if init:
            init = False
            if rungmd:
                gmd.set_initState(False)
            if runebeam:
                ebunch.set_initState(False)
            if runvls:
                spect.set_initState(False)
            if runhsd:
                for key in chans.keys():
                    port[key].set_initState(False)
        eventnum += 1

        
    with h5py.File(outname,'w') as f:
        print('writing to %s'%outname)
        if runhsd:
            Port.update_h5(f,port,hsdEvents,chans)
        if runvls:
            Vls.update_h5(f,spect,vlsEvents)
        if runebeam:
            Ebeam.update_h5(f,ebunch,ebeamEvents)
        if rungmd:
            Gmd.update_h5(f,gmd,gmdEvents)

    print("Hello, I'm done now!")
    return

if __name__ == '__main__':
    main()
