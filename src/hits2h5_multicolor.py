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


def main():
        ############################################
        ###### Change this to your output dir ######
        ############################################
    #scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output/h5files'
    #scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_2022/h5files'
    scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_multicolorhack/h5files'

    rng = np.random.default_rng(seed = int(time.time()%1*1e6))

    if len(sys.argv)<3:
        print('syntax: ./hits2h5_multicolor.py <nshots> <expname> <list of run numbers>')
    expname = sys.argv[1]
    nshots = int(sys.argv[2])
    runnums = [int(r) for r in sys.argv[3:]]

    print('starting analysis exp %s for runs '%(expname))
    print(' '.join([str(r) for r in runnums]) )
    cfgname = '%s/%s.hsdconfig.h5'%(scratchdir,expname)
    params = fillconfigs(cfgname)
    chans = params['chans']
    t0s = params['t0s']
    logicthresh = params['logicthresh']

    nr_expand = params['expand']

    spect = Vls()
    ebunch = Ebeam()
    port = {} 
    scale = int(1) # to better fill 16 bit int
    inflate = params['inflate']
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
    hsd = None
    vls = None 
    ebeam = None
    if runhsd and checkdet(runs,'hsd'):
        hsd = runs[0].Detector('hsd')
    if runvls and checkdet(runs,'andor'):
        vls = runs[0].Detector('andor')
    if runebeam and checkdet(runs,'ebeam'):
        ebeam = runs[0].Detector('ebeam')

    wv = {}
    wv_logic = {}
    v = [] # vls data matrix
    vc = [] # vls centroids vector
    vs = [] # vls sum is I think not used, maybe for normalization or used to be for integration and PDF sampling
    l3 = [] # e-beam l3 (linac 3) in GeV.

    init = True 
    vsize = 0

    vlsEvents = []
    hsdEvents = []


    
    print('chans: ',chans)
    for eventnum in range(nshots): # careful, going to pull shots as if from same event... so not processing evt by evt anymore
        evts = [next(r.events()) for r in runs]
        select = 1+int(rng.uniform()*len(runnums)-1)
        chooseevts = rng.choice(evts,select)

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
                spect.process_list([np.squeeze(vls.raw.value(evt)) for evt in chooseevts],max_len=len(runs))
                vlsEvents += [eventnum]
                #spect.print_v()
            except:
                print(eventnum,'skip per vls')
                continue

        if runebeam:
                ''' Ebeam specific section '''
                try:
                    ebunch.process_list([ebeam.raw.ebeamL3Energy(evt)+0.5 for evt in chooseevts],max_len=len(runs))
                except:
                    print(eventnum,'skipping ebeam, skip per l3')
                    continue

        if runhsd:
    
            ''' HSD-Abaco section '''
            for key in chans.keys(): # here key means 'port number'
                #try:
                ss = [np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=np.int16) for evt in chooseevts]
                port[key].process_list(ss,max_len=len(runs))


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

        if init:
            init = False
            if runebeam:
                ebunch.set_initState(False)
            if runvls:
                spect.set_initState(False)
            if runhsd:
                for key in chans.keys():
                    port[key].set_initState(False)
        eventnum += 1

    runstrings = ['%03i'%i for i in runnums]
    outname = '%s/hits.%s.runs_'%(scratchdir,expname) + '-'.join(runstrings) + '.h5'
    print('writing to %s'%outname)
        
    with h5py.File(outname,'w') as f:
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

    print("Hello, I'm done now!")
    return

if __name__ == '__main__':
    main()
