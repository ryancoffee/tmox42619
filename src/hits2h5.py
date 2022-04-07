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


def PWRspectrum(wv):
    return np.power(abs(np.fft.fft(wv).real),int(2))

def rollon(vec,n):
    vec[:int(n)] = vec[:int(n)]*np.arange(int(n),dtype=float)/float(n)
    return vec

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
    scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_2022/h5files'
    expname = 'tmox42619'
    runnum = 62 
    nshots = 100
    if len(sys.argv)>2:
        expname = sys.argv[1]
        runnum = int(sys.argv[2])

    if len(sys.argv)>3:
        nshots = int(sys.argv[3])

    print('starting analysis exp %s for run %i'%(expname,int(runnum)))
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

    ds = psana.DataSource(exp=expname,run=runnum)

    #for run in ds.runs():
    run = next(ds.runs())
        #np.savetxt('%s/waveforms.%s.%i.%i.dat'%(scratchdir,expname,runnum,key),wv[key],fmt='%i',header=headstring)
    for i in range(1):
        print(run.detnames)
        eventnum = 0
        runhsd=True
        runvls=False
        runebeam=False
        hsd = run.Detector('hsd')
        vls = run.Detector('andor')
        ebeam = run.Detector('ebeam')
        wv = {}
        wv_logic = {}
        v = [] # vls data matrix
        vc = [] # vls centroids vector
        vs = [] # vls sum is I think not used, maybe for normalization or used to be for integration and PDF sampling
        l3 = [] # e-beam l3 (linac 3) in GeV.


        init = True 
        vsize = 0

        print('chans: ',chans)
        for evt in run.events():
            if eventnum > nshots:
                break

            if runvls:
                ''' VLS specific section, do this first to slice only good shots '''
                try:
                    vlswv = np.squeeze(vls.raw.value(evt))
                    vlswv = vlswv-int(np.mean(vlswv[1900:])) # this subtracts baseline
                    if np.max(vlswv)<300:  # too little amount of xrays
                        print(eventnum,'skip per weak vls')
                        #eventnum += 1
                        continue
                    spect.process(vlswv)
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

                if eventnum<100:
                    if eventnum%10<2: 
                        print('working event %i, nedges = %s'%(eventnum,[port[k].getnedges() for k in chans.keys()] ))
                elif eventnum<1000:
                    if eventnum%100<2: 
                        print('working event %i, nedges = %s'%(eventnum,[port[k].getnedges() for k in chans.keys()] ))
                else:
                    if eventnum%1000<2: 
                        print('working event %i, nedges = %s'%(eventnum,[port[k].getnedges() for k in chans.keys()] ))
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
        if runvls:
            grpvls = f.create_group('vls')
            grpvls.create_dataset('data',data=spect.v,dtype=np.int16)
            grpvls.create_dataset('centroids',data=spect.vc,dtype=np.int16)
            grpvls.create_dataset('sum',data=spect.vs,dtype=np.uint64)
            grpvls.attrs.create('size',data=spect.vsize,dtype=np.int32)
        if runebeam:
            grpebeam = f.create_group('ebeam')
            grpebeam.create_dataset('l3energy',data=ebunch.l3,dtype=np.uint16)
        f.close()

    print("Hello, I'm done now!")
    return

if __name__ == '__main__':
    main()
