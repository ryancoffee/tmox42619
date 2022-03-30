#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.2.5/bin/python3

            #&0,3,TMO:MPOD:01:M2:C0,TMO:MPOD:01:M9:C0,TMO:MPOD:01:M5:C0,TMO:MPOD:01:M0:C0,TMO:MPOD:01:M1:C0,TMO:MPOD:01:M6:C0,TMO:MPOD:01:M7:C0,TMO:MPOD:01:M3:C0,
            #&1,9,TMO:MPOD:01:M2:C1,TMO:MPOD:01:M9:C1,TMO:MPOD:01:M5:C1,TMO:MPOD:01:M0:C1,TMO:MPOD:01:M1:C1,TMO:MPOD:01:M6:C1,TMO:MPOD:01:M7:C1,TMO:MPOD:01:M3:C1,
            #&2,11,TMO:MPOD:01:M2:C2,TMO:MPOD:01:M9:C2,TMO:MPOD:01:M5:C2,NA,NA,NA,NA,NA,
            #&4,10,TMO:MPOD:01:M2:C4,TMO:MPOD:01:M9:C4,TMO:MPOD:01:M5:C4,TMO:MPOD:01:M0:C4,TMO:MPOD:01:M1:C4,TMO:MPOD:01:M6:C4,TMO:MPOD:01:M7:C4,TMO:MPOD:01:M3:C4,
            #&5,12,TMO:MPOD:01:M2:C5,TMO:MPOD:01:M9:C5,TMO:MPOD:01:M5:C5,TMO:MPOD:01:M0:C5,TMO:MPOD:01:M1:C5,TMO:MPOD:01:M6:C5,TMO:MPOD:01:M7:C5,TMO:MPOD:01:M3:C5,
            #&12,5,TMO:MPOD:01:M2:C12,TMO:MPOD:01:M9:C12,TMO:MPOD:01:M5:C12,TMO:MPOD:01:M0:C12,TMO:MPOD:01:M1:C12,TMO:MPOD:01:M6:C12,TMO:MPOD:01:M7:C12,TMO:MPOD:01:M3:C12,
            #&13,6,TMO:MPOD:01:M2:C13,TMO:MPOD:01:M9:C13,TMO:MPOD:01:M5:C13,TMO:MPOD:01:M0:C13,TMO:MPOD:01:M1:C13,TMO:MPOD:01:M6:C13,TMO:MPOD:01:M7:C13,TMO:MPOD:01:M3:C13,
            #&14,8,TMO:MPOD:01:M2:C14,TMO:MPOD:01:M9:C14,TMO:MPOD:01:M5:C14,TMO:MPOD:01:M0:C14,TMO:MPOD:01:M1:C14,TMO:MPOD:01:M6:C14,TMO:MPOD:01:M7:C14,TMO:MPOD:01:M3:C14,
            #&15,2,TMO:MPOD:01:M2:C15,TMO:MPOD:01:M9:C15,TMO:MPOD:01:M5:C15,TMO:MPOD:01:M0:C15,TMO:MPOD:01:M1:C15,TMO:MPOD:01:M6:C15,TMO:MPOD:01:M7:C15,TMO:MPOD:01:M3:C15,
            #&16,13,TMO:MPOD:01:M2:C16,TMO:MPOD:01:M9:C6,TMO:MPOD:01:M5:C6,NA,NA,NA,NA,NA,

import psana
import numpy as np
import sys
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
    nr_expand = 4 
    chans = {0:3,1:9,2:11,4:10,5:12,12:5,13:6,14:8,15:2,16:13} # HSD to port number:hsd
    logicthresh = {0:-800000, 1:-800000, 2:-400000, 4:-800000, 5:-800000, 12:-800000, 13:-800000, 14:-800000, 15:-800000, 16:-400000}
    #slopethresh = {0:500,1:500,2:300,4:150,5:500,12:500,13:500,14:500,15:500,16:300}
    slopethresh = {0:100,1:100,2:60,4:100,5:100,12:100,13:100,14:100,15:100,16:60}
    #t0s = {0:109840,1:100456,2:99924,4:97180,5:99072,12:98580,13:98676,14:100348,15:106968,16:98028}
    #t0s = {0:109830,1:100451,2:99810,4:97180,5:99071,12:98561,13:98657,14:100331,15:106956,16:97330}
    t0s = {0:73227,1:66793,2:60000,4:64796,5:66054,12:65712,13:65777,14:66891,15:71312,16:60000} # final, based on inflate=4 nr_expand=4
    print('need to redo t0s and logicthresh\n%s\n%s'%(t0s,logicthresh))

    '''
    argon   prompt>300      proposed
    0       109500  109830
    1       100121  100451
    2       99480   99810
    4       96850   97180
    5       98741   99071
    12      98231   98561
    13      98327   98657
    14      100001  100331
    15      106626  106956
    16      97000   97330
    '''


    spect = Vls()
    ebunch = Ebeam()
    port = {} 
    scale = int(1) # to better fill 16 bit int
    inflate = int(4) 
    for key in logicthresh.keys():
        logicthresh[key] *= scale # inflating by factor of 4 since we are also scaling the waveforms by 4 in vertical to fill bit depth.

    for key in chans.keys():
        port[key] = Port(key,chans[key],t0=t0s[key],logicthresh=logicthresh[key],slopethresh=slopethresh[key],inflate=inflate,expand=nr_expand,scale=scale,nrolloff=10000)

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
                    try:
                        s = np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=np.int16) 
                        port[key].process(s)

                        if init:
                            init = False
                            ebunch.set_initState(False)
                            spect.set_initState(False)
                            for key in chans.keys():
                                port[key].set_initState(False)
                    except:
                        print(eventnum, 'failed hsd for some reason')
                        continue

                if eventnum<100:
                    if eventnum%10<2: 
                        print('working event %i'%eventnum)
                elif eventnum<1000:
                    if eventnum%100<2: 
                        print('working event %i'%eventnum)
                else:
                    if eventnum%1000<2: 
                        print('working event %i'%eventnum)
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
                for k in port[key].waves.keys():
                    wvgrp.create_dataset(k,data=port[key].waves[k],dtype=np.int32)
                g.attrs.create('inflate',data=port[key].inflate,dtype=np.uint8)
                g.attrs.create('expand',data=port[key].expand,dtype=np.uint8)
                g.attrs.create('t0',data=port[key].t0,dtype=float)
                g.attrs.create('slopethresh',data=port[key].slopethresh,dtype=np.uint64)
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
