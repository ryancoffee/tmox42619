#!/sdf/group/lcls/ds/ana/sw/conda2/inst/envs/ps-4.6.3/bin/python3
import h5py
import numpy as np
import sys

def main():
    chans = {0:3,1:9,2:11,4:10,5:12,12:5,13:6,14:8,15:2,16:13} # HSD to port number:hsd
    # second knee # logicthresh = {0:-2000000, 1:-1500000, 2:-800000, 4:-800000, 5:-2500000, 12:-3000000, 13:-2300000, 14:-2100000, 15:-2000000, 16:-3300000} # set by knee (log-log) in val histogram
    #logicthresh = {0:-(2**20), 1:-(2**20), 2:-(2**20), 4:-(2**20), 5:-(2**20), 12:-(2**20+2**18+2**17), 13:-(2**20+2**19), 14:-(2**20+2**18), 15:-(2**20), 16:-(2**21)} # set by 1st knee (log-log) in val histogram
    #logicthresh = {0:-1*((1<<18)), 1:-1*((1<<18)), 2:-1*((1<<18)+(1<<18)), 4:-1*((1<<18)), 5:-1*((1<<18)), 12:-1*((1<<18)), 13:-1*((1<<18)), 14:-1*((1<<18)), 15:-1*((1<<18)), 16:-1*((1<<18)+(1<<17))} # set by 1st knee (log-log) in val histogram
    #logicthresh = {0:-1*((1<<23)+(1<<17)), 1:-1*((1<<23)+(1<<17)), 2:-1*((1<<23)+(1<<17)), 4:-1*((1<<23)+(1<<17)), 5:-1*((1<<23)+(1<<17)), 12:-1*((1<<23)+(1<<18)), 13:-1*((1<<23)+(1<<18)), 14:-1*((1<<23)+(1<<18)), 15:-1*((1<<23)+(1<<17)), 16:-1*((1<<23))} # set by 1st knee (log-log) in val histogram


# for use with fftLogic_f16
    logicthresh = {0:-1*(1<<15), 1:-1*(1<<15), 2:-1*(1<<15), 4:-1*(1<<15), 5:-1*(1<<15), 12:-1*(1<<15), 13:-1*(1<<15), 14:-1*(1<<15), 15:-1*(1<<15), 16:-1*(1<<15)} # set by 1st knee (log-log) in val histogram
    _= [offsets.update({k:[0]*4}) for k in logicthresh.keys()]

    for k in logicthresh.keys():
        logicthresh[k] = logicthresh[k]>>2

    vlsthresh = 1000
    vlswin = (1024,2048)
    l3offset = 5100
    t0s = {0:4577,1:4186,2:4323,4:4050,5:4128,12:4107,13:4111,14:4180,15:4457,16:4085} # these are not accounting for the expand nor inflate, digitizer units, 6GSps, so 6k = 1usec

    if len(sys.argv)<2:
        print('I need an output filename to write configs to')
        return

    cfgfile = sys.argv[1]
    with h5py.File(cfgfile,'w') as f:
        f.attrs.create('expand',4) # expand controls the fractional resolution for scanedges by scaling index values and then zero crossing round to intermediate integers.
        f.attrs.create('inflate',2) # inflate pads the DCT(FFT) with zeros, artificially over sampling the waveform
        f.attrs.create('vlsthresh',data=vlsthresh)
        f.attrs.create('vlswin',data=vlswin)
        f.attrs.create('l3offset',data=l3offset)
        for k in chans.keys():
            key = 'port_%i'%int(k)
            c = f.create_group(key)
            c.attrs.create('hsd',data=chans[k])
            c.attrs.create('t0',data=t0s[k])
            c.attrs.create('logicthresh',data=logicthresh[k])
            c.attrs.create('offsets',data=offsets[k])
    return

if __name__ == '__main__':
    main()

