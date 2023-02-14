#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.2.5/bin/python3
import h5py
import numpy as np
import sys

def main(path,expname):
    chans = {0:3,1:9,2:11,4:10,5:12,12:5,13:6,14:8,15:2,16:13} # HSD to port number:hsd
    logicthresh = {0:((1<<23)+(1<<17)), 1:((1<<23)+(1<<17)), 2:((1<<23)+(1<<17)), 4:((1<<23)+(1<<17)), 5:((1<<23)+(1<<17)), 12:((1<<23)+(1<<18)), 13:((1<<23)+(1<<18)), 14:((1<<23)+(1<<18)), 15:((1<<23)+(1<<17)), 16:((1<<23))} # set by 1st knee (log-log) in val histogram
    for k in logicthresh.keys():
        logicthresh[k] = logicthresh[k]>>2


    cfgfile = '%s/%s.configs.h5'%(path,expname)
    with h5py.File(cfgfile,'w') as f:
        f.attrs.create('expand',4) # expand controls the fractional resolution for scanedges by scaling index values and then zero crossing and random round to intermediate integers.
        f.attrs.create('inflate',2) # inflate pads the DCT(FFT) with zeros, artificially over sampling the waveform
        for k in chans.keys():
            key = 'port_%i'%int(k)
            c = f.create_group(key)
            c.attrs.create('hsd',data=chans[k])
            c.attrs.create('logicthresh',data=logicthresh[k])
    return

if __name__ == '__main__':
    if len(sys.argv)<3:
        print('syntax:set_configs.py <path> <expname=tmox42619>')
    else:
        main(sys.argv[1],sys.argv[2])

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
