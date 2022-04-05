#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.2.5/bin/python3
import h5py
import numpy as np
import sys

def main():
    chans = {0:3,1:9,2:11,4:10,5:12,12:5,13:6,14:8,15:2,16:13} # HSD to port number:hsd
    # second knee # logicthresh = {0:-2000000, 1:-1500000, 2:-800000, 4:-800000, 5:-2500000, 12:-3000000, 13:-2300000, 14:-2100000, 15:-2000000, 16:-3300000} # set by knee (log-log) in val histogram
    logicthresh = {0:-(2**20), 1:-(2**20), 2:-(2**20), 4:-(2**20), 5:-(2**20), 12:-(2**20+2**18), 13:-(2**20+2**19), 14:-(2**20+2**18), 15:-(2**20), 16:-(2**21)} # set by 1st knee (log-log) in val histogram
    t0s = {0:73246,1:67008,2:68300,4:64922,5:66075,12:65762,13:65827,14:66906,15:71400,16:65391} # based on latest in repo ryan-dev (inflate=4 nr_expand=4)
    # Ideally we would measure the logicthresh knee for different delay windows, as the low energy hits might have lower carrier cascade in MCPs

    if len(sys.argv)>1:
        cfgfile = sys.argv[1]
        with h5py.File(cfgfile,'w') as f:
            for k in chans.keys():
                key = 'port_%i'%int(k)
                c = f.create_group(key)
                c.attrs.create('hsd',data=chans[k])
                c.attrs.create('t0',data=t0s[k])
                c.attrs.create('logicthresh',data=logicthresh[k])
    return

if __name__ == '__main__':
    main()


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


    # old values logicthresh = {0:-800000, 1:-800000, 2:-400000, 4:-800000, 5:-800000, 12:-800000, 13:-800000, 14:-800000, 15:-800000, 16:-400000}
    # old values t0s = {0:73227,1:66793,2:60000,4:64796,5:66054,12:65712,13:65777,14:66891,15:71312,16:60000} # final, based on inflate=4 nr_expand=4

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
