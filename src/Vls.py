import numpy as np
import typing
from typing import List
import h5py
import sys

def getcentroid(data,pct=.8):
    csum = np.cumsum(data.astype(float))
    s = float(csum[-1])*pct
    csum /= csum[-1]
    inds = np.where((csum>(.5-pct/2.))*(csum<(.5+pct/2.)))
    tmp = np.zeros(data.shape,dtype=float)
    tmp[inds] = data[inds].astype(float)
    num = np.sum(tmp*np.arange(data.shape[0],dtype=float))
    return (num/s,np.uint64(s))

class Vls:
    def __init__(self,thresh) -> None:
        self.v = []
        self.vsize = int(0)
        self.vc = [[]]
        self.vs = [[]]
        self.initState = True
        self.vlsthresh = thresh
        self.winstart = 0
        self.winstop = 1<<11
        return

    @classmethod
    def update_h5(cls,f,spect,vlsEvents):
        grpvls = None
        if 'vls' in f.keys():
            grpvls = f['vls']
        else:
            grpvls = f.create_group('vls')

        grpvls.create_dataset('data',data=spect.v,dtype=int)
        grpvls.create_dataset('centroids',data=spect.vc,dtype=np.float16)
        grpvls.create_dataset('sum',data=spect.vs,dtype=np.uint64)
        grpvls.attrs.create('size',data=spect.vsize,dtype=np.int32)
        grpvls.create_dataset('events',data=vlsEvents)
        return

    def setthresh(self,x):
        self.vlsthresh = x
        return self

    def process_list(self, vlswvs,max_len):
        nums = [np.sum([i*vlswv[i] for i in range(len(vlswv))]) for vlswv in vlswvs ]
        dens = [np.sum(vlswv) for vlswv in vlswvs]
        if self.initState:
            self.v = [np.sum(vlswvs,axis=0).astype(np.int16)]
            self.vsize = len(self.v)
            self.vc = [[np.uint16(nums[i]/dens[i]) for i in range(len(nums))] + [0 for i in range(max_len-len(nums))] ]
            self.vs = [[np.uint64(d) for d in dens] + [0 for i in range(max_len-len(nums))] ]
            self.initState = False
        else:
            self.v += [np.sum(vlswvs,axis=0).astype(np.int16)]
            self.vc += [[np.uint16(nums[i]/dens[i]) for i in range(len(nums))] + [0 for i in range(max_len-len(nums))] ]
            self.vs += [[np.uint64(d) for d in dens] + [0 for i in range(max_len-len(nums))] ]
        return self

    def setwin(self,low,high):
        self.winstart = int(low)
        self.winstop = int(high)
        return self

    def process(self, vlswv):
        mean = np.int16(0)
        try:
            mean = np.int16(np.mean(vlswv[1800:])) # this subtracts baseline
        except:
            print('Damnit, Vls!')
            return False
        else:
            if (np.max(vlswv)-mean)<self.vlsthresh:
                return False
            d = np.copy(vlswv-mean).astype(np.int16)
            c,s = getcentroid(d[self.winstart:self.winstop],pct=0.8)
            if self.initState:
                self.v = [d]
                self.vsize = len(self.v)
                self.vc = [np.float16(c)]
                self.vs = [np.uint64(s)]
                self.initState = False
            else:
                self.v += [d]
                self.vc += [np.float16(c)]
                self.vs += [np.uint64(s)]
        return True

    def set_initState(self,state: bool):
        self.initState = state
        return self

    def print_v(self):
        print(self.v[:10])
        return self

