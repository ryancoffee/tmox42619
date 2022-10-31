import numpy as np
import typing
from typing import List

IntArray = List[int]


class Vls:
    def __init__(self) -> None:
        self.v = [[]]
        self.vsize = int(0)
        self.vc = [[]]
        self.vs = [[]]
        self.initState = True
        return

    def process_list(self, vlswvs: List[IntArray]):
        nums = [np.sum([i*vlswv[i] for i in range(len(vlswv))]) for vlswv in vlswvs ]
        dens = [np.sum(vlswv) for vlswv in vlswvs]
        if self.initState:
            self.v = [[np.sum(vlswvs,axis=1).astype(np.int16)]]
            self.vsize = len(self.v)
            self.vc = [[np.uint16(nums[i]/dens[i]) for i in range(len(nums))]]
            self.vs = [[np.uint64(d) for d in dens]]
            self.initState = False
        else:
            self.v += [[np.sum(vlswvs,axis=1).astype(np.int16)]]
            self.vc += [[np.uint16(nums[i]/dens[i]) for i in range(len(nums))]]
            self.vs += [[np.uint64(d) for d in dens]]
        return self

    def process(self, vlswv: IntArray):
        mean = int(np.mean(vlswv[1900:])) # this subtracts baseline
        vlswv -= mean #vlswv-int(np.mean(vlswv[1900:])) # this subtracts baseline
        #print("processing vls",vlswv.shape[0])
        num = np.sum(np.array([i*vlswv[i] for i in range(len(vlswv))]))
        den = np.sum(vlswv)
        if self.initState:
            self.v = [vlswv.astype(np.int16)]
            self.vsize = len(self.v)
            self.vc = [np.uint16(num/den)]
            self.vs = [np.uint64(den)]
            self.initState = False
        else:
            self.v += [vlswv.astype(np.int16)]
            self.vc += [np.uint16(num/den)]
            self.vs += [np.uint64(den)]
        return self

    def set_initState(self,state: bool):
        self.initState = state
        return self

    def print_v(self):
        print(self.v[:10])
        return self

