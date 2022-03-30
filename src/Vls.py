import numpy as np


class Vls:
    def __init__(self):
        self.v = []
        self.vsize = int(0)
        self.vc = []
        self.vs = []
        self.initState = True
        return

    def process(self, vlswv):
        #print("processing vls",vlswv.shape[0])
        num = np.sum(np.array([i*vlswv[i] for i in range(len(vlswv))]))
        den = np.sum(vlswv)
        if self.initState:
            self.v = [vlswv.astype(np.int16)]
            self.vsize = len(self.v)
            self.vc = [np.uint16(num/den)]
            self.vs = [np.uint64(den)]
        else:
            self.v += [vlswv.astype(np.int16)]
            self.vc += [np.uint16(num/den)]
            self.vs += [np.uint64(den)]
        return self

    def set_initState(self,state):
        self.initState = state
        return self

    def print_v(self):
        print(self.v[:10])
        return self

