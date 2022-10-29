import numpy as np

class Ebeam:
    def __init__(self):
        self.l3 = [[]]
        self.initState = True
        return 
    def process_list(self,l3list):
        if self.initState:
            self.l3 = [[np.uint16(l3) for l3 in l3list]]
        else:
            self.l3 += [[np.uint16(l3) for l3 in l3list]]
        return self
    def process(self,l3in):
        if self.initState:
            self.l3 = [l3in]
        else:
            self.l3 += [np.uint16(l3in)]
        return self
    def set_initState(self,state):
        self.initState = state
        return self
