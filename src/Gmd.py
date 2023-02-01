import numpy as np

class Gmd:
    def __init__(self):
        self.en = [[]]
        self.initState = True
        return 

    @classmethod
    def update_h5(cls,f,gmd,gmdEvents):
        grpgmd = None
        if 'gmd' in f.keys():
            grpgmd = f['gmd']
        else:
            grpgmd = f.create_group('gmd')
        grpgmd.create_dataset('gmdenergy',data=gmd.en,dtype=np.uint16)
        grpgmd.create_dataset('events',data=gmdEvents)
        return

    def process_list(self,enlist,max_len):
        if self.initState:
            self.en = [[np.uint16(e*1000) for e in enlist] + [np.uint16(0) for i in range(max_len-len(enlist))]]
            self.initState = False
        else:
            self.en += [[np.uint16(e*1000) for e in enlist] + [np.uint16(0) for i in range(max_len-len(enlist))]]
        return self

    def process(self,e):
        if e<0:
            return False
        if self.initState:
            self.en = [np.uint16(e*1000)]
        else:
            self.en += [np.uint16(e*1000)]
        return True

    def set_initState(self,state):
        self.initState = state
        return self
