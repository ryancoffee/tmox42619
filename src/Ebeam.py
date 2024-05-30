import numpy as np

class Ebeam:
    def __init__(self):
        self.l3 = []
        self.l3offset=5100
        self.initState = True
        return 

    @classmethod
    def slim_update_h5(cls,f,ebunch,ebeamEvents):
        grpebeam = None
        if 'ebeam' in f.keys():
            grpebeam = f['ebeam']
        else:
            grpebeam = f.create_group('ebeam')
        d = grpebeam.create_dataset('l3energy',data=ebunch.l3,dtype=np.float16)
        d.attrs.create('l3offset',ebunch.l3offset,dtype=np.uint16)
        return

    @classmethod
    def update_h5(cls,f,ebunch,ebeamEvents):
        grpebeam = None
        if 'ebeam' in f.keys():
            grpebeam = f['ebeam']
        else:
            grpebeam = f.create_group('ebeam')
        d = grpebeam.create_dataset('l3energy',data=ebunch.l3,dtype=np.float16)
        d.attrs.create('l3offset',ebunch.l3offset,dtype=np.uint16)
        return

    def setoffset(self,x):
        self.l3offset = int(x)
        return self

    def test(self,l3in):
        if type(l3in)==type(None):
            return False
        try:
            d = np.float16(float(l3in)-float(self.l3offset))
        except:
            print('Damnit, Ebeam!')
            return False
        return True

    def process(self,l3in):
        d = np.float16(float(l3in)-float(self.l3offset))
        if self.initState:
            self.l3 = [d]
        else:
            self.l3 += [d]
        return True

    def set_initState(self,state):
        self.initState = state
        return self
