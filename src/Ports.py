import numpy as np
from scipy.fftpack import fft,ifft
from utils import mypoly,tanhInt,tanhFloat,randomround
import h5py
import time
from typing import Type,List

    
def fftLogic(s,inflate=1,nrolloff=128):
    sz = s.shape[0]
    result = np.zeros(sz*inflate,dtype=np.int32)
    rolloff_vec = 0.5*(1.+np.cos(np.arange(nrolloff<<1,dtype=float)*2*np.pi/float(nrolloff))) # careful, operating on left and right of middle indices in one go...2pi now not pi.
    rolloff_vec *= 8.0
    rolloff_int_vec = np.array([np.int32(v) for v in rolloff_vec])
    S = ((1<<3)*fft(s,axis=0))
    Sreal = np.copy(S.real).astype(np.int32)
    Simag = np.copy(S.imag).astype(np.int32)
    Sreal[sz//2-nrolloff:sz//2+nrolloff] *= rolloff_int_vec
    Simag[sz//2-nrolloff:sz//2+nrolloff] *= rolloff_int_vec
    Sreal[:sz//2-nrolloff] <<=3
    Simag[:sz//2-nrolloff] <<=3
    Sreal[-sz//2+nrolloff:] <<=3
    Simag[-sz//2+nrolloff:] <<=3
    if inflate>1:
        Sreal = np.concatenate((Sreal[:sz//2],np.zeros(2*sz*(inflate-1),dtype=np.int32),Sreal[sz//2:]))
        Simag = np.concatenate((Simag[:sz//2],np.zeros(2*sz*(inflate-1),dtype=np.int32),Simag[sz//2:]))
    Sdyimag = np.zeros(Sreal.shape,dtype=np.int32)
    Sdyreal = np.zeros(Sreal.shape,dtype=np.int32)
    Sdyimag[:sz//2] = Sreal[:sz//2]*np.arange(sz//2,dtype=np.int32)
    Sdyimag[:sz//2] = Sreal[:sz//2]*np.arange(sz//2,dtype=np.int32)
    Sdyreal[-sz//2:] = -1*Simag[-sz//2:]*np.arange(-sz//2,0,1,dtype=np.int32)
    Sdyreal[-sz//2:] = -1*Simag[-sz//2:]*np.arange(-sz//2,0,1,dtype=np.int32)

    y = (ifft( Sreal.astype(complex) + 1j*Simag.astype(complex) ,axis=0).real[:(inflate*sz)]).astype(np.int32)
    dy = (ifft((Sdyreal.astype(complex) + 1j*Sdyimag.astype(complex)).astype(complex),axis=0).real[:(inflate*sz)]).astype(np.int32)
    retvec = (y*dy)>>6
    return retvec


class Port:
    # Don't forget to multiply by inflate, also, these look to jitter by up to 1 ns
    # hard coded the x4 scale-up for the sake of filling int16 dynamic range with the 12bit vls data and finer adjustment with adc offset correction

    def __init__(self,portnum,hsd,nadcs=4,baselim=(1<<8),logicthresh=-1*(1<<20),inflate=(1<<1),expand=(1<<2),nrolloff=(1<<7)): # exand is for sake of Newton-Raphson
        self.rng = np.random.default_rng( time.time_ns()%(1<<8) )
        self.portnum = portnum
        self.hsd = hsd
        self.nadcs = nadcs
        self.baselim = baselim
        self.logicthresh = logicthresh
        self.initState = True
        self.inflate = inflate
        self.expand = expand
        self.nrolloff = nrolloff
        self.sz = 0
        self.tofs = []
        self.slopes = []
        self.addresses = []
        self.nedges = []
        self.raw = {}
        self.waves = {}
        self.logics = {}
        self.shot = int(0)

    @classmethod
    def update_h5(cls,f,port,hsdEvents,chans):
        for key in chans.keys(): # remember key == port number
            g = None
            if 'port_%i'%(key) in f.keys():
                g = f['port_%i'%(key)]
                rawgrp = g['raw']
                wvgrp = g['waves']
                lggrp = g['logics']
            else:
                g = f.create_group('port_%i'%(key))
                rawgrp = g.create_group('raw')
                wvgrp = g.create_group('waves')
                lggrp = g.create_group('logics')
            g.create_dataset('tofs',data=port[key].tofs,dtype=np.uint32) 
            g.create_dataset('slopes',data=port[key].slopes,dtype=np.int64) 
            g.create_dataset('addresses',data=port[key].addresses,dtype=np.uint64)
            g.create_dataset('nedges',data=port[key].nedges,dtype=np.uint64)
            for k in port[key].waves.keys():
                rawgrp.create_dataset(k,data=port[key].raw[k].astype(np.uint16),dtype=np.uint16)
                wvgrp.create_dataset(k,data=port[key].waves[k].astype(np.int16),dtype=np.int16)
                lggrp.create_dataset(k,data=port[key].logics[k].astype(np.int16),dtype=np.int16)
            g.attrs.create('inflate',data=port[key].inflate,dtype=np.uint8)
            g.attrs.create('expand',data=port[key].expand,dtype=np.uint8)
            g.attrs.create('logicthresh',data=port[key].logicthresh,dtype=np.int32)
            g.attrs.create('hsd',data=port[key].hsd,dtype=np.uint8)
            g.attrs.create('size',data=port[key].sz*port[key].inflate,dtype=np.uint64) ### need to also multiply by expand #### HERE HERE HERE HERE
            g.create_dataset('events',data=hsdEvents)
        return 

    def addeverysample(self,o,w,l):
        eventnum = len(self.addresses)
        self.raw.update( {'shot_%i'%eventnum:np.copy(o)} )
        self.waves.update( {'shot_%i'%eventnum:np.copy(w)} )
        self.logics.update( {'shot_%i'%eventnum:np.copy(l)} )
        return self


    def scanedges_simple(self,d):
        tofs = []
        slopes = []
        sz = d.shape[0]
        i:int = int(10)
        while i < sz-10:
            while d[i] < self.logicthresh:
                i += 1
                if i==sz-10: return tofs,slopes,len(tofs)
            while i<sz-10 and d[i]>0:
                i += 1
            stop = i
            ''' dx / (Dy) = dx2/dy2 ; dy2*dx/Dy - dx2 ; x2-dx2 = stop - dy2*1/Dy'''
            x0 = float(stop) - float(d[stop])/float(d[stop]-d[stop-1])
            i += 1
            v = float(self.expand)*float(x0)
            tofs += [np.uint32(randomround(v,self.rng))] 
            slopes += [d[stop]-d[stop-1]] 
        return tofs,slopes,np.uint32(len(tofs))


    def process(self,s):
        e:List[np.int32] = []
        de = []
        ne = 0
        if type(s) == type(None):
            self.addeverysample(np.zeros((2,),np.uint16),np.zeros((2,),np.int16),np.zeros((2,),np.float16))
            e:List[np.int32] = []
            de = []
            ne = 0
        else:
            s_orig = np.copy(s[:((1<<13)+(1<<8))])
            for adc in range(self.nadcs):
                b = np.mean(s[adc:self.baselim+adc:self.nadcs])
                s[adc::self.nadcs] = (s[adc::self.nadcs] ) - np.int16(b)
            logic = fftLogic(s,inflate=self.inflate,nrolloff=self.nrolloff) #produce the "logic vector"
            e,de,ne = self.scanedges_simple(logic) # scan the logic vector for hits
            self.addeverysample(s_orig,s,logic)

        if self.initState:
            self.sz = s.shape[0]*self.inflate*self.expand
            self.addresses = [np.uint64(0)]
            self.nedges = [np.uint64(ne)]
            if ne>0:
                self.tofs += e
                self.slopes += de
        else:
            self.addresses += [np.uint64(len(self.tofs))]
            self.nedges += [np.uint64(ne)]
            if ne>0:
                self.tofs += e
                self.slopes += de
        return True

    def set_initState(self,state=True):
        self.initState = state
        return self

    def print_tofs(self):
        print(self.tofs)
        print(self.slopes)
        return self

    def getnedges(self):
        if len(self.nedges)==0:
            return 0
        return self.nedges[-1]

