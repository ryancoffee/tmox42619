import numpy as np
from scipy.fftpack import dct,dst,rfft,irfft,fft,ifft
from utils import mypoly,tanhInt,tanhFloat,randomround
import h5py
import time
from typing import Type,List

    
#def dctLogic_windowed(s,inflate=1,nrolloff=0,winsz=256,stride=128):

def cfdLogic(s,invfrac=1<<2,offset=10):
    sz = s.shape[0]
    result = np.zeros(sz,dtype=np.int32)
    result[:-offset] = (1<<invfrac)*s[offset:]-s[:-offset]
    return result

def fftLogic_f16(s,inflate=1,nrolloff=128):
    sz = s.shape[0]
    result = np.zeros(sz*inflate,dtype=np.int16)
    rolloff_vec = (1<<3)*(1.+np.cos(np.arange(nrolloff<<1,dtype=float)*2*np.pi/float(nrolloff))) # careful, operating on left and right of middle indices in one go...2pi now not pi.
    smirror = np.append(s,np.flip(s,axis=0)).astype(np.int16)
    S = fft(smirror,axis=0)
    SR = np.copy(S.real).astype(np.int16)
    SI = np.copy(S.imag).astype(np.int16)
    SR[sz-nrolloff:sz+nrolloff] *= rolloff_vec.astype(np.int16)
    SI[sz-nrolloff:sz+nrolloff] *= rolloff_vec.astype(np.int16)
    SR[:sz-nrolloff] <<= 4
    SR[sz+nrolloff:] <<= 4
    SI[:sz-nrolloff] <<= 4
    SI[sz+nrolloff:] <<= 4
    S = SR + 1j*SI
    if inflate>1:
        S = np.concatenate((S[:sz],np.zeros(2*sz*(inflate-1),dtype=complex),S[sz:]))
    Sy = np.copy(S)*sz
    S[:sz] *= 1j*np.arange(sz,dtype=np.float16)/(sz)
    S[-sz:] *= np.flip(-1j*np.arange(sz,dtype=np.float16)/(sz),axis=0)
    y = ifft(Sy,axis=0).real[:(inflate*sz)].astype(np.int16)
    dy = ifft(S,axis=0).real[:(inflate*sz)].astype(np.int16)
    return -((y>>4)*(dy>>4)).astype(np.int16)


def fftLogic_fex(s,inflate=1,nrollon=8,nrolloff=32):
    sz = s.shape[0]
    result = np.zeros(sz*inflate,dtype=np.int32)
    rollon_vec = 0.5*(1.-np.cos(np.arange(nrollon,dtype=float)*np.pi/float(nrollon))) 
    rolloff_vec = 0.5*(1.+np.cos(np.arange(nrolloff<<1,dtype=float)*np.pi/float(nrolloff))) # careful, operating on left and right of middle indices in one go...2pi now not pi.
    smirror = np.append(s,np.flip(s,axis=0)).astype(float)
    smirror[:nrollon] *= rollon_vec
    smirror[-nrollon-1:] *= np.flip(rollon_vec,axis=0)
    smirror[(sz>>1)-nrolloff:(sz>>1)+nrolloff] *= rolloff_vec
    
    S = fft(smirror,axis=0)
    #S[sz-nrolloff:sz+nrolloff] *= rolloff_vec
    if inflate>1:
        S = np.concatenate((S[:sz],np.zeros(2*sz*(inflate-1),dtype=complex),S[sz:]))
    Sy = np.copy(S)*sz
    S[:sz] *= 1j*np.arange(sz,dtype=float)/(sz)
    S[-sz:] *= np.flip(-1j*np.arange(sz,dtype=float)/(sz),axis=0)
    y = ifft(Sy,axis=0).real[:(inflate*sz)]
    dy = ifft(S,axis=0).real[:(inflate*sz)]
    return -y*dy

def fftLogic(s,inflate=1,nrollon=64,nrolloff=128):
    sz = s.shape[0]
    result = np.zeros(sz*inflate,dtype=np.int32)
    rolloff_vec = 0.5*(1.+np.cos(np.arange(nrolloff<<1,dtype=float)*2*np.pi/float(nrolloff))) # careful, operating on left and right of middle indices in one go...2pi now not pi.
    smirror = np.append(s,np.flip(s,axis=0)).astype(float)
    S = fft(smirror,axis=0)
    S[sz-nrolloff:sz+nrolloff] *= rolloff_vec
    if inflate>1:
        S = np.concatenate((S[:sz],np.zeros(2*sz*(inflate-1),dtype=complex),S[sz:]))
    Sy = np.copy(S)*sz
    S[:sz] *= 1j*np.arange(sz,dtype=float)/(sz)
    S[-sz:] *= np.flip(-1j*np.arange(sz,dtype=float)/(sz),axis=0)
    y = ifft(Sy,axis=0).real[:(inflate*sz)]
    dy = ifft(S,axis=0).real[:(inflate*sz)]
    return -y*dy

"""
wv = hsd.raw.waveforms(evt)[1][0]
wvx = np.arange(wv.shape[0])
y = [hsd.raw.peaks(evt)[1][0][1][i] for i in range(len(hsd.raw.peaks(evt)[1][0][1]))]
x = [np.arange(hsd.raw.peaks(evt)[1][0][0][i],hsd.raw.peaks(evt)[1][0][0][i]+len(hsd.raw.peaks(evt)[1][0][1][i])) for i in range(len(y))]
plt.plot(wv)
_=[plt.plot(x[i],y[i]) for i in range(len(y))]
plt.show()
"""

class Port:
    # Note that t0s are aligned with 'prompt' in the digitizer logic signal
    # Don't forget to multiply by inflate, also, these look to jitter by up to 1 ns
    # hard coded the x4 scale-up for the sake of filling int16 dynamic range with the 12bit vls data and finer adjustment with adc offset correction

    def __init__(self,portnum,hsd,t0=0,nadcs=4,baselim=1000,logicthresh=-1*(1<<20),inflate=1,expand=1,nrollon=256,nrolloff=256): # exand is for sake of Newton-Raphson
        self.rng = np.random.default_rng( time.time_ns()%(1<<8) )
        self.portnum = portnum
        self.hsd = hsd
        self.t0 = t0
        self.nadcs = nadcs
        self.baselim = baselim
        self.logicthresh = logicthresh
        self.initState = True
        self.inflate = inflate
        self.expand = expand
        self.nrollon = nrollon
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
        self.processAlgos = 'fex2hits' # add as method to set the Algo for either of 'fex2hits', 'fex2coeffs', or just 'wave'

        self.e:List[np.uint32] = []
        self.de:List[np.int32] = []
        self.ne = 0
        self.r = []


    @classmethod
    def slim_update_h5(cls,f,port,hsdEvents,chans):
        for key in chans.keys(): # remember key == port number
            g = None
            if 'port_%i'%(key) in f.keys():
                g = f['port_%i'%(key)]
            else:
                g = f.create_group('port_%i'%(key))
            g.create_dataset('tofs',data=port[key].tofs,dtype=np.uint64) 
            g.create_dataset('slopes',data=port[key].slopes,dtype=np.int64) 
            g.create_dataset('addresses',data=port[key].addresses,dtype=np.uint64)
            g.create_dataset('nedges',data=port[key].nedges,dtype=np.uint64)
            g.attrs.create('inflate',data=port[key].inflate,dtype=np.uint8)
            g.attrs.create('expand',data=port[key].expand,dtype=np.uint8)
            g.attrs.create('t0',data=port[key].t0,dtype=float)
            g.attrs.create('logicthresh',data=port[key].logicthresh,dtype=np.int32)
            g.attrs.create('hsd',data=port[key].hsd,dtype=np.uint8)
            g.attrs.create('size',data=port[key].sz*port[key].inflate,dtype=np.uint64) ### need to also multiply by expand #### HERE HERE HERE HERE
            g.create_dataset('events',data=hsdEvents)
        return 

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
            g.create_dataset('tofs',data=port[key].tofs,dtype=np.uint64) 
            g.create_dataset('slopes',data=port[key].slopes,dtype=np.int64) 
            g.create_dataset('addresses',data=port[key].addresses,dtype=np.uint64)
            g.create_dataset('nedges',data=port[key].nedges,dtype=np.uint64)
            for k in port[key].waves.keys():
                rawgrp.create_dataset(k,data=port[key].raw[k].astype(np.uint16),dtype=np.uint16)
                wvgrp.create_dataset(k,data=port[key].waves[k].astype(np.int16),dtype=np.int16)
                lggrp.create_dataset(k,data=port[key].logics[k].astype(np.int32),dtype=np.int32)
            g.attrs.create('inflate',data=port[key].inflate,dtype=np.uint8)
            g.attrs.create('expand',data=port[key].expand,dtype=np.uint8)
            g.attrs.create('t0',data=port[key].t0,dtype=float)
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

    def addsample(self,o,w,l):
        eventnum = len(self.addresses)
        if eventnum<100:
            if eventnum%10<10: 
                self.raw.update( {'shot_%i'%eventnum:np.copy(o)} )
                self.waves.update( {'shot_%i'%eventnum:np.copy(w)} )
                self.logics.update( {'shot_%i'%eventnum:np.copy(l)} )
        elif eventnum<1000:
            if eventnum%100<10: 
                self.raw.update( {'shot_%i'%eventnum:np.copy(o)} )
                self.waves.update( {'shot_%i'%eventnum:np.copy(w)} )
                self.logics.update( {'shot_%i'%eventnum:np.copy(l)} )
        elif eventnum<10000:
            if eventnum%1000<10: 
                self.raw.update( {'shot_%i'%eventnum:np.copy(o)} )
                self.waves.update( {'shot_%i'%eventnum:np.copy(w)} )
                self.logics.update( {'shot_%i'%eventnum:np.copy(l)} )
        else:
            if eventnum%10000<10: 
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
            while d[i] > self.logicthresh:
                i += 1
                if i==sz-10: return tofs,slopes,len(tofs)
            while i<sz-10 and d[i]<0:
                i += 1
            stop = i
            ''' dx / (Dy) = dx2/dy2 ; dy2*dx/Dy - dx2 ; x2-dx2 = stop - dy2*1/Dy'''
            x0 = float(stop) - float(d[stop])/float(d[stop]-d[stop-1])
            i += 1
            v = float(self.expand)*float(x0)
            tofs += [np.uint32(randomround(v,self.rng))] 
            slopes += [d[stop]-d[stop-1]] 
        return tofs,slopes,np.uint64(len(tofs))

    def scanedges(self,d):
        tofs = []
        slopes = []
        sz = d.shape[0]
        newtloops = 6
        order = 3 # this should stay fixed, since the logic zeros crossings really are cubic polys
        i = 10
        while i < sz-10:
            while d[i] > self.logicthresh:
                i += 1
                if i==sz-10: return tofs,slopes,len(tofs)
            while i<sz-10 and d[i]<d[i-1]:
                i += 1
            start = i-1
            i += 1
            while i<sz-10 and d[i]>d[i-1]:
                i += 1
            stop = i
            i += 1
            if (stop-start)<(order+1):
                continue
            x = np.arange(stop-start,dtype=float) # set x to index values
            y = d[start:stop] # set y to vector values
            x0 = float(stop)/2. # set x0 to halfway point
            #y -= (y[0]+y[-1])/2. # subtract average (this gets rid of residual DC offsets)
    
            theta = np.linalg.pinv( mypoly(np.array(x).astype(float),order=order) ).dot(np.array(y).astype(float)) # fit a polynomial (order 3) to the points
            for j in range(newtloops): # 3 rounds of Newton-Raphson
                X0 = np.array([np.power(x0,int(k)) for k in range(order+1)])
                x0 -= theta.dot(X0)/theta.dot([i*X0[(k+1)%(order+1)] for k in range(order+1)]) # this seems like maybe it should be wrong
            tofs += [float(start + x0)] 
            #X0 = np.array([np.power(x0,int(i)) for k in range(order+1)])
            #slopes += [np.int64(theta.dot([i*X0[(i+1)%(order+1)] for i in range(order+1)]))]
            slopes += [float((theta[1]+x0*theta[2])/2**18)] ## scaling to reign in the obscene derivatives... probably shoul;d be scaling d here instead
        return tofs,slopes,np.uint32(len(tofs))

    def test(self,s):
        if type(s) == type(None):
            return False
        return True

    def process(self,s,x=0):
        if self.processAlgo =='fex2coeffs':
            return process_vfex2coeffs(s,x)
        elif self.processAlgo == 'fex2hits':
            return process_fex2hits(s,x)
        return process_wave(s,x=0)

    def process_fex2coeffs(self,s,x):
        print('HERE HERE HERE HERE')
        return True

    def advance_event(self):
        self.e = []
        self.de = []
        self.ne = 0
        self.r = []
        return self

    def process_fex2hits(self,slist,xlist):
        e = []
        de = []
        ne = 0
        r = []
        if type(s) == type(None):
            return False
        else:
            for i,s in enumerate(slist):
                if len(self.addresses)%100==0:
                    self.r = list(np.copy(s).astype(np.int16))
                ## no longer needing to correct for the adc offsets. ##
                logic = fftLogic_fex(s,inflate=self.inflate,nrollon=self.nrollon,nrolloff=self.nrolloff) #produce the "logic vector"
                e,de,ne = self.scanedges_simple(logic) # scan the logic vector for hits

            self.e += e
            self.de += de
            self.ne += ne

            print('NOT DONE HERE')
            if len(self.addresses)%100==0:
                self.addsample(r,s,logic)

        return True


    def process_wave(self,slist,xlist=[0]):
        e:List[np.int32] = []
        de = []
        ne = 0
        r = []
        s = slist[0]
        x = xlist[0]
        if type(s) == type(None):
            #self.addsample(np.zeros((2,),np.int16),np.zeros((2,),np.float16))
            e:List[np.int32] = []
            de = []
            ne = 0
            return False
        else:
            if len(self.addresses)%100==0:
                r = np.copy(s).astype(np.uint16)
            for adc in range(self.nadcs): # correcting systematic baseline differences for the four ADCs.
                b = np.mean(s[adc:self.baselim+adc:self.nadcs])
                s[adc::self.nadcs] = (s[adc::self.nadcs] ) - np.int32(b)
            #logic = fftLogic(s,inflate=self.inflate,nrolloff=self.nrolloff) #produce the "logic vector"
            logic = fftLogic_f16(s,inflate=self.inflate,nrolloff=self.nrolloff) #produce the "logic vector"
            e,de,ne = self.scanedges_simple(logic) # scan the logic vector for hits

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
        if len(self.addresses)%100==0:
            self.addsample(r,s,logic)
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
    def setRollOn(self,n):
        self.nrollon = n
        return self
    def setRollOff(self,n):
        self.nrolloff = n
        return self

