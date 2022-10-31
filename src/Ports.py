import numpy as np
from scipy.fftpack import dct,dst
from utils import mypoly,tanhInt
    
#def dctLogic_windowed(s,inflate=1,nrolloff=0,winsz=256,stride=128):
rng = np.random.default_rng()

def dctLogicInt(s,inflate=1,nrolloff=128):
    '''
    Cutting off the tail of the y before back transforming and multiplying
    This is to reduce distortion of the dy signal when suplressing noise by multiplying by signal.
    '''
    sz = s.shape[0]
    result = np.zeros(sz*inflate,dtype=np.int32)
    ampscale = 2**8
    rolloff_vec = (ampscale*(1.+np.cos(np.arange(nrolloff,dtype=float)*np.pi/float(nrolloff)))).astype(np.int32)
    sc = np.append(s,np.flip(s,axis=0)).astype(np.int32)
    ss = np.append(s,np.flip(-1*s,axis=0)).astype(np.int32)
    wc = dct(sc,type=2,axis=0).astype(np.int64)
    ws = dst(ss,type=2,axis=0).astype(np.int64)
    wc[-nrolloff:] *= rolloff_vec
    ws[-nrolloff:] *= rolloff_vec
    wc[:-nrolloff] *= ampscale # scaling since we are keeping to int32
    ws[:-nrolloff] *= ampscale
    wy = np.copy(wc)
    if inflate>1: # inflating seems to increase the aliasing... so keeping to inflate=1 for the time being.
        wc = np.append(wc,np.zeros((inflate-1)*wc.shape[0],dtype=np.int64)) # adding zeros to the end of the transfored vector
        ws = np.append(ws,np.zeros((inflate-1)*ws.shape[0],dtype=np.int64)) # adding zeros to the end of the transfored vector
        wy = np.append(wy,np.zeros((inflate-1)*wy.shape[0],dtype=np.int64)) # adding zeros to the end of the transfored vector
    wc[:s.shape[0]] *= np.arange(s.shape[0],dtype=np.int64) # producing the transform of the derivative
    ws[:s.shape[0]] *= np.arange(s.shape[0],dtype=np.int64) # producing the transform of the derivative
    dsc = (dst(wc,type=3)[:inflate*sz]//(4*sz**2)).astype(np.int64)
    dcs = (dct(ws,type=3)[:inflate*sz]//(4*sz**2)).astype(np.int64)
    dy = (dcs-dsc)
    y = (dct(wy,type=3,axis=0)[:inflate*sz]//(4*sz)).astype(np.int64)
    result = dy*tanhInt(-y,bits=8)   # constructing the logic waveform 
    return result

def dctLogic(s,inflate=1,nrolloff=128):
    result = np.zeros(s.shape,dtype=np.float32)
    if nrolloff>winsz:
        print('rolloff larger than windowed signal vec')
        return result
    if nrolloff!=0:
        print('rolloff is non-zero... dont bother with that')
        return result
    
    rolloff_vec = 0.5*(1.+np.cos(np.arange(nrolloff,dtype=float)*np.pi/float(nrolloff)))
    sz_roll = rolloff_vec.shape[0] 
    sz = s.shape[0]
    Yc = dct(np.append(s,np.flip(s,axis=0)),type=2)
    Ys = dst(np.append(s,np.flip(s,axis=0)),type=2)
    Yc[-sz_roll:] *= rolloff_vec
    Ys[-sz_roll:] *= rolloff_vec
    if inflate>1: # inflating seems to increase the aliasing... so keeping to inflate=1 for the time being.
        Yc = np.append(Yc,np.zeros((inflate-1)*Yc.shape[0])) # adding zeros to the end of the transfored vector
        Ys = np.append(Ys,np.zeros((inflate-1)*Ys.shape[0])) # adding zeros to the end of the transfored vector
    DYc = np.copy(Yc)
    DYs = np.copy(Ys)
    DYc[:s.shape[0]] *= np.arange(s.shape[0],dtype=float)/s.shape[0] # producing the transform of the derivative
    DYs[:s.shape[0]] *= np.arange(s.shape[0],dtype=float)/s.shape[0] # producing the transform of the derivative
    dys = dst(DYc,type=3)[:inflate*sz]/(4*sz**2)
    dyc = dct(DYs,type=3)[:inflate*sz]/(4*sz**2)
    dy = (dyc-dys)
    dy[:-1] /= np.cos(np.pi*np.arange(inflate*sz)/2.)[:-1]
    y = dct(Yc,type=3)[:inflate*sz]
    result = y*dy   # constructing the sig*deriv waveform 
    return result



def scanedges_simple(d,minthresh,expand=1):
    tofs = []
    slopes = []
    sz = d.shape[0]
    i = 10
    while i < sz-10:
        while d[i] > minthresh:
            i += 1
            if i==sz-10: return tofs,slopes,len(tofs)
        while i<sz-10 and d[i]<0:
            i += 1
        stop = i
        x0 = stop - 1./float(d[stop]-d[stop-1])*d[stop] 
        i += 1
        v = expand*float(x0)
        tofs += [np.int32(v) + np.int32(rng.random()<v%1)] 
        slopes += [d[stop]-d[stop-1]] ## scaling to reign in the obscene derivatives... probably shoul;d be scaling d here instead
    return tofs,slopes,len(tofs)

def scanedges(d,minthresh,expand=4):
    tofs = []
    slopes = []
    sz = d.shape[0]
    newtloops = 6
    order = 3 # this should stay fixed, since the logic zeros crossings really are cubic polys
    i = 10
    while i < sz-10:
        while d[i] > minthresh:
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
        X0 = np.array([np.power(x0,int(i)) for k in range(order+1)])
        #slopes += [np.int64(theta.dot([i*X0[(i+1)%(order+1)] for i in range(order+1)]))]
        slopes += [float((theta[1]+x0*theta[2])/2**18)] ## scaling to reign in the obscene derivatives... probably shoul;d be scaling d here instead
    return tofs,slopes,len(tofs)

class Port:
    # Note that t0s are aligned with 'prompt' in the digitizer logic signal
    # Don't forget to multiply by inflate, also, these look to jitter by up to 1 ns
    # hard coded the x4 scale-up for the sake of filling int16 dynamic range with the 12bit vls data and finer adjustment with adc offset correction

    def __init__(self,portnum,hsd,t0=0,nadcs=4,baselim=1000,logicthresh=-2400,scale=1,inflate=1,expand=1,nrolloff=256): # exand is for sake of Newton-Raphson
        self.portnum = portnum
        self.hsd = hsd
        self.t0 = t0
        self.nadcs = nadcs
        self.baselim = baselim
        self.logicthresh = logicthresh
        self.initState = True
        self.scale = scale
        self.inflate = inflate
        self.expand = expand
        self.nrolloff = nrolloff
        self.sz = 0
        self.tofs = []
        self.slopes = []
        self.addresses = []
        self.nedges = []
        self.waves = {}
        self.logics = {}
        self.shot = int(0)

    def addsample(self,w,l):
        eventnum = len(self.addresses)
        if eventnum<100:
            if eventnum%10<10: 
                self.waves.update( {'shot_%i'%eventnum:np.copy(w)} )
                self.logics.update( {'shot_%i'%eventnum:np.copy(l)} )
        elif eventnum<1000:
            if eventnum%100<10: 
                self.waves.update( {'shot_%i'%eventnum:np.copy(w)} )
                self.logics.update( {'shot_%i'%eventnum:np.copy(l)} )
        elif eventnum<10000:
            if eventnum%1000<10: 
                self.waves.update( {'shot_%i'%eventnum:np.copy(w)} )
                self.logics.update( {'shot_%i'%eventnum:np.copy(l)} )
        else:
            if eventnum%10000<10: 
                self.waves.update( {'shot_%i'%eventnum:np.copy(w)} )
                self.logics.update( {'shot_%i'%eventnum:np.copy(l)} )
        return self

    def process_list(self,ss):
        e = []
        de = []
        ne = 0
        if type(ss[0]) == type(None):
            e = []
            de = []
            ne = 0
        else:
            for s in ss:
                for adc in range(self.nadcs):
                    b = np.mean(s[adc:self.baselim+adc:self.nadcs])
                    s[adc::self.nadcs] = (s[adc::self.nadcs] * self.scale ) - int(self.scale*b)
                logic = dctLogicInt(s,inflate=self.inflate,nrolloff=self.nrolloff)
                edene = scanedges_simple(logic,self.logicthresh,self.expand)
                if edene[2]>0:
                    e += edene[0]
                    de += edene[1]
                    ne += edene[2]

        if self.initState:
            self.sz = ss[0].shape[0]*self.inflate*self.expand
            self.tofs = [0]
            if ne<1:
                self.addresses = [int(0)]
                self.nedges = [int(0)]
            else:
                self.addresses = [int(1)]
                self.nedges = [int(ne)]
                self.tofs += e
                self.slopes += de
        else:
            if ne<1:
                self.addresses += [int(0)]
                self.nedges += [int(0)]
            else:
                self.addresses += [int(len(self.tofs))]
                self.nedges += [int(ne)]
                self.tofs += e
                self.slopes += de
        return self

    def process(self,s):
        e = []
        de = []
        ne = 0
        if type(s) == type(None):
            e = []
            de = []
            ne = 0
        else:
            for adc in range(self.nadcs):
                b = np.mean(s[adc:self.baselim+adc:self.nadcs])
                s[adc::self.nadcs] = (s[adc::self.nadcs] * self.scale) - int(self.scale*b)
            logic = dctLogicInt(s,inflate=self.inflate,nrolloff=self.nrolloff) #produce the "logic vector"
            e,de,ne = scanedges_simple(logic,self.logicthresh,self.expand) # scan the logic vector for hits
            self.addsample(s,logic)

        if self.initState:
            self.sz = s.shape[0]*self.inflate*self.expand
            self.tofs = [0]
            if ne<1:
                self.addresses = [int(0)]
                self.nedges = [int(0)]
                self.tofs += []
                self.slopes += []
            else:
                self.addresses = [int(1)]
                self.nedges = [int(ne)]
                self.tofs += e
                self.slopes += de
        else:
            if ne<1:
                self.addresses += [int(0)]
                self.nedges += [int(0)]
                self.tofs += []
                self.slopes += []
            else:
                self.addresses += [int(len(self.tofs))]
                self.nedges += [int(ne)]
                self.tofs += e
                self.slopes += de
        return self

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



'''
def dctLogic(s,inflate=1,nrolloff=128):
    rolloff_vec = 0.5*(1.+np.cos(np.arange(nrolloff,dtype=float)*np.pi/float(nrolloff)))
    sz_roll = rolloff_vec.shape[0] 
    sz = s.shape[0]
    wave = np.append(s,np.flip(s,axis=0))
    WAVE = dct(wave,type=2)
    #WAVE = dct(wave)
    #WAVE = rollon(WAVE,4)
    WAVE[-sz_roll:] *= rolloff_vec
    if inflate>1: # inflating seems to increase the aliasing... so keeping to inflate=1 for the time being.
        WAVE = np.append(WAVE,np.zeros((inflate-1)*WAVE.shape[0])) # adding zeros to the end of the transfored vector
    DWAVE = np.copy(WAVE) # preparing to also make a derivative
    DWAVE[:s.shape[0]] *= np.arange(s.shape[0],dtype=float)/s.shape[0] # producing the transform of the derivative
    return dct(WAVE,type=3)[:inflate*sz]*dct(DWAVE,type=4)[:inflate*sz]/(4*sz**2) # constructing the sig*deriv waveform 
    #return dct(WAVE,type=3)[:inflate*sz]*dst(DWAVE,type=3)[:inflate*sz]/(4*sz**2) # constructing the sig*deriv waveform 
'''
