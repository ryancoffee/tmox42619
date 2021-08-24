#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.2.5/bin/python3

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

import psana
import numpy as np
import sys
import h5py
from scipy.fftpack import dct,idct,idst

from utils import mypoly

def PWRspectrum(wv):
	return np.power(abs(np.fft.fft(wv).real),int(2))

def rollon(vec,n):
	vec[:int(n)] = vec[:int(n)]*np.arange(int(n),dtype=float)/float(n)
	return vec



###########################################
########### Class definitions #############
###########################################

class Ebeam:
	def __init__(self):
		self.l3 = []
		self.initState = True
		return 
	def process(self,l3in):
		if self.initState:
			self.l3 = [l3in]
		else:
			self.l3 += [np.uint16(l3in)]
		return self
	def set_initState(self,state):
		self.initState = state
		return self


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

	
def dctLogic(s,inflate=4):
	sz = s.shape[0]
	wave = np.append(s,np.flip(s,axis=0))
	WAVE = dct(wave)
	WAVE = rollon(WAVE,10)
	WAVE = np.append(WAVE,np.zeros((inflate-1)*WAVE.shape[0])) # adding zeros to the end of the transfored vector
	DWAVE = np.copy(WAVE) # preparing to also make a derivative
	DWAVE[:s.shape[0]] *= np.arange(s.shape[0],dtype=float)/s.shape[0] # producing the transform of the derivative
	return idct(WAVE)[:inflate*sz]*idst(DWAVE)[:inflate*sz]/(4*sz**2) # constructing the sig*deriv waveform 

def scanedges(d,minthresh,expand=4):
	tofs = []
	slopes = []
	sz = d.shape[0]
	newtloops = 3
	order = 3
	i = 1
	while i < sz-10:
		while d[i] > minthresh:
			i += 1
			if i==sz-10: return tofs,slopes,len(tofs)
		while i<sz-10 and d[i]<d[i-1]:
			i += 1
		start = i
		i += 1
		while i<sz-10 and d[i]>d[i-1]:
			i += 1
		stop = i
		if stop-start<4:
			continue
		x = expand*np.arange(stop-start,dtype=float) # set x to index values
		y = d[start:stop] # set y to vector values
		x0 = float(expand*stop)/2. # set x0 to halfway point
		y -= (y[0]+y[-1])/2. # subtract average (this gets rid of residual DC offsets)
		theta = np.linalg.pinv( mypoly(np.array(x).astype(float),order=order) ).dot(np.array(y).astype(float)) # fit a polynomial (order 3) to the points
		for j in range(newtloops): # 3 rounds of Newton-Raphson
			X0 = np.array([np.power(x0,int(i)) for i in range(order+1)])
			x0 -= theta.dot(X0)/theta.dot([i*X0[(i+1)%(order+1)] for i in range(order+1)]) # this seems like maybe it should be wrong
		tofs += [np.uint32(start + x0)] ###### CAREFUL  THIS IS NEW!  (4 am idea)... expand here is further subdividing the infated indices beacaus of Newton-Raphson.
		X0 = np.array([np.power(x0,int(i)) for i in range(order+1)])
		#slopes += [np.int32(theta.dot([i*X0[(i+1)%(order+1)] for i in range(order+1)]))]
		slopes += [np.int16((theta[1]+x0*theta[2])/2**18)] ## scaling to reign in the obscene derivatives... probably shoul;d be scaling d here instead
	return tofs,slopes,len(tofs)

class Port:
	# Note that t0s are aligned with 'prompt' in the digitizer logic signal
	# Don't forget to multiply by inflate, also, these look to jitter by up to 1 ns
	# hard coded the x4 scale-up for the sake of filling int16 dynamic range with the 12bit vls data and finer adjustment with adc offset correction

	def __init__(self,portnum,hsd,t0=0,nadcs=4,baselim=1000,logicthresh=-24000,slopethresh=500,scale=4,inflate=4,expand=4): # exand is for sake of Newton-Raphson
		self.portnum = portnum
		self.hsd = hsd
		self.t0 = t0
		self.nadcs = nadcs
		self.baselim = baselim
		self.logicthresh = logicthresh
		self.slopethresh = slopethresh
		self.initState = True
		self.scale = scale
		self.inflate = inflate
		self.expand = expand
		self.sz = 0
		self.tofs = []
		self.slopes = []
		self.addresses = []
		self.nedges = []
		self.waves = {}
		self.shot = int(0)

	def process(self,s):
		if type(s) == type(None):
			e = []
			de = []
			ne = 0
		else:
			s *= self.scale
			#print(s[:10])
			for adc in range(self.nadcs):
				b = np.mean(s[adc:self.baselim:self.nadcs])
				s[adc::self.nadcs] -= np.int16(b)
			logic = dctLogic(s,self.inflate) #produce the "logic vector"
			e,de,ne = scanedges(logic,self.logicthresh,self.expand) # scan the logic vector for hits
                        # the expand here is how much we subdivide the pixels in the already dct expanded digitizer steps (sake of Newton-Raphson root resolution)
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

def main():
        ############################################
        ###### Change this to your output dir ######
        ############################################
	scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output/h5files'
	expname = 'tmox42619'
	runnum = 62 
	nshots = 100
	if len(sys.argv)>2:
		expname = sys.argv[1]
		runnum = int(sys.argv[2])

	if len(sys.argv)>3:
		nshots = int(sys.argv[3])

	print('starting analysis exp %s for run %i'%(expname,int(runnum)))
	nr_expand = 4
	chans = {0:3,1:9,2:11,4:10,5:12,12:5,13:6,14:8,15:2,16:13} # HSD to port number:hsd
	logicthresh = {0:-8000, 1:-8000, 2:-400, 4:-8000, 5:-8000, 12:-8000, 13:-8000, 14:-8000, 15:-8000, 16:-8000}
	slopethresh = {0:500,1:500,2:300,4:150,5:500,12:500,13:500,14:500,15:500,16:300}
	#t0s = {0:109840,1:100456,2:99924,4:97180,5:99072,12:98580,13:98676,14:100348,15:106968,16:98028}
	t0s = {0:109830,1:100451,2:99810,4:97180,5:99071,12:98561,13:98657,14:100331,15:106956,16:97330}
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


	spect = Vls()
	ebunch = Ebeam()
	port = {} 
	scale = 4
	inflate = 4 # this determines the oversampling
	for key in logicthresh.keys():
		logicthresh[key] *= scale # inflating by factor of 4 since we are also scaling the waveforms by 4 in vertical to fill bit depth.

	for key in chans.keys():
		port[key] = Port(key,chans[key],t0=t0s[key],logicthresh=logicthresh[key],slopethresh=slopethresh[key],inflate=inflate,expand=nr_expand,scale=scale)

	ds = psana.DataSource(exp=expname,run=runnum)

	#for run in ds.runs():
	run = next(ds.runs())
		#np.savetxt('%s/waveforms.%s.%i.%i.dat'%(scratchdir,expname,runnum,key),wv[key],fmt='%i',header=headstring)
	for i in range(1):
		print(run.detnames)
		eventnum = 0
		runhsd=True
		runvls=False
		hsd = run.Detector('hsd')
		vls = run.Detector('andor')
		ebeam = run.Detector('ebeam')
		wv = {}
		wv_logic = {}
		v = [] # vls data matrix
		vc = [] # vls centroids vector
		vs = [] # vls sum is I think not used, maybe for normalization or used to be for integration and PDF sampling
		l3 = [] # e-beam l3 (linac 3) in GeV.


		init = True 
		vsize = 0

		for evt in run.events():
			if eventnum > nshots:
				break

			''' VLS specific section, do this first to slice only good shots '''
			try:
				vlswv = np.squeeze(vls.raw.value(evt))
				vlswv = vlswv-int(np.mean(vlswv[1900:])) # this subtracts baseline
				if np.max(vlswv)<300:  # too little amount of xrays
					print(eventnum,'skip per weak vls')
					eventnum += 1
					continue
				spect.process(vlswv)
				#spect.print_v()

			except:
				print(eventnum,'skip per vls')
				eventnum += 1
				continue

			''' Ebeam specific section '''
			try:
				thisl3 = ebeam.raw.ebeamL3Energy(evt)
				thisl3 += 0.5
				ebunch.process(thisl3)
			except:
				print(eventnum,'skip per l3')
				eventnum += 1
				continue


			''' HSD-Abaco section '''
			for key in chans.keys(): # here key means 'port number'
				s = np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=np.int16) 
				port[key].process(s)

			if init==True:
				init = False
				ebunch.set_initState(False)
				spect.set_initState(False)
				for key in chans.keys():
					port[key].set_initState(False)

			if eventnum%50==0: 
				print(eventnum)
			eventnum += 1

		f = h5py.File('%s/hits.%s.run%i.hires.h5'%(scratchdir,expname,runnum),'w') 
                # use f.create_group('port_%i'%i,portnum)
		#_ = [print(key,chans[key]) for key in chans.keys()]
		for key in chans.keys(): # remember key == port number
			g = f.create_group('port_%i'%(key))
			g.create_dataset('tofs',data=port[key].tofs,dtype=np.uint32) 
			g.create_dataset('slopes',data=port[key].slopes,dtype=np.int64) 
			g.create_dataset('addresses',data=port[key].addresses,dtype=np.uint64)
			g.create_dataset('nedges',data=port[key].nedges,dtype=np.uint16)
			g.attrs.create('inflate',data=port[key].inflate,dtype=np.uint8)
			g.attrs.create('expand',data=port[key].expand,dtype=np.uint8)
			g.attrs.create('t0',data=port[key].t0,dtype=np.uint32)
			g.attrs.create('slopethresh',data=port[key].slopethresh,dtype=np.uint16)
			g.attrs.create('hsd',data=port[key].hsd,dtype=np.uint8)
			g.attrs.create('size',data=port[key].sz*port[key].inflate,dtype=np.uint32) ### need to also multiply by expand #### HERE HERE HERE HERE
		grpvls = f.create_group('vls')
		grpvls.create_dataset('data',data=spect.v,dtype=np.int16)
		grpvls.create_dataset('centroids',data=spect.vc,dtype=np.int16)
		grpvls.create_dataset('sum',data=spect.vs,dtype=np.uint64)
		grpvls.attrs.create('size',data=spect.vsize,dtype=np.int32)
		grpebeam = f.create_group('ebeam')
		grpebeam.create_dataset('l3energy',data=ebunch.l3,dtype=np.uint16)
		f.close()

	print("Hello, I'm done now!")
	return

if __name__ == '__main__':
	main()
