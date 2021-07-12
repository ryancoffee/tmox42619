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

def PWRspectrum(wv):
	return np.power(abs(np.fft.fft(wv).real),int(2))

def rollon(vec,n):
	vec[:int(n)] = vec[:int(n)]*np.arange(int(n),dtype=float)/float(n)
	return vec


def mypoly(x,order=4):
	result = np.ones((x.shape[0],order+1),dtype=float)
	result[:,1] = x.copy()
	if order < 2:
		return result
	for p in range(2,order+1):
		result[:,p] = np.power(result[:,1],int(p))
	return result

def processEbeam(thisl3,l3,initState):
	if initState:
		l3 = [np.uint16(thisl3)]
	else:
		l3 += [np.uint16(thisl3)]
	return l3

def processVls(evt,vls,vlswv,v,vc,vs,vsize,initState):
	num = np.sum(np.array([i*vlswv[i] for i in range(len(vlswv))]))
	den = np.sum(vlswv)

	if initState:
		v = [vlswv.astype(np.int16)]
		vsize = len(v)
		vc = [np.uint16(num/den)]
		vs = [np.uint64(den)]
	else:
		v += [vlswv.astype(np.int16)]
		vc += [np.uint16(num/den)]
		vs += [np.uint64(den)]
	return v,vc,vs,vsize


def dctLogic(s,inflate=4):
	sz = s.shape[0]
	wave = np.append(s,np.flip(s,axis=0))
	WAVE = dct(wave)
	WAVE = rollon(WAVE,10)
	WAVE = np.append(WAVE,np.zeros((inflate-1)*WAVE.shape[0]))
	DWAVE = np.copy(WAVE)
	DWAVE[:s.shape[0]] *= np.arange(s.shape[0],dtype=float)/s.shape[0]
	return idct(WAVE)[:inflate*sz]*idst(DWAVE)[:inflate*sz]/(4*sz**2) # constructing the sig*deriv waveform 

def scanedges(d,minthresh):
	tofs = []
	sz = d.shape[0]
	newtloops = 3
	order = 3
	i = 1
	while i < sz-10:
		while d[i] > minthresh:
			i += 1
			if i==sz-10: return tofs,len(tofs)
		while i<sz-10 and d[i]<d[i-1]:
			i += 1
		start = i
		i += 1
		while i<sz-10 and d[i]>d[i-1]:
			i += 1
		stop = i
		if stop-start<4:
			continue
		x = np.arange(stop-start,dtype=float)
		y = d[start:stop]
		x0 = float(stop)/2.
		y -= (y[0]+y[-1])/2.
		theta = np.linalg.pinv( mypoly(np.array(x).astype(float),order=order) ).dot(np.array(y).astype(float))
		for j in range(newtloops): # 3 rounds of Newton-Raphson
			X0 = np.array([np.power(x0,int(i)) for i in range(order+1)])
			x0 -= theta.dot(X0)/theta.dot([i*X0[(i+1)%(order+1)] for i in range(order+1)]) # this seems like maybe it should be wrong
		tofs += [start + x0]
	return tofs,len(tofs)

def processPort(key,s,tofs,nedges,addresses,initState):

	if type(s) == type(None):
		e = []
		ne = 0
	else:
		nadcs = 4
		for adc in range(nadcs):
			base = np.mean(s[adc:baselim:nadcs])
			s[adc::4] -= np.int16(base) # this now needs to be a signed int
		logic = dctLogic(s,inflate)
		e,ne = scanedges(logic,logicthresh[key]) # logic waveform in and edges out
	if initState: 
		sz[key] = s.shape[0]
		tofs[key] = [0] # setting the 0'th address of tofs[key] to catch all addresses for nedges == 0 case
		if ne<1:
			addresses[key] = [int(0)]
			nedges[key] = [int(0)]
			tofs[key] += [] 
		else:
			addresses[key] = [int(1)] 
			nedges[key] = [int(ne)]
			tofs[key] += e
	else:
		if ne<1:
			addresses[key] += [int(0)]
			nedges[key] += [int(0)]
		else:
			addresses[key] += [int(len(tofs[key]))] 
			nedges[key] += [int(ne)]
			tofs[key] += e
	return tofs,nedges,addresses


def main():
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

	ds = psana.DataSource(exp=expname,run=runnum)

	for run in ds.runs():
		#np.savetxt('%s/waveforms.%s.%i.%i.dat'%(scratchdir,expname,runnum,key),wv[key],fmt='%i',header=headstring)
		print(run.detnames)
		eventnum = 0
		runhsd=True
		runvls=False
		hsd = run.Detector('hsd')
		vls = run.Detector('andor')
		ebeam = run.Detector('ebeam')
		wv = {}
		wv_logic = {}
		v = []
		vc = []
		vs = []
		l3 = []
		tofs = {}
		addresses = {}
		nedges = {}
		sz = {}
		# Note that t0s are aligned with 'prompt' in the digitizer logic signal
		# Don't forget to multiply by inflate, also, these look to jitter by up to 1 ns
		chans = {0:3,1:9,2:11,4:10,5:12,12:5,13:6,14:8,15:2,16:13} # HSD to port number
		t0s = {0:4585,1:4206,2:4166,4:4055,5:4139,12:4139,13:4133,14:4185,15:4460,16:4096}
		logicthresh = {0:-8000, 1:-8000, 2:-400, 4:-8000, 5:-8000, 12:-8000, 13:-8000, 14:-8000, 15:-8000, 16:-8000}
		# hard coded the x4 scale-up for the sake of filling int16 dynamic range with the 12bit vls data and finer adjustment with adc offset correction


		init = True
		vsize = 0
		inflate = 4
		baselim = 1000 
		for evt in run.events():
			if eventnum > nshots:
				break
			try:
				vlswv = np.squeeze(vls.raw.value(evt))
				vlswv = vlswv-int(np.mean(vlswv[1900:]))
				if np.max(vlswv<300): 
					if eventnum%50<1: 
						print(eventnum)
					eventnum += 1
					continue
			except:
				print('skip per vls')
				continue

			try:
				thisl3 = ebeam.raw.ebeamL3Energy(evt)
				thisl3 += 0.5
			except:
				print('skip per l3')
				continue


			''' Ebeam specific section '''
			l3 = processEbeam(thisl3,l3,initState=init)	

			''' VLS specific section, do this first to slice only good shots '''

			v,vc,vs,vsize = processVls(vlswv,v,vc,vs,vsize,initState=init)

			''' HSD-Abaco section '''
			for key in chans.keys():
				# hard coding the scale inflation for accounting the 4 different ADC offsets.

				# wrap in callable method from here to...
				# inputs would be int16 waveform, and the inflation factor also include the scaling by 4 to use 14 bits depth instead of only 12
				s = 4*np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=np.int16) 

				tofs,nedges,addresses = processPort(key,s,tofs,nedges,addresses,initState=init)

			if init and len(v)>0: init = False

			if eventnum%50==0: 
				print(eventnum)
			eventnum += 1

		f = h5py.File('%s/hits.%s.run%i.h5'%(scratchdir,expname,runnum),'w') 
                # use f.create_group('port_%i'%i,portnum)
		#_ = [print(key,chans[key]) for key in chans.keys()]
		print(chans.keys(),tofs.keys())
		for key in chans.keys():
			g = f.create_group('port_%i'%(key))
			g.create_dataset('tofs',data=tofs[key],dtype=np.uint32) 
			g.create_dataset('addresses',data=addresses[key],dtype=np.uint64)
			g.create_dataset('nedges',data=nedges[key],dtype=np.uint16)
			g.attrs.create('inflate',data=inflate,dtype=np.uint8)
			g.attrs.create('t0',data=t0s[key]*inflate,dtype=np.uint8)
			g.attrs.create('hsd',data=chans[key],dtype=np.uint8)
			g.attrs.create('size',data=sz[key]*inflate,dtype=np.uint8)
		grpvls = f.create_group('vls')
		grpvls.create_dataset('data',data=v,dtype=np.int16)
		grpvls.create_dataset('centroids',data=vc,dtype=np.int16)
		grpvls.create_dataset('sum',data=vs,dtype=np.uint64)
		grpvls.attrs.create('size',data=vsize,dtype=np.int32)
		grpebeam = f.create_group('ebeam')
		grpebeam.create_dataset('l3energy',data=l3,dtype=np.uint16)
		f.close()

	print("Hello, I'm done now!")
	return

if __name__ == '__main__':
	main()
