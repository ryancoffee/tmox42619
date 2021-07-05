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

def scanedges(d,minthresh,inflate=1):
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
		f = h5py.File('%s/hits.%s.run%i.h5'%(scratchdir,expname,runnum),'w') #HERE HERE HERE HERE HERE
		print(run.detnames)
		eventnum = 0
		runhsd=True
		runvls=False
		hsd = run.Detector('hsd')
		vls = run.Detector('andor')
		ebeam = run.Detector('ebeam')
		wv = {}
		#wv_dct = {}
		#wv_dct_deriv = {}
		wv_logic = {}
		v = []
		vc = []
		vs = []
		l3 = []
		tofs = {}
		addresses = {}
		nedges = {}
		chans = {0:3,1:9,2:11,4:10,5:12,12:5,13:6,14:8,15:2,16:13}
		# Note that t0s are aligned with 'prompt' in the digitizer logic signal
		# Don't forget to multiply by inflate, also, these look to jitter by up to 1 ns
		t0s = {0:4585,1:4206,2:4166,4:4055,5:4139,12:4139,13:4133,14:4185,15:4460,16:4096}

		logicthresh = {0:-2000, 1:-2000, 2:-100, 4:-2000, 5:-2000, 12:-2000, 13:-2000, 14:-2000, 15:-2000, 16:-2000}

		init = True
		sz = 0
		inflate = 4
		baselim = 1000
		for evt in run.events():
			if eventnum > nshots:
				break
			'''
			if eventnum < 22800:
				eventnum += 1
				if eventnum%500<1: 
					print(eventnum)
				continue
			'''

			''' VLS specific section, do this first to slice only good shots '''
			#print(vlswv.shape)
			try:
				vlswv = np.squeeze(vls.raw.value(evt))
				baseline = vlswv[1900:]#/148.
			except:
				print('skip per vls')
				continue

			thisl3 = ebeam.raw.ebeamL3Energy(evt)
			if type(thisl3) == type(None):
				print('skip per l3')
				continue
			else:
				thisl3 += 0.5
			bs = np.sum(baseline)/len(baseline)
			vlscrop = vlswv[950:1250]-bs
			num = np.sum(np.array([i*vlscrop[i] for i in range(len(vlscrop))]))
			den = np.sum(vlscrop)
			if type(den) == type(None):
				print('skip per vlscrop')
				continue
			if den<5000: 
				if eventnum%50<1: print(eventnum)
				eventnum += 1
				continue

			if init:
				v = [vlscrop]
				vc = [np.uint16(num/den)]
				vs = [np.uint64(den)]
				l3 = [np.uint16(thisl3)]
			else:
				v += [vlscrop]
				vc += [np.uint16(num/den)]
				vs += [np.uint64(den)]
				l3 += [np.uint16(thisl3)]


			''' HSD-Abaco section '''
			sz = hsd.raw.waveforms(evt)[ chans[0] ][0].shape[0]

			for key in chans:
				s = np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=np.int16)

				if type(s) == type(None):
					e = []
					ne = 0
				else:

					base = np.mean(s[:baselim])
					s -= np.int16(base)
					wave = np.append(s,np.flip(s,axis=0))
					WAVE = dct(wave)
					WAVE = rollon(WAVE,10)
					WAVE = np.append(WAVE,np.zeros((inflate-1)*WAVE.shape[0]))
					DWAVE = np.copy(WAVE)
					DWAVE[:s.shape[0]] *= np.arange(s.shape[0],dtype=float)/s.shape[0]
					logic = idct(WAVE)[:inflate*sz]*idst(DWAVE)[:inflate*sz]/(4*sz**2)
					e,ne = scanedges(logic,logicthresh[key],inflate)
				if init: 
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


			if init and len(v)>0: init = False

			if eventnum%50<1: print(eventnum); eventnum += 1

		for key in chans:
			f.create_dataset('port_%i_tofs'%(key),data=tofs[key],dtype=np.uint32)
			f.create_dataset('port_%i_addresses'%(key),data=addresses[key],dtype=np.uint64)
			f.create_dataset('port_%i_nedges'%(key),data=nedges[key],dtype=np.uint16)
		f.create_dataset('vls',data=v,dtype=np.int16)
		f.create_dataset('vls_centroids',data=vc,dtype=np.int16)
		f.create_dataset('vls_sum',data=vs,dtype=np.uint64)
		f.create_dataset('l3energy',data=l3,dtype=np.uint16)
		f.close()

	print("Hello, I'm done now!")
	return

if __name__ == '__main__':
	main()
