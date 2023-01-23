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

def main():
	scratchdir = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output'
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
		#h5py.File('%s/waveforms.%s.run%i'%()) #HERE HERE HERE HERE HERE
		print(run.detnames)
		eventnum = 0
		runhsd=True
		runvls=False
		hsd = run.Detector('hsd')
		vls = run.Detector('andor')
		wv = {}
		wv_dct = {}
		wv_dct_deriv = {}
		wv_logic = {}
		chans = {0:3,1:9,2:11,4:10,5:12,12:5,13:6,14:8,15:2,16:13}

		init = True
		sz = 0
		inflate = 4
		for evt in run.events():
			if eventnum > nshots:
				break
			if init == True:
				sz = hsd.raw.waveforms(evt)[ chans[0] ][0].shape[0]
				for key in chans:
					wv[key] = np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=int)
					wave = np.append(np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=int),np.flip(np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=int),axis=0))
					WAVE = dct(wave)
					#WAVE[(3*sz)//4:] = 0
					#WAVE[0] = 0;WAVE[1] *= .1;WAVE[2] *= .25;WAVE[3] *= .5;WAVE[4] *= .75;WAVE[5] *= .9;
					WAVE = rollon(WAVE,10)
					WAVE = np.append(WAVE,np.zeros((inflate-1)*WAVE.shape[0]))
					DWAVE = np.copy(WAVE)
					DWAVE *= np.array(WAVE.shape[0],dtype=float)/WAVE.shape[0]
					wv_dct[key] = idct(WAVE)[:inflate*sz]/(2*sz)
					wv_dct_deriv[key] = idst(DWAVE)[:inflate*sz]/(2*sz)
					wv_logic[key] = (wv_dct[key] * wv_dct_deriv[key])
				init = False
			else:
				for key in chans:
					wv[key] = np.column_stack((wv[key],np.array( hsd.raw.waveforms(evt)[ chans[key] ][0], dtype=int) ) )
					wave = np.append(np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=int),np.flip(np.array(hsd.raw.waveforms(evt)[ chans[key] ][0] , dtype=int),axis=0))
					WAVE = dct(wave)
					#WAVE[(3*sz)//4:] = 0
					#WAVE[0] = 0;WAVE[1] *= .1;WAVE[2] *= .25;WAVE[3] *= .5;WAVE[4] *= .75;WAVE[5] *= .9;
					WAVE = rollon(WAVE,10)
					WAVE = np.append(WAVE,np.zeros((inflate-1)*WAVE.shape[0]))
					DWAVE = np.copy(WAVE)
					DWAVE *= np.array(WAVE.shape[0],dtype=float)/WAVE.shape[0]
					wv_dct[key] = np.column_stack((
								wv_dct[key],
								idct(WAVE)[:inflate*sz]/(2*sz)
								))
					wv_dct_deriv[key] = np.column_stack((
								wv_dct_deriv[key],
								idst(DWAVE)[:inflate*sz]/(2*sz)
								))
					wv_logic[key] = np.column_stack((
								wv_logic[key],
								(wv_dct[key][:,-1]*wv_dct_deriv[key][:,-1])
								))
				
				if eventnum%5<1: print(eventnum) 
			
			eventnum += 1

		for key in chans:
			print(wv[key].shape)
			headstring = 'port_%i'%(key)
			np.savetxt('%s/waveforms.%s.%i.%i.dat'%(scratchdir,expname,runnum,key),wv[key],fmt='%i',header=headstring)
			#np.savetxt('%s/waveforms_dct.%s.%i.%i.dat'%(scratchdir,expname,runnum,key),wv_dct[key],fmt='%i',header=headstring)
			#np.savetxt('%s/waveforms_dct_deriv.%s.%i.%i.dat'%(scratchdir,expname,runnum,key),wv_dct_deriv[key],fmt='%i',header=headstring)
			np.savetxt('%s/waveforms_logic.%s.%i.%i.dat'%(scratchdir,expname,runnum,key),wv_logic[key],fmt='%i',header=headstring)

	print("Hello, I'm done now!")
	return

if __name__ == '__main__':
	main()
