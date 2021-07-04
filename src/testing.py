#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.2.5/bin/python3

import psana
import numpy as np

def PWRspectrum(wv):
	return np.power(abs(np.fft.fft(wv).real),int(2))

def main():
	ds = psana.DataSource(exp='tmoc00118',run=479)
	for run in ds.runs():
		print(run.detnames)
		hsd = run.Detector('hsd')
		printkeys = True
		eventnum = 0
		for evt in run.events():
			if eventnum > 100:
				break

			if printkeys == True:
				print(hsd.raw.waveforms(evt).keys())
				wv0 = hsd.raw.waveforms(evt)[3][0]	
				pwr0 = np.zeros(shape=wv0.shape,dtype=float)
				pwr1 = np.zeros(shape=wv0.shape,dtype=float)
				pwr2 = np.zeros(shape=wv0.shape,dtype=float)
				pwr4 = np.zeros(shape=wv0.shape,dtype=float)
				pwr5 = np.zeros(shape=wv0.shape,dtype=float)
				pwr12 = np.zeros(shape=wv0.shape,dtype=float)
				pwr13 = np.zeros(shape=wv0.shape,dtype=float)
				pwr14 = np.zeros(shape=wv0.shape,dtype=float)
				pwr15 = np.zeros(shape=wv0.shape,dtype=float)
				pwr16 = np.zeros(shape=wv0.shape,dtype=float)

			wv0 = hsd.raw.waveforms(evt)[3][0]	
			#&0,3,TMO:MPOD:01:M2:C0,TMO:MPOD:01:M9:C0,TMO:MPOD:01:M5:C0,TMO:MPOD:01:M0:C0,TMO:MPOD:01:M1:C0,TMO:MPOD:01:M6:C0,TMO:MPOD:01:M7:C0,TMO:MPOD:01:M3:C0,
			wv1 = hsd.raw.waveforms(evt)[9][0]	
			#&1,9,TMO:MPOD:01:M2:C1,TMO:MPOD:01:M9:C1,TMO:MPOD:01:M5:C1,TMO:MPOD:01:M0:C1,TMO:MPOD:01:M1:C1,TMO:MPOD:01:M6:C1,TMO:MPOD:01:M7:C1,TMO:MPOD:01:M3:C1,
			wv2 = hsd.raw.waveforms(evt)[11][0]	
			#&2,11,TMO:MPOD:01:M2:C2,TMO:MPOD:01:M9:C2,TMO:MPOD:01:M5:C2,NA,NA,NA,NA,NA,
			wv4 = hsd.raw.waveforms(evt)[10][0]	
			#&4,10,TMO:MPOD:01:M2:C4,TMO:MPOD:01:M9:C4,TMO:MPOD:01:M5:C4,TMO:MPOD:01:M0:C4,TMO:MPOD:01:M1:C4,TMO:MPOD:01:M6:C4,TMO:MPOD:01:M7:C4,TMO:MPOD:01:M3:C4,
			wv5 = hsd.raw.waveforms(evt)[12][0]	
			#&5,12,TMO:MPOD:01:M2:C5,TMO:MPOD:01:M9:C5,TMO:MPOD:01:M5:C5,TMO:MPOD:01:M0:C5,TMO:MPOD:01:M1:C5,TMO:MPOD:01:M6:C5,TMO:MPOD:01:M7:C5,TMO:MPOD:01:M3:C5,
			wv12 = hsd.raw.waveforms(evt)[5][0]	
			#&12,5,TMO:MPOD:01:M2:C12,TMO:MPOD:01:M9:C12,TMO:MPOD:01:M5:C12,TMO:MPOD:01:M0:C12,TMO:MPOD:01:M1:C12,TMO:MPOD:01:M6:C12,TMO:MPOD:01:M7:C12,TMO:MPOD:01:M3:C12,
			wv13 = hsd.raw.waveforms(evt)[6][0]	
			#&13,6,TMO:MPOD:01:M2:C13,TMO:MPOD:01:M9:C13,TMO:MPOD:01:M5:C13,TMO:MPOD:01:M0:C13,TMO:MPOD:01:M1:C13,TMO:MPOD:01:M6:C13,TMO:MPOD:01:M7:C13,TMO:MPOD:01:M3:C13,
			wv14 = hsd.raw.waveforms(evt)[8][0]	
			#&14,8,TMO:MPOD:01:M2:C14,TMO:MPOD:01:M9:C14,TMO:MPOD:01:M5:C14,TMO:MPOD:01:M0:C14,TMO:MPOD:01:M1:C14,TMO:MPOD:01:M6:C14,TMO:MPOD:01:M7:C14,TMO:MPOD:01:M3:C14,
			wv15 = hsd.raw.waveforms(evt)[2][0]	
			#&15,2,TMO:MPOD:01:M2:C15,TMO:MPOD:01:M9:C15,TMO:MPOD:01:M5:C15,TMO:MPOD:01:M0:C15,TMO:MPOD:01:M1:C15,TMO:MPOD:01:M6:C15,TMO:MPOD:01:M7:C15,TMO:MPOD:01:M3:C15,
			wv16 = hsd.raw.waveforms(evt)[13][0]	
			#&16,13,TMO:MPOD:01:M2:C16,TMO:MPOD:01:M9:C6,TMO:MPOD:01:M5:C6,NA,NA,NA,NA,NA,

			pwr0 += PWRspectrum(wv0)
			pwr1 += PWRspectrum(wv1)
			pwr2 += PWRspectrum(wv2)
			pwr4 += PWRspectrum(wv4)
			pwr5 += PWRspectrum(wv5)
			pwr12 += PWRspectrum(wv12)
			pwr13 += PWRspectrum(wv13)
			pwr14 += PWRspectrum(wv14)
			pwr15 += PWRspectrum(wv15)
			pwr16 += PWRspectrum(wv16)
			eventnum += 1

		np.savetxt('pwerspect.dat',np.column_stack((pwr0,pwr1,pwr2,pwr4,pwr5,pwr12,pwr13,pwr14,pwr15,pwr16)))

	print("Hello")
	return

if __name__ == '__main__':
	main()
