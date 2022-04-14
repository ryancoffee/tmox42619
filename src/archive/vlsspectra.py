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

def PWRspectrum(wv):
	return np.power(abs(np.fft.fft(wv).real),int(2))

def main():
	expname = 'tmox42619'
	runnum = 17
	if len(sys.argv)>2:
		expname = sys.argv[1]
		runnum = int(sys.argv[2])

	print('starting analysis exp %s for run %i'%(expname,int(runnum)))

	ds = psana.DataSource(exp=expname,run=runnum)

	for run in ds.runs():
		print(run.detnames)
		eventnum = 0
		runhsd=False
		runvls=True
		if runhsd:
			hsd = run.Detector('hsd')
			initaccum = True
			chans = {0:3,1:9,2:11,4:10,5:12,12:5,13:6,14:8,15:2,16:13}
			pwr = {}

			for evt in run.events():
				if eventnum > 10000:
					break

				if initaccum == True:
					print(hsd.raw.waveforms(evt).keys())
					tmp = hsd.raw.waveforms(evt)[3][0]	
					for key in chans:
						pwr.update({key:np.zeros(shape=tmp.shape,dtype=float)})
					initaccum = False

			
				for key in pwr:
					pwr[key] += PWRspectrum(hsd.raw.waveforms(evt)[chans[key]][0])

				eventnum += 1

			for key in pwr:
				scale = np.mean(pwr[key][5900:6000])
				pwr[key] /= float(scale)
			headstring = "port0\tport1\tport2\tport4\tport5\tport12\tport13\tport14\tport15\tport16"
			np.savetxt('pwrspect.%s.%i.noisenorm.dat'%(expname,runnum),np.column_stack((pwr[0],pwr[1],pwr[2],pwr[4],pwr[5],pwr[12],pwr[13],pwr[14],pwr[15],pwr[16])),header=headstring)
		if runvls:
			vls = run.Detector('andor')
			initaccum = True
			sp = {}
			for evt in run.events():
				if eventnum>1000:
					break
				sp.update({eventnum:vls.raw.value(evt)})
				eventnum += 1


			headstring = '#'
			out = []
			for key in sp.keys():
				if initaccum:
					out = sp[key]
					initaccum = False
				else:
					out = np.row_stack((out,sp[key]))

				headstring += '%i'%key
			np.savetxt('vlsspectrum.%s.%i.dat'%(expname,runnum),out.T,fmt='%i',header=headstring)


	print("Hello, I'm done now!")
	return

if __name__ == '__main__':
	main()
