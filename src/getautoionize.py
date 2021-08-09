#!/usr/bin/python3
import numpy as np
import sys
import re

def main():
    if len(sys.argv) <2:
        print('give me a file to process')
        return
    fname = sys.argv[1]
    m = re.search('(.+)/Neon_port_(\d+).(enhist_\w+).dat',fname)
    if m:
        print(m.group(3))
        path = m.group(1)
        port = m.group(2)
        nametail = m.group(3)

        data = np.loadtxt(fname)
        print(data.shape)

        sp875=np.sum(data[:,310:320],axis=1)
        sp880=np.sum(data[:,330:340],axis=1)
        outfile = '%s/Neon_875.port_%d.%s.spectrum.dat'%(path,int(port),nametail)
        np.savetxt(outfile,sp875,fmt='%i')
        outfile = '%s/Neon_880.port_%d.%s.spectrum.dat'%(path,int(port),nametail)
        np.savetxt(outfile,sp880,fmt='%i')


    return

if __name__ == '__main__':
    main()
