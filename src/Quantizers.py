#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.2.5/bin/python3

import numpy as np
import h5py
import sys
import re
import matplotlib.pyplot as plt
from typing import List

class Quantizer:
    def __init__(self,style='nonuniform',nbins=1024):
        self.style:str = style
        self.nbins:np.uint32 = nbins
        self.qbins:List[float] = []

    def setbins(self,data):
        if self.style=='nonuniform':
            ubins = np.arange(np.min(data),np.max(data)+1)
            csum = np.cumsum( np.histogram(data,bins=ubins)[0] )
            yb = np.arange(0,csum[-1],step=np.float(csum[-1])/(self.nbins+1))
            self.qbins = np.interp(yb,csum,(ubins[:-1]+ubins[1:])/2.)
        else:
            self.qbins = np.arange(np.min(data),np.max(data)+1)
        return self

    def histogram(self,data):
        return np.histogram(data,bins=self.qbins)[0]
    def bincenters(self):
        return (self.qbins[:-1] + self.qbins[1:])/2.0
    def binedges(self):
        return self.qbins
    def binwidths(self):
        return self.qbins[1:]-self.qbins[:-1]
    def saveH5(self,fname):
        with h5py.File(fname,'w') as f:
            ds = f.create_dataset('qbins',data=self.qbins,dtype=float)
            ds.attrs.create('nbins',self.nbins,dtype=np.uint32)
            ds.attrs.create('style',self.style,dtype=str)
        return self
    def copybins(self,bins):
        self.qbins = bins
        return self

    @classmethod
    def loadH5(cls,fname):
        with h5py.File(fname,'r') as f:
            q = cls(style=f['qbins'].attrs['style'],nbins=f['qbins'].attrs['nbins'])
            q.copybins(f['qbins'][()])
        return q

def main():
    if len(sys.argv)<2:
        print('syntax: Quantizer.py <nbins> <fname> <opt plotting? [Tt]rue>')
        return
    plotting = False
    if len(sys.argv)>2:
        if (sys.argv[-1] == 'True' or sys.argv[-1]=='true'):
            plotting = True

    fname = sys.argv[2]
    data = {} 
    quants = {}
    with h5py.File(fname,'r') as f:
        portkeys = [k for k in f.keys() if re.search('port',k)]
        for k in portkeys:
            quants[k] = Quantizer(style='nonuniform',nbins=np.uint32(sys.argv[1]))
            quants[k].setbins(f[k]['tofs'][()])
            print(len(quants[k].bincenters()))
            data[k] = quants[k].histogram(f[k]['tofs'][()])

    if plotting:
        for k in list(quants.keys())[:1]:
            #plt.plot(data[k],'.')
            #plt.plot(quants[k].bincenters(),data[k]/quants[k].binwidths(),'.')
            plt.stem(quants[k].bincenters(),1000./quants[k].binwidths(),linefmt='b-',markerfmt=' ')
        plt.show()
        
    return

if __name__ == '__main__':
    main()
