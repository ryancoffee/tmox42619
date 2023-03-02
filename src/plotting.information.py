#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import re

plotting = False
plottingchans = True

def shannon(s):
    ps = np.copy(s).astype(float)
    ps /= np.sum(s).astype(type(ps))
    inds = np.where(ps>0)
    return np.sum(-ps[inds]*np.log2(ps[inds]))

def main():
    
    if plotting:
        with h5py.File(sys.argv[1],'r') as f:
            portkeys = [k for k in f.keys() if (re.search('port',k) and not re.search('_16',k) and not re.search('_2',k))] # keeping the bare MCP ports 2 and 16 here
            for k in portkeys:
                if k != 'port_12':
                    continue
                print(k,f[k].keys())
                d = f[k]['hist'][()]
                s = np.sum(d,axis=1)
                w = (f[k]['qbins'][()][1:]-f[k]['qbins'][()][:-1])/6./8.
                b = 0.5*(f[k]['qbins'][()][1:]+f[k]['qbins'][()][:-1])/6./8.
                fig = plt.figure(figsize = (8,4))
                axq = fig.add_subplot(1,2,1)
                axs = fig.add_subplot(1,2,2)
                q = np.arange(len(s))
                axq.step(q,s)
                axq.plot(q,np.max(s)/np.max(2./w)*1./w)
                axq.set_title('no bin info')
                #axq.vlines(np.arange(len(s)),np.zeros(len(s)),s/w)
                axq.set_xlabel('quant bins')
                axq.set_ylabel('hits/bin')
                axq.legend(['counts','1/binwidth'])
                axs.set_title('with bin info')
                axs.step(b-b[0],s/w)
                axs.plot(b-b[0],np.max(s/w)/np.max(2./w)/w)
                #axs.vlines(b-b[0],np.zeros(len(s)),s/w)
                axs.set_xlabel('ToF [ns]')
                axs.set_ylabel('hits/ns')
                axs.set_xlim(0,250)
                axs.yaxis.set_label_position('right')
                axs.yaxis.tick_right()
                axs.legend(['counts/width','1/binwidth'])
                #plt.savefig('%s.png'%(sys.argv[1]))
                plt.show()
            print('%s\tShannon entropy = %.3f'%(k,shannon(s)))
            print('%s\tShannon entropy (1/widths) = %.3f'%(k,shannon(1/w)))
            print('%s\tRatio of entropy (quantized/binwidths) = %.3f'%(k,(shannon(s)/shannon(1/w))))

    if plottingchans:
        portorder = {'port_4':0,'port_12':1,'port_5':2,'port_13':3,'port_14':4,'port_15':5,'port_1':6,'port_0':7}
        fnames = sys.argv[1:]
        m = re.search('(^.*)/h5files.*\.h5',fnames[0])
        if not m:
            print('failed fname match')
            return
        path = m.group(1)
        nports = len(portorder.keys())
        nbins = 1<<8
        ncycles = 1<<3
        sumn = 1<<3
        with h5py.File(fnames[0],'r') as f:
            nbins = f['port_12']['hist'][()].shape[0]
        quantHist = np.zeros((nbins,ncycles*len(fnames)*nports),dtype=float)
        for i,fname in enumerate(fnames):
            with h5py.File(fname,'r') as f:
                for k,a in portorder.items():
                    for n in range(ncycles):
                        d = f[k]['hist'][()]
                        #s = np.sum(d,axis=1)
                        s = np.sum(d[:,n*sumn:n*sumn+sumn],axis=1)
                        w = (f[k]['qbins'][()][1:]-f[k]['qbins'][()][:-1])/6./8.
                        quantHist[:,i*nports*len(fnames)+a*nports+n] = np.log2(1+s/w)

        fig = plt.figure(figsize = (10,5))
        axq = fig.add_subplot(1,1,1)
        im = axq.pcolor(quantHist)
        axq.set_xlabel('time (or photon energy... or what?)')
        axq.set_ylabel('qbins [Who knows what the bins mean?]')
        fig.colorbar(im,ax=axq,label = 'log of values/secret widths')
        plt.savefig('%s/figs/mergingNeonRuns.png'%path)
        plt.show()
    return

if __name__ == '__main__':
    if len(sys.argv)>1:
        main()
    else:
        print('syntax: ./src/plotting.information.py <quantized input file>')

