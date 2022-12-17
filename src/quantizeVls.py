#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3
import numpy as np
import sys
import h5py
import re
import matplotlib.pyplot as plt
from Quantizers import Quantizer

def main():
    if len(sys.argv)<3:
        print('syntax: quantizeHits.py <ntofbins> <nvlsbins> <fnames>')
        return

    donorm = False
    fnames = sys.argv[3:]
    ntofbins = np.uint32(sys.argv[1])
    nvlsbins = np.uint32(sys.argv[2])
    vlsoffset = 1000
    tofs = {} 
    addresses = {} 
    nedges = {} 
    quants = {}
    hist = {}
    vlscenters = []
    vlssums = []
    portkeys = []
    for fname in fnames:
        with h5py.File(fname,'r') as f:
            portkeys = [k for k in f.keys() if (re.search('port',k) and not re.search('_16',k) and not re.search('_2',k))]
            if len(quants.keys())==0:
                for k in portkeys:
                    quants[k] = Quantizer(style='nonuniform',nbins=ntofbins)
                    tofs[k] = list(f[k]['tofs'][()])
                    addresses[k] = list(f[k]['addresses'][()])
                    nedges[k] = list(f[k]['nedges'][()])
                    hist[k] = np.zeros((nvlsbins,ntofbins),dtype=float)
                vlscenters = list(f['vls']['centroids'][()])
                vlssums = list(f['vls']['sum'][()])
            else:
                for k in portkeys:
                    offsetTofs = len(tofs[k])
                    addresses[k] += [int(offsetTofs+a) for a in f[k]['addresses'][()]]
                    nedges[k] += list(f[k]['nedges'][()])
                    tofs[k] += list(f[k]['tofs'][()])
                vlscenters += list(f['vls']['centroids'][()])
                vlssums += list(f['vls']['sum'][()])
    _=[print('%s\t%i\t%i'%(k,len(tofs[k]),addresses[k][-1]+nedges[k][-1])) for k in portkeys]
    for k in portkeys:
        quants[k].setbins(data=tofs[k])

    vlsbins = [np.uint32(max(0,min(int(v-vlsoffset)//2,nvlsbins-1))) for v in vlscenters]
    vlsnorm = np.zeros(nvlsbins)
    for shot,vlsbin in enumerate(vlsbins):
        vlsnorm += vlssums[shot]
        for k in portkeys:
            a = addresses[k][shot]
            n = nedges[k][shot]
            try:
                #hist[k][vlsbin,:] += quants[k].histogram(tofs[k][a:a+n]).astype(float)
                hist[k][vlsbin,:] += quants[k].histogram(tofs[k][a:a+n]).astype(float)/quants[k].binwidths()
            except:
                print(vlsbin,a,n)

    crop=50
    fig1,ax = plt.subplots(2,4,figsize=(18,9))
    for i,k in enumerate(portkeys):
        X,Y = np.meshgrid(quants[k].binedges()[1:-30],np.arange(nvlsbins+1))
        if donorm:
            for v in range(nvlsbins):
                if vlsnorm[v]>0:
                    hist[k][v,:] /= vlsnorm[v]
                else:
                    hist[k][v,:] *= 0
        ax[i//4,i%4].pcolor(X,Y,np.log2(hist[k][:,1:-30]))#,origin='lower')
        ax[i//4,i%4].set_title('%s'%k)
        ax[i//4,i%4].set_xlabel('tofs')
        ax[i//4,i%4].set_ylabel('vls')
    plt.savefig('Figure_1_tofs_oneto%i.png'%crop)
    plt.show()

    fig2,ax = plt.subplots(2,4,figsize=(18,9))
    for i,k in enumerate(portkeys):
        X,Y = np.meshgrid(np.arange(len(quants[k].binedges())-31),np.arange(nvlsbins+1))
        if donorm:
            for v in range(nvlsbins):
                if vlsnorm[v]>0:
                    hist[k][v,:] /= vlsnorm[v]
                else:
                    hist[k][v,:] *= 0
        ax[i//4,i%4].pcolor(X,Y,np.log2(hist[k][:,1:-30]))#,origin='lower')
        ax[i//4,i%4].set_title('%s'%k)
        ax[i//4,i%4].set_xlabel('qbins')
        ax[i//4,i%4].set_ylabel('vls')
    plt.savefig('Figure_1_qbins_oneto%i.png'%crop)
    plt.show()
        #outname = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_vernier/ascii/test_%s_hist.dat'%(k)


    for k in ['port_0','port_1','port_14','port_5','port_12']:
        fg,ax = plt.subplots(1,1,figsize=(5,4))
        X,Y = np.meshgrid(quants[k].binedges()[1:-crop],np.arange(nvlsbins+1))
        ax.pcolor(X,Y,hist[k][:,1:-crop])
        ax.set_title('%s'%k)
        ax.set_ylabel('vls [arb]')
        ax.set_xlabel('ToF [arb]')
        plt.savefig('%s_tofs_oneto%i.png'%(k,crop))
        #plt.savefig('%s_tofs_oneto%i_raw.png'%(k,crop))
        plt.show()

    for k in ['port_0','port_1','port_14','port_5','port_12']:
        fg,ax = plt.subplots(1,1,figsize=(5,4))
        X,Y = np.meshgrid(np.arange(len(quants[k].binedges())-1-crop),np.arange(nvlsbins+1))
        ax.pcolor(X,Y,hist[k][:,1:-crop])
        ax.set_title('%s'%k)
        ax.set_ylabel('vls [arb]')
        ax.set_xlabel('Qbins [index]')
        plt.savefig('%s_qbins_oneto%i.png'%(k,crop))
        #plt.savefig('%s_qbins_oneto%i_raw.png'%(k,crop))
        plt.show()

    return

if __name__=='__main__':
    main()

