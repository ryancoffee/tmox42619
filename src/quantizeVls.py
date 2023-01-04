#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3
import numpy as np
import sys
import h5py
import re
import matplotlib.pyplot as plt
from Quantizers import Quantizer

def main():
    if len(sys.argv)<3:
        print('syntax: quantizeVls.py <ntofbins> <nvlsbins> <fnames>')
        return

    donorm = False
    fnames = sys.argv[3:]
    ntofbins = np.uint32(sys.argv[1])
    nvlsbins = np.uint32(sys.argv[2])
    vlsoffset = 256 
    tofs = {} 
    addresses = {} 
    nedges = {} 
    quants = {}
    hist = {}
    vlscenters = []
    vlssums = []
    portkeys = []
    for fname in fnames:
        vlspitchcorrect = 0
        m = re.search('run_(\d+)',fname)
        if m:
            vlspitchcorrect = 142 if int(m.group(1))<87 and int(m.group(1))>81 else 0
        with h5py.File(fname,'r') as f:
            portkeys = [k for k in f.keys() if (re.search('port',k) and not re.search('_16',k) and not re.search('_2',k))]
            if len(quants.keys())==0:
                for k in portkeys:
                    quants[k] = Quantizer(style='nonuniform',nbins=ntofbins)
                    tofs[k] = list(f[k]['tofs'][()])
                    addresses[k] = list(f[k]['addresses'][()].astype(np.uint64))
                    nedges[k] = list(f[k]['nedges'][()])
                    hist[k] = np.zeros((nvlsbins,ntofbins),dtype=float)
                vlscenters = list(f['vls']['centroids'][()]+vlspitchcorrect)
                vlssums = list(f['vls']['sum'][()])
            else:
                for k in portkeys:
                    offsetTofs = np.uint64(len(tofs[k]))
                    addresses[k] += [offsetTofs+np.uint64(a) for a in f[k]['addresses'][()]]
                    nedges[k] += list(f[k]['nedges'][()])
                    tofs[k] += list(f[k]['tofs'][()])
                vlscenters += list(f['vls']['centroids'][()]+vlspitchcorrect)
                vlssums += list(f['vls']['sum'][()])

    if len(tofs[k])>1:
        _=[print('%s\t%i\t%i'%(k,len(tofs[k]),addresses[k][-1]+nedges[k][-1])) for k in portkeys]

    for k in portkeys:
        if len(tofs[k])>1:
            quants[k].setbins(data=tofs[k])

    h,b = np.histogram(vlscenters,100)
    plt.plot(b[:-1],h,'.')
    plt.show()


    vlsbins = [np.uint32(max(0,min(int(v-vlsoffset)>>1,nvlsbins-1))) for v in vlscenters]
    print('len(vlsbins) = %i'%len(vlsbins))
    vlsnorm = np.zeros(nvlsbins)
    for shot,vlsbin in enumerate(vlsbins):
        vlsnorm += vlssums[shot]
        for k in portkeys:
            try:
                a = addresses[k][shot]
                n = nedges[k][shot]
                #hist[k][vlsbin,:] += quants[k].histogram(tofs[k][a:a+n]).astype(float)
                hist[k][vlsbin,:] += quants[k].histogram(tofs[k][a:a+n]).astype(float)/quants[k].binwidths()
            except:
                print('%i vlsbin failed'%vlsbin)

    crop=2
    ycrop=1
    fig1,ax = plt.subplots(2,4,figsize=(18,9))
    for i,k in enumerate(portkeys):
        if len(tofs[k])>1:
            X,Y = np.meshgrid(quants[k].binedges()[1:-crop-1],np.arange(nvlsbins+1-ycrop))
            if donorm:
                for v in range(nvlsbins):
                    if vlsnorm[v]>0:
                        hist[k][v,:] /= vlsnorm[v]
                    else:
                        hist[k][v,:] *= 0
            ax[i//4,i%4].pcolor(X,Y,np.log2(hist[k][:-ycrop,1:-crop]))#,origin='lower')
            ax[i//4,i%4].set_title('%s'%k)
            ax[i//4,i%4].set_xlabel('tofs')
            ax[i//4,i%4].set_ylabel('vls')
    plt.savefig('Figure_1_tofs_oneto%i.png'%crop)
    plt.show()

    fig2,ax = plt.subplots(2,4,figsize=(18,9))
    for i,k in enumerate(portkeys):
        if len(tofs[k])>1:
            X,Y = np.meshgrid(np.arange(len(quants[k].binedges())-crop-1),np.arange(nvlsbins+1-ycrop))
            if donorm:
                for v in range(nvlsbins):
                    if vlsnorm[v]>0:
                        hist[k][v,:] /= vlsnorm[v]
                    else:
                        hist[k][v,:] *= 0
            ax[i//4,i%4].pcolor(X,Y,np.log2(hist[k][:-ycrop,1:-crop]))#,origin='lower')
            ax[i//4,i%4].set_title('%s'%k)
            ax[i//4,i%4].set_xlabel('qbins')
            ax[i//4,i%4].set_ylabel('vls')
    plt.savefig('Figure_1_qbins_oneto%i.png'%crop)
    plt.show()
        #outname = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_vernier/ascii/test_%s_hist.dat'%(k)


    '''
    for k in ['port_0','port_1','port_14','port_5','port_12']:
        fg,ax = plt.subplots(1,1,figsize=(5,4))
        X,Y = np.meshgrid(quants[k].binedges()[1:-crop],np.arange(nvlsbins+1-ycrop))
        ax.pcolor(X,Y,hist[k][:-ycrop,1:-crop])
        ax.set_title('%s'%k)
        ax.set_ylabel('vls [arb]')
        ax.set_xlabel('ToF [arb]')
        plt.savefig('%s_tofs_oneto%i.png'%(k,crop))
        #plt.savefig('%s_tofs_oneto%i_raw.png'%(k,crop))
        plt.show()

    for k in ['port_0','port_1','port_14','port_5','port_12']:
        fg,ax = plt.subplots(1,1,figsize=(5,4))
        X,Y = np.meshgrid(np.arange(len(quants[k].binedges())-1-crop),np.arange(nvlsbins+1-ycrop))
        ax.pcolor(X,Y,hist[k][:-ycrop,1:-crop])
        ax.set_title('%s'%k)
        ax.set_ylabel('vls [arb]')
        ax.set_xlabel('Qbins [index]')
        plt.savefig('%s_qbins_oneto%i.png'%(k,crop))
        #plt.savefig('%s_qbins_oneto%i_raw.png'%(k,crop))
        plt.show()
        '''

    return

if __name__=='__main__':
    main()

