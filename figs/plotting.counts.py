#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3

import matplotlib.pyplot as plt
import h5py
import numpy as np
import sys
import re


def main(path,fname,runstring):
    print('%s/%s'%(path,fname))
    #path = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_vernier_1000vlsthresh/h5files/'
    #fname = 'hits.tmox42619.run_188.h5.counts.qtofs.qgmd.h5'
    gmdinds = [i for i in range(1,32,5)]
    documsum = False
    
    run188lims = {}
    run188lims['O photo'] ={
        'port_12':(80,100,.25),
        'port_14':(70,90,.25),
        'port_0':(75,95,.25)
        }
    run188lims['N photo'] = {
        'port_12':(30,50,1),
        'port_14':(32,50,1),
        'port_0':(50,75,1)
        }
    run188lims['N Auger'] = {
        'port_12':(12,30,1),
        'port_14':(15,32,1),
        'port_0':(20,50,1)
    }
    run188lims['O Auger'] = {
        'port_12':(2,12,1),
        'port_14':(2,15,1),
        'port_0':(2,20,1)
    }
    with h5py.File('%s/%s'%(path,fname),'r') as f:
        portkeys = [k for k in f.keys() if (re.search('port',k) and not re.search('_16',k) and not re.search('_2',k))] # keeping the bare MCP ports 2 and 16 here
        gmdscale = 0.5*(f['gmd']['bins'][()][:-1]+f['gmd']['bins'][()][1:])/f['gmd']['norm'][()]
        if documsum:
            for k in portkeys:
                _= [plt.plot(f[k]['quantbins'][()][:-1]/8./6.,np.cumsum(f[k]['hist'][()][i,:]*gmdscale[i])) for i in gmdinds]
                plt.xlabel('tof [ns]')
                plt.ylabel('mean cumulative counts/shot')
                plt.ylim((0,25))
                plt.xlim((690,1300))
                plt.title('%s %s'%(runstring,k))
                plt.grid()
                plt.legend(['%i uJ'%int(0.5*(f['gmd']['bins'][()][i]+f['gmd']['bins'][()][i+1])+0.5) for i in gmdinds],bbox_to_anchor=(1.05,1),loc='upper left')
                plt.tight_layout()
                plt.savefig('/cds/home/c/coffee/Downloads/%s_cumcounts_%s.png'%(runstring,k))
                #plt.show()
        for k in ['port_12','port_14','port_0']:
            XX,YY = np.meshgrid(f[k]['quantbins'][()][:-1]/8./6.,0.5*(f['gmd']['bins'][()][:-1]+f['gmd']['bins'][()][1:]))
            tofwidths = (f[k]['quantbins'][()][1:]-f[k]['quantbins'][()][:-1])/8./6.
            spect = np.ones(f[k]['hist'].shape,dtype=float)
            for limlabel in run188lims.keys():
                theselims = run188lims[limlabel]
                for i in range(spect.shape[0]):
                    spect[i,:] = f[k]['hist'][()][i,:]*gmdscale[i]/tofwidths
                plt.pcolor(XX[1:,theselims[k][0]:theselims[k][1]],
                    YY[1:,theselims[k][0]:theselims[k][1]],
                    spect[1:,theselims[k][0]:theselims[k][1]])
                plt.title('%s %s %s'%(runstring,limlabel,k))
                plt.colorbar(label='cnts./shot/ns')
                plt.clim((0,theselims[k][2]))
                plt.xlabel('ToF [ns]')
                plt.ylabel('Pulse Energy [uJ]')
                plt.show()


if __name__ == '__main__':
    if len(sys.argv)>1:
        m = re.search('(.*/h5files)/(hits.*(run_\d+)\.h5\.counts\..*\.h5)',sys.argv[1])
        if m:
            main(m.group(1),m.group(2),m.group(3))
        else:
            print('failed the path/fname match')
    else:
        print('syntax: figs/plotting.counts.py fnameQuantizedHist')
