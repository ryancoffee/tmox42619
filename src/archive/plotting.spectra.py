#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3

import matplotlib.pyplot as plt
import h5py
import numpy as np
import sys
import re
import utils

anglesdict = {'port_4':90,'port_12':-90,'port_5':67.5,'port_13':-67.5,'port_14':-45,'port_15':-22.5,'port_1':22.5,'port_0':0}
def main(path,fname,runstr):
    makeplots = True
    correction = 6
    with h5py.File('%s/%s'%(path,fname),'r') as f:
        portkeys = [k for k in f.keys() if (re.search('port',k) and not re.search('_16',k) and not re.search('_2',k))] # keeping the bare MCP ports 2 and 16 here
        gmdscale = 0.5*(f['gmd']['bins'][()][:-1]+f['gmd']['bins'][()][1:])/f['gmd']['norm'][()]
        if makeplots:
            fig = plt.figure(figsize=(8,10))
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
            portorder = [4,12,5,13,14,15,1,0]
            portlist = ['port_%i'%i for i in portorder]
            for i,k in enumerate(portlist):
            #for k in portkeys:
                tofwidths = (f[k]['quantbins'][()][1:]-f[k]['quantbins'][()][:-1])/8./6.
                logwidths = (
                    np.log2((f[k]['quantbins'][()][1:]-f[k]['quantbins'][()][4])/8./6.) - 
                    np.log2((f[k]['quantbins'][()][:-1]-f[k]['quantbins'][()][4])/8./6.)
                )
                tofbins = 0.5*(f[k]['quantbins'][()][1:]+f[k]['quantbins'][()][:-1])/8./6.
                spect = np.sum(f[k]['hist'][()][4:-4,:],axis=0)
                for j,tw in enumerate(tofwidths):
                    spect[j] /= tw
                ax1.step(tofbins-tofbins[5]+correction,4.*i+spect*1e-3,label=k)
                ax2.step(np.log2(tofbins-tofbins[5]+correction),4.*i+spect*1e-3,label='%.1fdeg'%(anglesdict[k]))
            ax1.legend(bbox_to_anchor=(1.02,1),loc = 'upper left')
            ax2.legend(bbox_to_anchor=(1.02,1),loc = 'upper left')
            ax1.set_xlim((-2,46))
            ax1.set_xlabel('ToF [ns]')
            ax1.set_ylabel('signal arb. offset [k cnts/ns]')
            ax2.set_xlabel('log_2(ToF)')
            ax2.set_ylabel('signal arb. offset [k cnts/ns]')
            ax2.set_xlim((2.25,5.75))
            fig.tight_layout()
            plt.savefig('%s/../figs/spect_%s.png'%(path,runstr))
            plt.show()

    return

if __name__ == '__main__':
    if len(sys.argv)>1:
        m = re.search('(.*/h5files)/(hits.*(run_\d+)\.h5\.counts\..*\.h5)',sys.argv[1])
        if m:
            main(path=m.group(1),fname=m.group(2),runstr=m.group(3))
        else:
            print('failed the path/fname match')
    else:
        print('syntax: figs/plotting.spectra.py fnameQuantizedHist')
