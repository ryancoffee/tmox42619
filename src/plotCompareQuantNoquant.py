#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3

import numpy as np
import h5py
import matplotlib.pyplot as plt
import sys
import re

verbose = True

def main(fname):
    with h5py.File(fname,'r') as q:
        if verbose: 
            print(q.keys())
        portkeys = [k for k in q.keys() if re.search('port',k)]
        if verbose: 
            _=[print(k) for k in portkeys]
        for i,k in enumerate(portkeys):
            h = np.sum(q[k]['hist'][()],axis=1)
            qb = q[k]['qbins'][()]
            qdiff = qb[1:]-qb[:-1]
            print(np.min(qdiff))
            plt.stairs((i*20.)+h/qdiff,(qb-qb[0])/8/6)
        plt.xlabel('ToF [ns]')
        #plt.xlim(40,60)
        plt.xlim(0,160)
        plt.legend([k for k in portkeys])
        plt.tight_layout()
        plt.savefig('./figures/qbinsRecovered.png')
        plt.show()
        fig,axs = plt.subplots(1,2)
        axs[0].pcolor(q['port_12']['hist'][()][:,::32].T,cmap='Greys',vmin=0,vmax=3)
        axs[1].pcolor(q['port_0']['hist'][()][:,::32].T,cmap='Greys',vmin=0,vmax=3)
        axs[0].set_xlabel('qbins')
        axs[1].set_xlabel('qbins')
        axs[0].set_ylabel('shot number')
        axs[0].set_title('horiz')
        axs[1].set_title('vert')
        plt.savefig('./figures/qbinsSnow_ports_12_0.png')
        plt.show()


    return


if __name__ == '__main__':
    if len(sys.argv)<2:
        print('./src/plotCOmpareQuantNonquant.py <quantizename>')
    else:
        main(sys.argv[1])

