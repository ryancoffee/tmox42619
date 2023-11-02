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
        labellist = []
        portkeys = ['port_13','port_12','port_0','port_4']
        i = 0
        for k in portkeys:
            plotk = False
            h = np.sum(q[k]['hist'][()],axis=1)
            qb = q[k]['qbins'][()]
            qdiff = qb[1:]-qb[:-1]
            print(np.min(qdiff))
            if re.search('port_4',k):
                plotk=True
                labellist += ['horiz (north)']
            if re.search('port_0',k):
                plotk=True
                labellist += ['vert']
            if re.search('port_12',k):
                plotk=True
                labellist += ['horiz (south)']
            if re.search('port_13',k):
                plotk=True
                labellist += ['\'13\'']
            if plotk:
                plt.stairs((i*100.)+h/qdiff,(qb-qb[0])/8/6)
                i += 1
        plt.xlabel('ToF [ns]')
        #plt.xlim(40,60)
        plt.xlim(20,120)
        plt.legend(labellist)
        plt.tight_layout()
        plt.savefig('./figures/qbinsRecovered_13_12_0_4.png')
        plt.show()
        fig,axs = plt.subplots(1,4)
        axs[3].pcolor(q['port_12']['hist'][()][:,::32].T,cmap='Greys',vmin=0,vmax=3)
        axs[1].pcolor(q['port_0']['hist'][()][:,::32].T,cmap='Greys',vmin=0,vmax=3)
        axs[2].pcolor(q['port_4']['hist'][()][:,::32].T,cmap='Greys',vmin=0,vmax=3)
        axs[0].pcolor(q['port_13']['hist'][()][:,::32].T,cmap='Greys',vmin=0,vmax=3)
        axs[0].set_xlabel('qbins')
        axs[1].set_xlabel('qbins')
        axs[2].set_xlabel('qbins')
        axs[3].set_xlabel('qbins')
        axs[0].set_ylabel('shot number')
        axs[3].set_title('horiz (south)')
        axs[1].set_title('vert')
        axs[2].set_title('horiz (north)')
        axs[0].set_title('\'13\'')
        plt.savefig('./figures/qbinsSnow_ports_13_12_0_4.png')
        plt.show()

        '''
        N2O:
        0V retardation on port 13.
        25V retardation on port 5.
        150V retardation on port 12.
        175V retardation on port 4.
        325V retardation on port 0.
        350V retardation on port 1.
        425V retardation on port 15.
        450V retardation on port 14.
        '''
        portkeys = ['port_13','port_5','port_12','port_4','port_0','port_1','port_15','port_14']
        rets = [0,25,150,175,325,350,425,450]
        labellist = []
        scale=1.5
        for i,k in enumerate(portkeys):
            if i >2:
                scale=3
            if i >5:
                scale=10
            h = np.sum(q[k]['hist'][()],axis=1)
            qb = q[k]['qbins'][()]
            qdiff = qb[1:]-qb[:-1]
            plt.stairs((i*100.)+(scale)*h/qdiff,(qb-qb[0])/8/6)
            labellist += ['%s,ret=%i'%(k,rets[i])]
        plt.plot([29.5,29.5,30,30],[50,75,75,100],'-',color='k')
        plt.plot([30,30,34,34],[150,175,175,200],'-',color='k')
        plt.plot([34,34,35,35],[250,275,275,300],'-',color='k')
        plt.plot([35,35,45,45],[350,375,375,400],'-',color='k')
        plt.plot([45,45,49,49],[450,475,475,500],'-',color='k')
        plt.plot([49,49,64,64],[550,575,575,600],'-',color='k')
        plt.plot([64,64,75,75],[650,675,675,700],'-',color='k')
        plt.plot([35,36],[50,100],'-',color='k')
        plt.plot([36,43],[150,200],'-',color='k')
        plt.plot([43,45],[250,300],'-',color='k')
        plt.plot([45,80],[350,400],'-',color='k')
        plt.xlabel('ToF [ns]')
        plt.xlim(20,120)
        plt.tight_layout()
        plt.legend(labellist)
        plt.savefig('./figures/qbinsRecovered_retorder.png')
        plt.show()


    return


if __name__ == '__main__':
    if len(sys.argv)<2:
        print('./src/plotCOmpareQuantNonquant.py <quantizename>')
    else:
        main(sys.argv[1])

