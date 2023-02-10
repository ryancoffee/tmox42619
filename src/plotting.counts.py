#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3

import matplotlib.pyplot as plt
import h5py
import numpy as np
import sys
import re
import utils


def main(path,fname,runstring):
    print('%s/%s'%(path,fname))
    #path = '/reg/data/ana16/tmo/tmox42619/scratch/ryan_output_vernier_1000vlsthresh/h5files/'
    #fname = 'hits.tmox42619.run_188.h5.counts.qtofs.qgmd.h5'
    gmdinds = [i for i in range(1,32,5)]
    documsum = False
    makeplots = True
    dofitting = False
    
    run188lims = {}
    run188lims['O photo'] ={
        'port_12':(160,200,.25,10),
        'port_5':(155,190,.25,10),
        'port_14':(140,180,.25,10),
        'port_15':(125,200,.25,10),
        'port_0':(135,190,.25,10)
        }
    run188lims['N photo'] = {
        'port_12':(61,85,1,5),
        'port_5':(59,95,1,5),
        'port_14':(64,90,1,5),
        'port_15':(98,115,1,5),
        'port_0':(94,135,1,5)
        }
    run188lims['N Auger'] = {
        'port_12':(25,52,.5,5),
        'port_5':(20,55,.5,5),
        'port_14':(28,60,.5,5),
        'port_15':(35,90,.5,5),
        'port_0':(37,92,.5,5)
    }
    run188lims['O Auger'] = {
        'port_12':(2,24,1,5),
        'port_5':(2,24,1,5),
        'port_14':(2,30,1,5),
        'port_15':(3,35,1,5),
        'port_0':(3,35,1,5)
    }
    '''
    run199calibtofs={
        'port_12':((,),(,),(,),(,),(,))
        'port_5':((,),(,),(,),(,),(,))
        'port_14':((,),(,),(,),(,),(,))
        'port_15':((,),(,),(,),(,),(,))
        'port_0':((,),(,),(,),(,),(,))
        'port_1':((,),(,),(,),(,),(,))
        'port_13':((,),(,),(,),(,),(,))
        'port_4':((,),(,),(,),(,),(,))
    }
    run197calibtofs={# below threshold
        'port_12':((,),(,),(,),(,),(,))
        'port_5':((,),(,),(,),(,),(,))
        'port_14':((,),(,),(,),(,),(,))
        'port_15':((,),(,),(,),(,),(,))
        'port_0':((,),(,),(,),(,),(,))
        'port_1':((,),(,),(,),(,),(,))
        'port_13':((,),(,),(,),(,),(,))
        'port_4':((,),(,),(,),(,),(,))
    }
    '''
    run195calibtofs={# just above threshold (says 555eV, but really 10V off, I think it's really 545eV actual)
        'port_12':((223,3),(74,136.4),(81,132.4),(34,372),(7,503)),
        'port_5':((218,3),(70,136.4),(78,132.4),(33,372),(7,503)),
        'port_14':((202,3),(74,136.4),(80,132.4),(36,372),(7,503)),
        'port_15':((222,3),(106,136.4),(110,132.4),(46,372),(7,503)),
        'port_0':((202,3),(107,136.4),(111,132.4),(48,372),(7,503)),
        'port_1':((192,3),(86,136.4),(89,132.4),(39,372),(6,503)),
        'port_13':((199,3),(62,136.4),(70,132.4),(30,372),(6,503)),
        'port_4':((222,3),(70,136.4),(76,132.4),(30,372),(6,503))
    }
    '''
    run193calibtofs={# (560eV 550eV actual)
        'port_12':((,),(,),(,),(,),(,))
        'port_5':((,),(,),(,),(,),(,))
        'port_14':((,),(,),(,),(,),(,))
        'port_15':((,),(,),(,),(,),(,))
        'port_0':((,),(,),(,),(,),(,))
        'port_1':((,),(,),(,),(,),(,))
        'port_13':((,),(,),(,),(,),(,))
        'port_4':((,),(,),(,),(,),(,))
    }
    run190calibtofs={# (565eV, 555eV actual))
        'port_12':((,),(,),(,),(,),(,))
        'port_5':((,),(,),(,),(,),(,))
        'port_14':((,),(,),(,),(,),(,))
        'port_15':((,),(,),(,),(,),(,))
        'port_0':((,),(,),(,),(,),(,))
        'port_1':((,),(,),(,),(,),(,))
        'port_13':((,),(,),(,),(,),(,))
        'port_4':((,),(,),(,),(,),(,))
    }
    '''
    run188calibtofs={# (575, 565eV actual) Nitrogen energy is 421.1eV(hnu) - 8.5/12.5 (Central/Terminal) is binding energy From Hoshino et al., J. Phys. B: At. Mol. Opt. Phys. 51 (2018) 065402
                                                            # binding 412.6eV central : 408.6eV terminal, so 565-binding = (152.4/156.4 Central/Terminal)
                        # Angle‐resolved photoelectron spectroscopy of the core levels of N2O M. Schmidbauer, The Journal of Chemical Physics 94, 5299 (1991); doi: 10.1063/1.460514
                        # 372 eV is high energy head of Nc+Nt AUger, Paola Bolognesi et al., J. Chem. Phys. 125, 054306 (2006); https://doi.org/10.1063/1.2213254
                        # 503 eV is the high energy head of O Auger, F. P. Larkins: Auger spectra of nitrous oxide : J. Chem. Phys. 86, 3239 (1987); https://doi.org/10.1063/1.451982
        'port_12':[(180,23),(64,156.4),(71,152.4),(30,372),(7,503)],
        'port_5':[(180,23),(64,156.4),(71,152.4),(30,372),(7,503)],
        'port_14':[(180,23),(64,156.4),(71,152.4),(30,372),(7,503)],
        'port_15':[(180,23),(64,156.4),(71,152.4),(30,372),(7,503)],
        'port_0':[(180,23),(64,156.4),(71,152.4),(30,372),(7,503)],
        'port_1':[(150,23),(82,156.4),(86,152.4),(39,372),(7,503)],
        'port_13':[(155,23),(57,156.4),(63,152.4),(27,372),(7,503)],
        'port_4':[(177,23),(62,156.4),(67,152.4),(28,372),(7,503)]
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
                plt.show()
        if makeplots:
            for k in ['port_12','port_5','port_14','port_15','port_0']:
                XX,YY = np.meshgrid(f[k]['quantbins'][()][:-1]/8./6.,0.5*(f['gmd']['bins'][()][:-1]+f['gmd']['bins'][()][1:]))
                tofwidths = (f[k]['quantbins'][()][1:]-f[k]['quantbins'][()][:-1])/8./6.
                spect = np.ones(f[k]['hist'].shape,dtype=float)
                X = [(f[k]['quantbins'][()][v[0]]-f[k]['quantbins'][()][0])/8./6. for v in run188calibtofs[k][1:]]
                Y = [v[1] for v in run188calibtofs[k][1:]]
                XLOG2 = np.array(np.log2(X))
                YLOG2 = np.array(np.log2(Y))
                x0,theta = utils.fitpoly(XLOG2,YLOG2,order=2)
                for limlabel in run188lims.keys():
                    fig = plt.figure(figsize=(12,4))
                    ax1 = fig.add_subplot(111)
                    ax2 = ax1.twiny()
                    #ax3 = ax1.twiny()
                    theselims = run188lims[limlabel]
                    for i in range(spect.shape[0]):
                        spect[i,:] = f[k]['hist'][()][i,:]*gmdscale[i]/tofwidths
                    pc = ax1.pcolor(XX[1:,theselims[k][0]:theselims[k][1]],
                        YY[1:,theselims[k][0]:theselims[k][1]],
                        spect[1:,theselims[k][0]:theselims[k][1]])
                    etickinds = [i for i in range(theselims[k][0],theselims[k][1],int(theselims[k][3]))]
                    eticklocations = [f[k]['quantbins'][()][i]/8./6. for i in etickinds]
                    etickinds = [i for i in range(theselims[k][0],theselims[k][1],int(theselims[k][3]))]
                    enticklabels = []
                    eticklabels = ['%i'%i for i in etickinds]
                    for i in etickinds:
                        log2x = np.log2((f[k]['quantbins'][()][i]-f[k]['quantbins'][()][0])/8./6.)
                        log2y = np.sum([theta[j]*(log2x-x0)**int(j) for j in range(len(theta))])
                        enticklabels += ['%.1f'%2**log2y] 
                    #ax1.set_title('%s %s %s'%(runstring,limlabel,k))
                    fig.colorbar(pc,label='cnts./shot/ns',ax=ax1)
                    #fig.tight_layout()
                    pc.set_clim((0,theselims[k][2]))
                    ax1.set_xlabel('ToF [ns]')
                    ax2.set_xlim(ax1.get_xlim())
                    ax2.set_xticks(eticklocations)
                    ax2.set_xticklabels(eticklabels)
                    ax2.set_xlabel('quant bins')
                    #ax3.set_xlim(ax1.get_xlim())
                    #ax3.set_xticks(eticklocations)
                    #ax3.set_xticklabels(enticklabels)
                    #ax3.set_xlabel('Energy [eV]')
                    ax1.set_ylabel('Pulse Energy [uJ]')
                    plt.savefig('%s/../figs/%s_%s_%s.png'%(path,runstring,limlabel,k))
                    plt.show()
        if dofitting:
            for k in ['port_12','port_5','port_14','port_15','port_0','port_1','port_4','port_13']:
                #plt.loglog(X,Y,'.')
                #plt.plot(XLOG2,YLOG2,'.')
                xcurve = np.arange(2,8,step=.1,dtype=float)
                ycurve = [np.sum([theta[i]*(x-x0)**int(i) for i in range(len(theta))]) for x in xcurve]
                #plt.plot(xcurve,ycurve)
                #plt.xlim((3,8))
                #plt.ylim((4,10))
                plt.title(k)
                plt.loglog(X,Y,'.')
                plt.loglog([2**x for x in xcurve],[2.**y for y in ycurve])
                plt.show()
    return


if __name__ == '__main__':
    if len(sys.argv)>1:
        m = re.search('(.*/h5files)/(hits.*(run_\d+)\.h5\.counts\..*\.h5)',sys.argv[1])
        if m:
            main(m.group(1),m.group(2),m.group(3))
        else:
            print('failed the path/fname match')
    else:
        print('syntax: figs/plotting.counts.py fnameQuantizedHist')
