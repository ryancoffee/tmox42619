#!/cds/sw/ds/ana/conda2/inst/envs/ps-4.5.7-py39/bin/python3

import sys
import numpy as np
import h5py
import re
import matplotlib.pyplot as plt

def main(fname):
    data = {}
    qbins = {}
    gmdwin = []
    portkeys = []
    x = np.arange(1<<8,dtype=float)
    valenceO={'port_0':'-','port_1':1930,'port_12':1477,'port_13':1308,'port_14':2534,'port_15':2323,'port_4':1502,'port_5':1317}
    valenceN={'port_0':'-','port_1':1855,'port_12':1442,'port_13':1280,'port_14':2365,'port_15':2203,'port_4':1471,'port_5':1292}
    augersN={'port_0':3712,'port_1':5127,'port_12':2000,'port_13':1627,'port_14':'-','port_15':'-','port_4':2057,'port_5':1663}
    augersO={'port_0':2099,'port_1':2241,'port_12':1594,'port_13':1383,'port_14':3406,'port_15':2915,'port_4':1618,'port_5':1406}
    photosO={'port_0':'-','port_1':'-','port_12':'-','port_13':4314,'port_14':'-','port_15':'-','port_4':'-','port_5':5477}
    photosN1={'port_0':'-','port_1':'-','port_12':4646,'port_13':2321,'port_14':'-','port_15':'-','port_4':'-','port_5':2444}
    photosN2={'port_0':'-','port_1':'-','port_12':4442,'port_13':2300,'port_14':'-','port_15':'-','port_4':'-','port_5':2416}
    with h5py.File(fname,'r') as f:
        portkeys = [k for k in f.keys() if re.search('port_',k)]
        print(portkeys)
        for k in portkeys:
            data.update({k:f[k]['hist'][()]})
            qbins.update({k:f[k]['qbins'][()]})
        x = np.arange(len(data[portkeys[0]]),dtype=float)
        gmdwin = list(f['gmdwin'][()])
    for i,k in enumerate(portkeys):
        qbincenters = 0.5*(qbins[k][:-1] + qbins[k][1:])
        qbinwidths = (qbins[k][1:] - qbins[k][:-1])
        xvals = qbins[k][:-1]-qbins[k][0]
        #plt.plot(i*.50+1./qbinwidths,label='%s, %i-%i uJ'%(k,gmdwin[0],gmdwin[1]))
        #plt.plot(xvals,i*50+np.sum(data[k],axis=1)/qbinwidths,label='%s, %i-%i uJ'%(k,gmdwin[0],gmdwin[1]))
        plt.plot(xvals,i*50+np.sum(data[k],axis=1)/qbinwidths,label='%s'%(k))
    plt.xlabel('Time-of-Flight [arb. units]')
    plt.ylabel('signal [counts]')
    plt.xlim((0,6000))
    plt.legend()
    plt.savefig('./figures/multiretardationNNO.png')
    plt.show()

    calib = {}
    xvals = {}
    yvals = {}
    for k in portkeys:
        calib.update({k:[]})
        if valenceO[k] != '-':
            calib[k] += [600.-41.6,valenceO[k]]
        if valenceN[k] != '-':
            calib[k] += [600.-37.3,valenceN[k]]
        if photosO[k] != '-':
            calib[k] += [600.-543.1,photosO[k]]
        if photosN1[k] != '-':
            calib[k] += [600.-409.9,photosN1[k]]
        if photosN2[k] != '-':
            calib[k] += [600.-405.9,photosN2[k]]
        if augersO[k] != '-':
            calib[k] += [501.,augersO[k]]
        if augersN[k] != '-':
            calib[k] += [373.,augersN[k]]
        print(calib[k])
        xvals.update({k:np.log2(calib[k][1::2])})
        yvals.update({k:np.log2(calib[k][0::2])})
        print(xvals[k],yvals[k])
        plt.plot(xvals[k],yvals[k],'-',label='%s'%(k))
    plt.legend()
    plt.xlabel('log2(ToF)')
    plt.ylabel('log2(photon Energy) [log2(eV)])')
    plt.savefig('./figures/multiretardationNNO_calibpoints.png')
    plt.show()




    print(gmdwin)
    return

if __name__ == '__main__':
    if len(sys.argv)<2:
        print('syntax: plotting.multiretardation.py <quantized fname>')
    else:
        main(sys.argv[1])

