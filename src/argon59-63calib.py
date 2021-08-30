#!/usr/bin/python3

import numpy as np
import h5py
import sys
import re
from utils import mypoly,pkey

def compute_spectrum(c,ret,f,t0s):
    print(ret)
    print(f[pkey(0)].keys())
    data = {}
    ebins = np.arange(0,1024,.075)
    ltbins = np.arange(10,16,.0025)
    data.update({'ebins':ebins[:-1]})
    data.update({'ltbins':ltbins[:-1]})
    portlist = [0,1,4,5,12,14,15]
    for p in portlist:
        log2tofs = [np.log2(float(t)) for t in f[pkey(p)]['tofs'][()]-float(t0s[pkey(p)]) if t > 2**11 and t<2**15]
        theta = np.array(c[ret][pkey(p)].attrs['theta'],dtype=float)
        xvals = np.array(log2tofs - c[ret][pkey(p)].attrs['x0'],dtype=float)
        yvals = mypoly(xvals,len(theta)-1).dot(theta)
        data.update( {pkey(p):np.histogram(np.power(2.,yvals),ebins)[0]} )
        data.update( {'log2t_%s'%pkey(p):np.histogram(log2tofs,ltbins)[0]} )
    return data

def main():
    fname = 'test.h5'
    path = './'
    ret = {'run7':'vret_0','run59':'vret_0'}#,'run61':'vret_50','run62':'vret_100','run63':'vret_150'}
    ret.update( {'run8':'vret_0','run9':'vret_0','run10':'vret_0','run11':'vret_0','run61':'vret_0','run62':'vret_0','run63':'vret_0'} ) # artificially using same vret_0 calib for even the retarded cases
    runstring = 'run59'
    if len(sys.argv)<2:
        print('I need an h5 calibration file and a data file to proces (run7, 59) for now not yet 61-63)')
        print('%s <calibfile> <datafile>'%sys.argv[0])
        return

    m = re.search('^(.+/)calib/(.+)\.h5',sys.argv[1])
    if m:
        cpath = m.group(1)
        cname = m.group(2)
    else:
        print('failed to match calib filename')
        return

    m = re.search('^(.+/)h5files/(.+)\.(run.+)\.h5',sys.argv[2])
    if m:
        fpath = m.group(1)
        fname = m.group(2)
        runstring = m.group(3)
    else:
        print('failed to match data h5 filename')
        return

    cfname = '%scalib/%s.h5'%(cpath,cname)
    ifname = '%sh5files/%s.%s.h5'%(fpath,fname,runstring)

    t0s = {pkey(0):73227,pkey(1):66973,pkey(2):66545,pkey(4):64796,pkey(5):66054,pkey(12):65712,pkey(13):65777,pkey(14):66891,pkey(15):71312,pkey(16):64887} # hard coding this for runs 7-63
    print(t0s)


    with h5py.File(ifname,'r') as f:
        with h5py.File(cfname,'r') as c:
            data = compute_spectrum(c,ret[runstring],f,t0s)
            for p in (0,1,4,5,12,14,15):
                np.savetxt('%sascii/%s.%s.%s.spec'%(fpath,fname,pkey(p),runstring),np.column_stack((data['ebins'],data[pkey(p)])),fmt='%.2f')
                np.savetxt('%sascii/%s.%s.%s.log2t'%(fpath,fname,pkey(p),runstring),np.column_stack((data['ltbins'],data['log2t_%s'%pkey(p)])),fmt='%.2f')



    return

if __name__ == '__main__':
    main()

