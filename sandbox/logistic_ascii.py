#!/usr/bin/python3

import numpy as np
import sys
import h5py
import re

def filldict(d,f):
    line = f.readline()
    nlines = 0 
    while line:
        nlines += 1
        if nlines>4096:
            break
        m = re.search('^(.+)\s+(.+)\s*$',line)
        if m:
            rstring = m.group(1)
            xstring = m.group(2)
            if rstring in d.keys():
                d[rstring] += [float(xstring)]
            else:
                d.update({rstring:[float(xstring)]})
        line = f.readline()
    return d

def main():
    nbins = 2**10
    if len(sys.argv)<2:
        print('Give a list of ascii files to process')
        return
    fnames = sys.argv[1:]
    print(fnames)
    accum = {}
    for fname in fnames:
        with open(fname,'r') as f:
            accum = filldict(accum,f)
    print(accum.keys())
    with h5py.File('logistic.h5','a') as f:
        for k in accum.keys():
            if k in f.keys():
                continue
            else:
                rgrp = f.create_group(k)
                h,b = np.histogram(accum[k],nbins)
                rgrp.create_dataset('accum',data=accum[k],dtype=np.float32)
                rgrp.create_dataset('hist',data=h,dtype=np.uint16)
                rgrp.create_dataset('bins',data=[b[0],b[-1]],dtype=np.float64)

    return

if __name__=='__main__':
    main()
