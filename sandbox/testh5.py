#!/usr/bin/python3

import numpy as np
import h5py
import sys

def hello():
    print('hello')
    return

def main():
    fname = 'test.h5'
    with h5py.File(fname,'w') as f:
        g = f.create_group('this_list')
        l =  [1,2,3,4,5,6]
        g.create_dataset('data',data=l)
        g.attrs['len']=len(l)

    print('wrote %s'%fname)
    return

if __name__ == '__main__':
    main()
