#!/usr/bin/python3

import numpy as np
import sys
import re

def main():
    filehead = './'
    filetail = 'hist.dat'
    s = np.zeros((150,),dtype=int)
    for fname in sys.argv[1:]:
        m = re.search('^(.+_)(\d+)(.hist.dat)',sys.argv[1])
        if m:
            filehead = m.group(1)
            port = m.group(2)
            filetail = m.group(3)
            data = np.loadtxt(fname)
            s += np.sum(data[2900:3050,555:562].astype(int),axis=1)
        else:
            print('syntax = main <list of files to build spectrum>')
            continue
    np.savetxt('%s.spec'%(filehead),s)
    return

if __name__ == '__main__':
    if len(sys.argv)>1:
        main()
    else:
        print('syntax = main <list of files to build spectrum>')
