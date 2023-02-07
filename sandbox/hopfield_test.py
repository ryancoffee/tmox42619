#!/usr/bin/python3
import numpy as np


def main():
    '''
    XDATA = -1.+2.*np.eye(10)
    print(XDATA)
    inds = np.where(XDATA<0)
    YDATA = np.copy(XDATA)
    MASK = np.random.random(YDATA.shape)-.05
    inds = np.where(MASK<0)
    YDATA[inds] = 1.

    W = []
    for i in [0,1,2,4,5]:
        W += [np.outer(XDATA[:,i],XDATA[:,i]) - np.eye(XDATA.shape[0])]
    Wtot = np.sum(np.array(W),axis=0)
    for i in [0,1,2,4,5]:
        print(i,YDATA[:,i],np.inner(YDATA[:,i],Wtot))

    print('Now moving to binary/int arrays')
    '''
    XDATA = np.eye(10)
    YDATA = np.copy(XDATA)
    MASK = np.random.random(YDATA.shape)-.25
    inds = np.where(MASK<0)
    YDATA[inds] = 1.

    W = []
    for i in range(0,XDATA.shape[0]//2):
        W += [-np.outer(YDATA[:,i],YDATA[:,i]) + np.eye(XDATA.shape[0])]
    Wtot = np.sum(np.array(W),axis=0)/len(W)
    for i in range(XDATA.shape[0]//2,XDATA.shape[0]):
        print(i,YDATA[:,i].astype(int),(np.inner(YDATA[:,i],Wtot)*10).astype(int))
        #print(i,YDATA[:,i].astype(int),(np.inner(YDATA[:,i],Wtot)+.5).astype(int))

    return

if __name__ == '__main__':
    main()
