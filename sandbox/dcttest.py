#!/usr/bin/python3

import numpy as np
import h5py
import math
from scipy.fftpack import dct,dst,idct,idst
import sys

def gauss(x,xc,w):
    y = np.exp(-np.power((x-xc)/w,int(2)))
    return y

def cossqfilt(x):
    f = np.power(np.cos(np.arange(len(x))*math.pi/16),int(2))
    res = x*f
    return res
def sinsqfilt(x):
    f = np.power(np.sin(np.arange(len(x))*math.pi/16),int(2))
    res = x*f
    return res

def thereandback(sig):
    s = sinsqfilt(sig)
    c = cossqfilt(sig)
    c = np.roll(c,-8)
    #smean = np.mean(s)
    #cmean = np.mean(c)
    #s -= smean
    #c -= cmean
    S = np.concatenate((dct(np.concatenate((s,np.flip(s,axis=0))),type=1), np.zeros(len(s)*6)))
    C = np.concatenate((dct(np.concatenate((c,np.flip(c,axis=0))),type=1), np.zeros(len(c)*6)))
    sback = dct(S,type=1)
    cback = dct(C,type=1)
    cback = np.roll(cback,32)
    #SS = dct(np.concatenate((s,-1*np.flip(s,axis=0))),type=1)
    #sc = dst(SC,type=1)/len(SC)*np.power(2/math.pi,3./2)
    #ss = dst(SS,type=1)/len(SC)*np.power(2/math.pi,3./2)
    #sback += smean
    #cback += cmean
    return sback[:len(S)//4],cback[:len(S)//4],sinsqfilt(np.ones(len(S)//4)),cossqfilt(np.ones(len(S)//4))

def deriv(s):
    s0 = s[0]
    s[-10:] *= np.flip(np.arange(10,dtype=float)/10.,axis=0)
    f = np.arange(len(s)*2,dtype=float)/len(s)/2.
    SS = dst(np.concatenate((s,-1*np.flip(s,axis=0))) )
    SC = dct(np.concatenate((s,np.flip(s,axis=0))) )
    ds = idct(-math.sqrt(2./math.pi)*s0 + 2.*math.pi*f*SS) 
    ds += idst(2.*math.pi*f*SC)
    '''
    sc_back = idkt( f*SC )
    ss_back = idct( f*SS )
    ds = sc_back - ss_back
    return math.sqrt(len(s))*ds[:len(s)]
    '''
    return math.sqrt(len(s))*ds

def main():
    scale = 1.
    if len(sys.argv)>1:
        scale = float(sys.argv[1])
    rng = np.random.default_rng(int(4)) # fixed seed for repeatability
    totlen = int(2**8)
    seclen = int(2**3)
    t = np.arange(totlen,dtype=float)/float(seclen)
    s = np.zeros(totlen,dtype=float)
    nrandpeaks = rng.poisson(scale*totlen/seclen)
    centers = rng.random(nrandpeaks)*float(totlen/seclen)
    widths = rng.poisson(2,nrandpeaks)/float(seclen)
    for i in range(nrandpeaks):
        s += gauss(t,float(centers[i]),float(widths[i]))
    ds = deriv(s)
    sback = np.cumsum(ds,axis=0)/math.sqrt(totlen)
    sc,ss,sq,cq = thereandback(s)
    np.savetxt('dcttest.dat',np.column_stack((np.concatenate((t,t+t[1]+t[-1])),np.concatenate((s,np.flip(s,axis=0))),ds,sback,sc,ss,sq,cq)),fmt='%.3f')
    return

if __name__ == '__main__':
    main()
