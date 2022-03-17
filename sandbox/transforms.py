#!/usr/bin/python3

from scipy.fftpack import dct,dst
import numpy as np

def show(x):
    offset = 5+np.abs(np.min(x))
    print(' '*(offset-1)+'|')
    _ = [print(' '*int(offset+v-1)+'*') for v in x]
    print(' '*(offset-1)+'|')
    return

def dctmat(sz):
    return 2*dct( np.identity(2*sz),axis=1 , type=2)[::2].T
def dstmat(sz):
    return 2*dst( np.identity(2*sz),axis=1 , type=2)[1::2].T
def idctmat(sz):
    return 1./2*dct( np.identity(2*sz),axis=0 , type=3)[::2].T
def idstmat(sz):
    return 1./2*dst( np.identity(2*sz),axis=0 , type=3)[1::2].T

def mydct(cmat,x):
    xin = np.append(x[::2],x[-1::-2])
    return np.inner(cmat,xin)[::2]

def mydst(smat,x):
    xin = np.append(x[1::2],-1*x[-2::-2])
    return np.inner(smat,xin)[1::2]

def myidct(Cmat,x):
    res = np.inner(Cmat,x)
    sz=len(x)
    xout = np.zeros(sz)
    xout[::2]=res[:sz//2]
    xout[1::2]=res[-1:-sz//2-1:-1]
    return xout.astype(int)//32

def myidst(Smat,x):
    res = np.inner(Smat,x)
    sz=len(x)
    xout = np.zeros(sz)
    xout[1::2]=res[:sz//2]
    xout[::2]=-1*res[-1:-sz//2-1:-1]
    return xout

def main():
    x = np.array([ 0, -4,  8,  16,  12,  4,  0, -4])
    x_sym = np.append(x,np.flip(x,axis=0))
    x_anti = np.append(x,np.flip(-x,axis=0))
    sz = len(x)
    show(x)

    c = dctmat(sz)
    s = dstmat(sz)
    C = idctmat(sz)
    S = idstmat(sz)
    print('mydct:')
    print(mydct(c,x))
    print(dct(x_sym)[::2])
    print('mydst:')
    print(mydst(s,x))
    print(dst(x_anti)[1::2])
    print('myidct(C,mydct(c,x)')
    show(myidct(C,mydct(c,x)).astype(int)//16)
    print(myidct(C,mydct(c,x)).astype(int))
    print('\n'*10)
    print( (100*(2*sz)*np.linalg.pinv((dct(np.identity(2*sz),type=2)))[:,0::2]).astype(int) )
    print( (100*0.50*dct(np.identity(2*sz),type=3)[:,::2]).astype(int) )
    print( (np.inner(C.T,np.inner(c,x)).astype(int)//32 ) )
    print(x)
    print( (np.inner(S.T,np.inner(s,x)).astype(int)//32 ) )
    print( (myidct(C,mydct(c,x)).astype(int) ) )

    return

if __name__ == '__main__':
    main()
else:
    x = np.array([ 0, -4,  8,  16,  12,  4,  0, -4])
    X = (dct(np.append(x,np.flip(x,axis=0)),type=2)).astype(int)
    xback=(dct(np.append(X,np.zeros(16),axis=0),type=3)/8).astype(int)
    cmat = dct( np.concatenate((np.identity(8),np.flip(np.identity(8),axis=1))) , type=2, axis = 0)
    d = np.zeros(x.shape[0]*2)
    d[::2]=1
    Cmat = dct( np.diag(d), axis=0,type=3 )
    Y = (dst(np.append(x,np.flip(-x,axis=0)),type=2)).astype(int)
    smat = dst( np.concatenate((np.identity(8),np.flip(-1*np.identity(8),axis=1))) , type=2, axis = 0)
