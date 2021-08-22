#!/usr/bin/python3

import h5py
import numpy as np
import sys

class Calib:
    def __init__(self):
        self.calib = {}
        self.rets = []

    def setrets(self,rets=[50,125]):
        self.rets = rets
        for r in self.rets:
            self.calib['vret_%i'%r] = {'runs':[]}
        return self

    def addrun(self,ret,run):
        self.calib['vret_%i'%ret]['runs'] += [run]
        return self

    def setports(self,ports=[0,1,2,4,5,12,13,14,15,16]):
        self.ports = ports
        self.t0s = {}
        for r in self.rets:
            for p in ports:
                self.calib['vret_%i'%r]['port_%i'%p] = {}
        return self

    def addpoints(self,ret,port,pointsdict):
        self.calib['vret_%i'%ret]['port_%i'%port].update(pointsdict)
        return self

    def sett0s(self,t0s):
        self.t0s.update( t0s )
        return self

    def h5out(self,fname):
        print('storing to {}'.format(fname))
        with h5py.File(fname,'w') as f:
            for r in self.rets:
                retgrp = f.create_group('vret_%i'%r) 
                retgrp.attrs['runs']=self.calib['vret_%i'%r]['runs']
                for p in self.ports:
                    #print(self.calib['vret_%i'%r]['port_%i'%p].keys())
                    prtgrp = retgrp.create_group('port_%i'%p)
                    prtgrp.attrs['t0'] = self.t0s['port_%i'%p]
                    prtgrp.create_dataset('energies',data=list(self.calib['vret_%i'%r]['port_%i'%p].keys()))
                    prtgrp.create_dataset('indices',data=list(self.calib['vret_%i'%r]['port_%i'%p].values()))
        return self


def main():

    if len(sys.argv)<2:
        print('give output .h5 filename as arguement')
    argoncalib = Calib().setrets(rets=[50,125]).setports(ports=[0,1,2,4,5,12,13,14,15,16])
    argoncalib.sett0s( {'port_0':109830,'port_1':100451,'port_2':99810,'port_4':97180,'port_5':99071,'port_12':98561,'port_13':98657,'port_14':100331,'port_15':106956,'port_16':97330} )



    argoncalib.addrun(ret=50,run=61)
    argoncalib.addpoints(ret=50,port=0, pointsdict = {634.2:3644,400.5:4621,155.5:7291} )
    argoncalib.addpoints(ret=50,port=1, pointsdict = {634.2:3647,400.5:4622,323.7:5153,155.5:7297} )
    argoncalib.addpoints(ret=50,port=2, pointsdict = {634.2:3647,400.5:4526,323.7:5005,155.5:6762} )
    argoncalib.addpoints(ret=50,port=4, pointsdict = {634.2:3647,400.5:4619,323.7:5153,155.5:7292} )
    argoncalib.addpoints(ret=50,port=5, pointsdict = {634.2:3647,400.5:4622,323.7:5157,155.5:7292} )
    argoncalib.addpoints(ret=50,port=12, pointsdict = {634.2:3647,400.5:4619,323.7:5153,155.5:7297} )
    argoncalib.addpoints(ret=50,port=13, pointsdict = {684.2:3533,450.5:4403,373.7:4863,205.5:6545} )
    argoncalib.addpoints(ret=50,port=14, pointsdict = {634.2:3647,400.5:4619,323.7:5149,155.5:7297} )
    argoncalib.addpoints(ret=50,port=15, pointsdict = {634.2:3647,400.5:4619,323.7:5153,155.5:7292} )
    argoncalib.addpoints(ret=50,port=16, pointsdict = {} )

    argoncalib.addrun(ret=50,run=39)
    argoncalib.addpoints(ret=50,port=0, pointsdict = {334.2:5015,320.7:5067,155.5:7291}) #,23.7:9022} )
    argoncalib.addpoints(ret=50,port=1, pointsdict = {334.2:5017,320.7:5080,155.5:7297})#,23.7:9022} )
    argoncalib.addpoints(ret=50,port=2, pointsdict = {370.7:4941,205.5:6754})#,73.7:8000} )
    argoncalib.addpoints(ret=50,port=4, pointsdict = {334.2:5017,320.7:5079,155.5:7299})#,23.7:9010} )
    argoncalib.addpoints(ret=50,port=5, pointsdict = {334.2:5015,320.7:5086,155.5:7302})#,23.7:9005} )
    argoncalib.addpoints(ret=50,port=12, pointsdict = {334.2:5017,320.7:5082,155.5:7302})#,23.7:9023} )
    argoncalib.addpoints(ret=50,port=13, pointsdict = {384.2:4749,370.7:4802,205.5:6551})#,73.7:7693} )
    argoncalib.addpoints(ret=50,port=14, pointsdict = {334.2:5015,320.7:5082,155.5:7299})#,23.7:9045} )
    argoncalib.addpoints(ret=50,port=15, pointsdict = {334.2:5017,320.7:5077,155.5:7299})#,23.7:9036} )
    argoncalib.addpoints(ret=50,port=16, pointsdict = {} )

    argoncalib.addrun(ret=125,run=44)
    argoncalib.addpoints(ret=125,port=0, pointsdict = {259.2:5534,245.7:5622,80.5:9442,25.5:15731} )
    argoncalib.addpoints(ret=125,port=1, pointsdict = {259.2:5537,245.7:5628,80.5:9467,25.5:15793} )
    argoncalib.addpoints(ret=125,port=2, pointsdict = {370.7:4953,205.5:7267,150.5:8000} )
    argoncalib.addpoints(ret=125,port=4, pointsdict = {259.2:5537,245.7:5625,80.5:9462,25.5:15908} )
    argoncalib.addpoints(ret=125,port=5, pointsdict = {259.2:5540,245.7:5622,80.5:9467,25.5:15917} )
    argoncalib.addpoints(ret=125,port=12, pointsdict = {259.2:5540,245.7:5619,80.5:9478,25.5:15835} )
    argoncalib.addpoints(ret=125,port=13, pointsdict = {384.2:4744,370.7:4792,205.5:6543,150.5:7679} )
    argoncalib.addpoints(ret=125,port=14, pointsdict = {259.2:5537,245.7:5622,80.5:9462,25.5:15873} )
    argoncalib.addpoints(ret=125,port=15, pointsdict = {259.2:5537,245.7:5622,80.5:9467,25.5:16060} )
    argoncalib.addpoints(ret=125,port=16, pointsdict = {} )


    argoncalib.h5out(sys.argv[1])

    return

if __name__ == '__main__':
    main()
