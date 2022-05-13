#!/usr/bin/python3

import h5py
import numpy as np
import sys
from utils import pkey

class Calib:
    def __init__(self):
        self.calib = {}
        self.rets = []

    def setrets(self,rets=[50,125]):
        self.rets = rets
        for r in self.rets:
            self.calib['vret_%i'%r] = {'runs':[]} 
            self.calib['vret_%i'%r].update( {'phens':[]} )
        return self

    def addphen(self,ret,phen):
        self.calib['vret_%i'%ret]['phens'] += [phen]
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
                retgrp.attrs['phens']=self.calib['vret_%i'%r]['phens']
                for p in self.ports:
                    #print(self.calib['vret_%i'%r]['port_%i'%p].keys())
                    prtgrp = retgrp.create_group('port_%i'%p)
                    prtgrp.attrs['t0'] = self.t0s['port_%i'%p]
                    prtgrp.create_dataset('energies',data=list(self.calib['vret_%i'%r]['port_%i'%p].keys()))
                    prtgrp.create_dataset('log2tofs',data=list(self.calib['vret_%i'%r]['port_%i'%p].values()))
        return self


def main():

    if len(sys.argv)<2:
        print('give output .h5 filename as arguement')
        return
    argoncalib = Calib().setrets(rets=[0]).setports(ports=[0,1,2,4,5,12,13,14,15,16])
    #argoncalib.sett0s( {'port_0':109830,'port_1':100451,'port_2':99810,'port_4':97180,'port_5':99071,'port_12':98561,'port_13':98657,'port_14':100331,'port_15':106956,'port_16':97330} )
    argoncalib.sett0s( {pkey(0):73227,pkey(1):66973,pkey(2):66545,pkey(4):64796,pkey(5):66054,pkey(12):65712,pkey(13):65777,pkey(14):66891,pkey(15):71312,pkey(16):64887} ) # final, based on inflate=4 nr_expand=4


    ###########################################
    ## this block from 08292021 hand fitting ##
    ###########################################
    ###########################################
    argoncalib.addrun(ret=0,run=7)
    argoncalib.addphen(ret=0,phen=600)
    argoncalib.addpoints(ret=0,port=0, pointsdict = {584.2:11.3167,              351.0:11.7049,              205.233:12.0905,203.485:12.0973} )
    argoncalib.addpoints(ret=0,port=1, pointsdict = {584.2:11.3185,570.7:11.3368,351.0:11.7066,273.7:11.8933,205.233:12.0929,203.485:12.0984} )
    argoncalib.addpoints(ret=0,port=4, pointsdict = {584.2:11.3166,570.7:11.3349,351.0:11.7034,273.7:11.8913,205.233:12.0908,203.485:12.0977} )
    argoncalib.addpoints(ret=0,port=5, pointsdict = {584.2:11.3191,570.7:11.3379,351.0:11.7051,273.7:11.8924,205.233:12.0922,203.485:12.0994} )
    argoncalib.addpoints(ret=0,port=12, pointsdict = {584.2:11.3180,570.7:11.3359,351.0:11.7036,273.7:11.8901,205.233:12.0922,203.485:12.0986} )
    argoncalib.addpoints(ret=0,port=13, pointsdict = {584.2:11.3164,570.7:11.3360,351.0:11.7032,273.7:11.8908,205.233:12.0905,203.485:12.0969} )
    argoncalib.addpoints(ret=0,port=14, pointsdict = {584.2:11.3191,570.7:11.3374,351.0:11.7063,273.7:11.8913,205.233:12.0930,203.485:12.0998} )
    argoncalib.addpoints(ret=0,port=15, pointsdict = {584.2:11.3161,570.7:11.3353,351.0:11.7034,273.7:11.8907,205.233:12.0905,203.485:12.0969} )

    argoncalib.addrun(ret=0,run=59)
    argoncalib.addphen(ret=0,phen=700)
    argoncalib.addpoints(ret=0,port=0, pointsdict = {684.2:11.2007,              451.0:11.5187,373.7:11.6688,206.233:12.0909,203.485:12.0973} )
    argoncalib.addpoints(ret=0,port=1, pointsdict = {684.2:11.2003,670.7:11.2165,451.0:11.5198,373.7:11.6622,205.233:12.0928,203.485:12.0987} )
    argoncalib.addpoints(ret=0,port=4, pointsdict = {684.2:11.1996,670.7:11.2150,451.0:11.5179,373.7:11.6625,205.233:12.0911,203.485:12.0979} )
    argoncalib.addpoints(ret=0,port=5, pointsdict = {684.2:11.2013,670.7:11.2163,451.0:11.5198,373.7:11.6623,205.233:12.0923,203.485:12.0992} )
    argoncalib.addpoints(ret=0,port=12, pointsdict = {684.2:11.2007,670.7:11.2152,451.0:11.5196,373.7:11.6619,205.233:12.0926,203.485:12.0986} )
    argoncalib.addpoints(ret=0,port=13, pointsdict = {684.2:11.1990,670.7:11.2143,451.0:11.5181,373.7:11.6609,205.233:12.0906,203.485:12.0966} )
    argoncalib.addpoints(ret=0,port=14, pointsdict = {684.2:11.2022,670.7:11.2171,451.0:11.5209,373.7:11.6646,205.233:12.0928,203.485:12.0996} )
    argoncalib.addpoints(ret=0,port=15, pointsdict = {684.2:11.1996,670.7:11.2141,451.0:11.5185,373.7:11.6615,205.233:12.0902,203.485:12.0970} )


    ###########################################
    ###########################################

    argoncalib.h5out(sys.argv[1])

    return

    ###########################################
    ###########################################
    ## this block from 08282021 hand fitting ##
    '''
    argoncalib.addrun(ret=0,run=59)
    argoncalib.addrun(ret=50,run=61)
    argoncalib.addrun(ret=100,run=62)
    argoncalib.addrun(ret=150,run=63)
    argoncalib.addpoints(ret=0,port=0, pointsdict = {201.1:12.693,203.25:12.685,205.23:12.678,205.63:12.677,207.23:12.672} )
    argoncalib.addpoints(ret=0,port=1, pointsdict = {201.1:12.694,203.25:12.685,205.23:12.679,205.63:12.677,207.23:12.672} )
    argoncalib.addpoints(ret=0,port=4, pointsdict = {201.1:12.694,203.25:12.686,205.23:12.679,205.63:12.678,207.23:12.672} )
    argoncalib.addpoints(ret=0,port=5, pointsdict = {201.1:12.695,203.25:12.686,205.23:12.680,205.63:12.678,207.23:12.673} )
    argoncalib.addpoints(ret=0,port=12, pointsdict = {201.1:12.694,203.25:12.685,205.23:12.679,205.63:12.677,207.23:12.672} )
    argoncalib.addpoints(ret=0,port=13, pointsdict = {201.1:12.692,203.25:12.684,205.23:12.678,205.63:12.676,207.23:12.671} )
    argoncalib.addpoints(ret=0,port=14, pointsdict = {201.1:12.694,203.25:12.686,205.23:12.679,205.63:12.678,207.23:12.673} )
    argoncalib.addpoints(ret=0,port=15, pointsdict = {201.1:12.694,203.25:12.685,205.23:12.679,205.63:12.677,207.23:12.672} )

    argoncalib.addpoints(ret=50,port=0, pointsdict = {151.1:12.8500,153.25:12.8395,155.23:12.8313,155.63:12.8293,157.23:12.8225} )
    argoncalib.addpoints(ret=50,port=1, pointsdict = {151.1:12.8513,153.25:12.8410,155.23:12.8327,155.63:12.8309,157.23:12.8242} )
    argoncalib.addpoints(ret=50,port=4, pointsdict = {151.1:12.8514,153.25:12.8411,155.23:12.8327,155.63:12.8312,157.23:12.8244} )
    argoncalib.addpoints(ret=50,port=5, pointsdict = {151.1:12.8518,153.25:12.8415,155.23:12.8333,155.63:12.8313,157.23:12.8248} )
    argoncalib.addpoints(ret=50,port=12, pointsdict = {151.1:12.8521,153.25:12.8411,155.23:12.8332,155.63:12.8313,157.23:12.8248} )
    argoncalib.addpoints(ret=50,port=14, pointsdict = {151.1:12.8520,153.25:12.8413,155.23:12.8335,155.63:12.8310,157.23:12.8246} )
    argoncalib.addpoints(ret=50,port=15, pointsdict = {151.1:12.8513,153.25:12.8405,155.23:12.8327,155.63:12.8308,157.23:12.8239} )

    argoncalib.addpoints(ret=100,port=0, pointsdict = {101.1:13.0763,103.25:13.0610,105.23:13.0501,105.63:13.0475,107.23:13.0378} )
    argoncalib.addpoints(ret=100,port=1, pointsdict = {101.1:13.0784,103.25:13.0631,105.23:13.0524,105.63:13.0501,107.23:13.0398} )
    argoncalib.addpoints(ret=100,port=4, pointsdict = {101.1:13.0779,103.25:13.0630,105.23:13.0527,105.63:13.0497,107.23:13.0399} )
    argoncalib.addpoints(ret=100,port=5, pointsdict = {101.1:13.0792,103.25:13.0641,105.23:13.0531,105.63:13.0507,107.23:13.0407} )
    argoncalib.addpoints(ret=100,port=12, pointsdict = {101.1:13.0801,103.25:13.0650,105.23:13.0539,105.63:13.0517,107.23:13.0413} )
    argoncalib.addpoints(ret=100,port=14, pointsdict = {101.1:13.0791,103.25:13.0640,105.23:13.0530,105.63:13.0508,107.23:13.0406} )
    argoncalib.addpoints(ret=100,port=15, pointsdict = {101.1:13.0786,103.25:13.0635,105.23:13.0526,105.63:13.0506,107.23:13.0402} )

    argoncalib.addpoints(ret=150,port=0, pointsdict = {51.1:13.4688,53.25:13.4406,55.23:13.4208,55.63:13.4163,57.23:13.3982} )
    argoncalib.addpoints(ret=150,port=1, pointsdict = {51.1:13.4728,53.25:13.4442,55.23:13.4246,55.63:13.4201,57.23:13.4019} )
    argoncalib.addpoints(ret=150,port=4, pointsdict = {51.1:13.4735,53.25:13.4443,55.23:13.4252,55.63:13.4198,57.23:13.4020} )
    argoncalib.addpoints(ret=150,port=5, pointsdict = {51.1:13.4744,53.25:13.4460,55.23:13.4259,55.63:13.4215,57.23:13.4035} )
    argoncalib.addpoints(ret=150,port=12, pointsdict = {51.1:13.4766,53.25:13.4479,55.23:13.4280,55.63:13.4231,57.23:13.4052} )
    argoncalib.addpoints(ret=150,port=14, pointsdict = {51.1:13.4749,53.25:13.4463,55.23:13.4265,55.63:13.4216,57.23:13.4035} )
    argoncalib.addpoints(ret=150,port=15, pointsdict = {51.1:13.4745,53.25:13.4459,55.23:13.4262,55.63:13.4212,57.23:13.4032} )
    '''
    ###########################################
    ###########################################


    ###########################################
    ###########################################
    # Argon Auger energies = 
    '''
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
    '''


if __name__ == '__main__':
    main()
