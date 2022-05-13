#!/usr/bin/python3

import numpy as np
import h5py
import sys

def main():
    if len(sys.argv)<2:
        print('Give me an output calibfname!')
        return

    calibfname = sys.argv[1]
    with h5py.File(calibfname,'w-') as lt2ecalib:

        lt2ecalib.attrs['desc'] = 'usage is the logt output as scaled by multirun2Dhist_Neon.py, e.g. also only fitting log2t for inds 2000:4000, exp( log(x).dot(coeffs) ) since these fits were with log(e) vs log( log2(x) )... dont ask, Im nuts'
        lt2ecalib.create_dataset('x0',data=6.7)
        lt2ecalib.create_group('port_0').create_dataset('lt2e_coeffs',data=[4.90616,-4.09146,-0.46745,1.7479814]) # setting equal to channel 5 since fit cuts through the data well.
        lt2ecalib.create_group('port_1').create_dataset('lt2e_coeffs',data=[4.9273,-4.02607,-0.586989,1.45948])
        lt2ecalib.create_group('port_2').create_dataset('lt2e_coeffs',data=[4.9,-4.1,-0.5,1.7]) # default
        lt2ecalib.create_group('port_4').create_dataset('lt2e_coeffs',data=[4.91518,-4.06502,-0.504313,1.66375])
        lt2ecalib.create_group('port_5').create_dataset('lt2e_coeffs',data=[4.90989,-4.07983,-0.468298,1.744])
        lt2ecalib.create_group('port_12').create_dataset('lt2e_coeffs',data=[4.91625,-4.07671,-0.476034,1.75343])
        lt2ecalib.create_group('port_13').create_dataset('lt2e_coeffs',data=[4.91871,-4.07965,-0.542224,1.66011])
        lt2ecalib.create_group('port_14').create_dataset('lt2e_coeffs',data=[4.93082,-4.03442,-0.62086,1.44578])
        lt2ecalib.create_group('port_15').create_dataset('lt2e_coeffs',data=[4.91596,-4.08016,-0.474128,1.76402])
        lt2ecalib.create_group('port_16').create_dataset('lt2e_coeffs',data=[4.9,-4.1,-0.5,1.7]) # default

    return

if __name__ == '__main__':
    main()
