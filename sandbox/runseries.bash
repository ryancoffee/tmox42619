#!/usr/bin/bash
rstart=2.75
rstop=3.75
rstep=0.01
for r in $(seq $rstart $rstep $rstop); do 
	./logistichist.py 16384 256 $r /media/coffee/9C33-6BBD/logistic/hist_2.25_3.75.h5; 
	echo "finished $r" 
done
