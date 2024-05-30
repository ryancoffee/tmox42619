#!/usr/bin/bash
# assume 5 hours for each run of about 50k shots

#SBATCH --partition=roma
#
#SBATCH --job-name=hits2h5_min
#SBATCH --output=output-%j.txt
#SBATCH --error=output-%j.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8g
#
#SBATCH --time=0-08:00:00
#
#SBATCH --gpus 0

source /sdf/group/lcls/ds/ana/sw/conda2/manage/bin/psconda.sh
export scratchpath=/sdf/data/lcls/ds/tmo/tmox42619/scratch/ryan_output_slurm/h5files
if ! [ -f $scratchpath ]; then
	mkdir -p $scratchpath
fi
export datapath=/sdf/data/lcls/ds/tmo/tmox42619/xtc
export expname=tmox42619
export nshots=100000
export configfile=${scratchpath}/${expname}.hsdconfig.h5
printf -v runstr "r%04d" $1
if [ -f ${datapath}/${expname}-${runstr}-s000-c000.xtc2 ]; then
	python3 ./src/set_configs.py ${configfile}
	python3 ./src/hits2h5_minimal.py $1
else
	echo "XTC2 file not found for run ${expname}:${runstr}"
fi

