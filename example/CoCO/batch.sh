#!/bin/bash
#PBS -j oe
#PBS -W umask=022
#PBS -l ncpus=32

cd $PBS_O_WORKDIR
source ~/miniconda3/etc/profile.d/conda.sh
ulimit -s unlimited
conda activate lmp
exec python3 run.py > out.txt
