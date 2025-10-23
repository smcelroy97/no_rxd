#!/bin/bash
#$ -cwd
#$ -N SMgridSearch
#$ -q cpu.q
#$ -pe smp 64
#$ -l h_vmem=475G
#$ -l h_rt=2:40:00
#$ -o /ddn/smcelroy97/no_rxd/singleSim.out
#$ -e /ddn/smcelroy97/no_rxd/singleSim.err

source ~/.bashrc
conda activate "netpyne_sm"
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
mpiexec -n $NSLOTS -hosts $(hostname) nrniv -python -mpi init.py