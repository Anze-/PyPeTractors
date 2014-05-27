#!/bin/sh
#PBS -V
#PBS -N mag_pendulum_job
#PBS -l nodes=120
#PBS -l alberto.anzellotti=10:00:00
#PBS -m bea
#PBS -M alberto.anzellotti@science.unitn.it
cd $PBS_O_WORKDIR
/usr/local/bin/mpirun -np 120 /home/alberto.anzellotti/mpi1.o
python collector.py
