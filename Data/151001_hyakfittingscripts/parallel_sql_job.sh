#!/bin/bash
#PBS -N test
#PBS -l nodes=1:ppn=16,mem=20gb,feature=16core 
#PBS -l walltime=48:00:00
#PBS -o /gscratch/rna/willy/151001_hyakfittingscripts/error_out_fittingscript_steadystate_100115
#PBS -j oe
#PBS -d /gscratch/rna/willy/151001_hyakfittingscripts
HYAK_SLOTS=`wc -l < $PBS_NODEFILE`
module load parallel_sql
module load anaconda_2.3
parallel_sql --sql -a parallel --exit-on-term -j $HYAK_SLOTS
