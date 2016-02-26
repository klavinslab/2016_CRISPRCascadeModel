#!/bin/bash
#PBS -N test
#PBS -l nodes=1:ppn=16,mem=20gb,feature=16core
#PBS -l walltime=48:00:00
#PBS -o /gscratch/rna/willy/HyakFittingScripts_CascadeAnalysis_20151209/error_out_fittingscript_steadystate_20151209_exp2
#PBS -j oe
#PBS -d /gscratch/rna/willy/HyakFittingScripts_CascadeAnalysis_20151209
HYAK_SLOTS=`wc -l < $PBS_NODEFILE`
module load parallel_sql
module load anaconda_2.3
parallel_sql --sql -a parallel --exit-on-term -j $HYAK_SLOTS