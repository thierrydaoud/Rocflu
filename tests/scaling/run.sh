#!/bin/bash
# Run the solver on a pre-processed scaling case.
# Usage (inside a job): run.sh <rundir> <nranks> <bindir>; SRUN_MPI sets the srun MPI plugin (default pmix).
set -e
rundir=$1; n=$2; bin=$3
cd "$rundir"
srun --mpi=${SRUN_MPI:-pmix} -n "$n" "$bin/rflump" -c rectshktb -v 1 > rflump.${SLURM_JOB_ID:-local}.log 2>&1
