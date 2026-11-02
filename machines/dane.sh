# Dane (LLNL). Sourced by rocflu.sh. openmpi has no Slurm PMIx support: pmi2.
# No burst QOS: standby jobs are cancelled when a normal job needs the nodes.
MACHINE_MODULES="StdEnv cmake/3.30.5 gcc/13.3.1 openmpi/4.1.2"
MACHINE_SRUN_MPI=pmi2
MACHINE_ACCOUNT=ufl
MACHINE_PARTITION=pbatch
MACHINE_TEST_PARTITION=pdebug
MACHINE_QOS_PRIMARY=normal
MACHINE_QOS_BURST=
