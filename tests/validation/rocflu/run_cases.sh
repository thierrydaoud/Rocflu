#!/bin/bash
# Runs the cases listed in a file, one after another, on the ranks of this allocation:
#   sbatch -n 8 -t 01:00:00 --wrap "bash $PWD/run_cases.sh <Rocflu bin directory> <root>/cases.txt"
# (--wrap keeps the script in place, so it finds the repository.) Each case: mesh, validate,
# rflumap/rflupart/rfluinit, rflump, rflupost, then result.txt in the plain columns of
# metrics/compare.py. A case with done.ok is skipped; a failed one gets the file "failed".
set -u
BIN=$(cd "$1" && pwd)
LIST=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
cd "$REPO" && source ./rocflu.sh load > /dev/null  # the machine's modules and MACHINE_SRUN_MPI (sets ROOT)
CASES=$(dirname "$LIST")
N=${SLURM_NTASKS:-1}
while read -r name; do
  [ -n "$name" ] || continue
  cd "$CASES/$name" || continue
  [ -f done.ok ] && { echo "$name done"; continue; }
  rm -f failed
  { "$REPO/tests/mesh_generator" < mesh.in > mesh.log 2>&1 && mv hass.hyb.asc rectshktb.hyb.asc \
    && (cd "$REPO" && ./rocflu.sh validate "$CASES/$name/rectshktb.inp") > validate.log 2>&1 \
    && "$BIN/rflumap" -c rectshktb -m 1 -p "$N" -r "$N" -v 1 > rflumap.log 2>&1 \
    && "$BIN/rflupart" -c rectshktb -v 1 > rflupart.log 2>&1 \
    && "$BIN/rfluinit" -c rectshktb -v 1 > rfluinit.log 2>&1 \
    && SECONDS=0 && srun --mpi="${MACHINE_SRUN_MPI:-pmix}" -n "$N" "$BIN/rflump" -c rectshktb -v 1 < /dev/null > rflump.log 2>&1 \
    && echo "$SECONDS" > solver.seconds \
    && "$BIN/rflupost" -c rectshktb -v 0 > rflupost.log 2>&1 \
    && python3 "$HERE/to_columns.py" . > result.txt \
    && touch done.ok; } || touch failed
  echo "$name $( [ -f done.ok ] && echo ok || echo FAILED )"
done < "$LIST"
