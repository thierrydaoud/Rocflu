#!/bin/bash
# Pre-process a scaling case: mesh, rflumap, rflupart, rfluinit, one process each.
# Usage (inside a job): prep.sh <rundir> <nregions> <bindir>
# Each step's wall time and peak memory go to <step>.time (/usr/bin/time -v).
set -e
rundir=$1; n=$2; bin=$3
here=${SCALING_DIR:?set SCALING_DIR to tests/scaling}
cd "$rundir"
/usr/bin/time -v -o mesh.time "$here/../mesh_generator" < mesh.in > mesh.log 2>&1
mv hass.hyb.asc rectshktb.hyb.asc
/usr/bin/time -v -o rflumap.time  "$bin/rflumap"  -c rectshktb -m 1 -p "$n" -r "$n" -v 1 > rflumap.log 2>&1
/usr/bin/time -v -o rflupart.time "$bin/rflupart" -c rectshktb -v 1 > rflupart.log 2>&1
/usr/bin/time -v -o rfluinit.time "$bin/rfluinit" -c rectshktb -v 1 > rfluinit.log 2>&1
grep -l "^ *ERROR" rflumap.log rflupart.log rfluinit.log && exit 1
for s in mesh rflumap rflupart rfluinit; do
  printf '%-9s %s  %s\n' $s "$(grep 'Elapsed (wall' $s.time | awk '{print $NF}')" \
    "$(grep 'Maximum resident' $s.time | awk '{printf "%.2f GB", $NF/1048576}')"
done
