#!/usr/bin/env bash
set -euo pipefail

execu=${EXECUTABLE:-./a.out}
np=${MPI_PROCS:-8}
prefix=${OUT_BASENAME:-test_EW_SF_Ebeam_}

cat > _template.in <<'EOF'
100
256
__ENU__
2
0.01d0
225.0d0
40.0
0
6  12
0
.true.
EOF

printf "# Enu(MeV)    sigma\n" > xsec_vs_enu.txt
printf "Enu_MeV,sigma\n" > xsec_vs_enu.csv

for enu in $(seq 200 100 2000); do
  sed "s/__ENU__/${enu}.0/" _template.in | mpirun -np "$np" "$execu" > /dev/null
  outfile="${prefix}${enu}.out"
  sigma=$(head -n1 "$outfile" | awk '{print $1}')
  printf "%6d        %.16e\n" "$enu" "$sigma" >> xsec_vs_enu.txt
  echo "[scan] Enu=$enu MeV -> sigma=$sigma"
done
