#!/bin/bash

# ============GENERATE A FORCEFIELD====================
if [ -e output.ff ]; then
    rm -r output.ff
fi
mkdir output.ff
python prm2gmx.py


# ============GENERATE AN ITP====================
#echo "BCL    BCL" > residuetypes.dat
#cp dat/forcefield.itp output.ff/.
#pdb2gmx -ff output -f dat/bchl.gro -o output.ff/bcl.gro -p output.ff/bcl.top -i output.ff/bcl_posre.itp
#tail -n+21 output.ff/bcl.top | head -n-8 >> output.ff/bcl.itp
#rm residuetypes.dat


# ============BUILD FMO TOPOLOGY====================
if [ -e mddir ]; then
    rm -r mddir
fi
mkdir mddir
cp -r dat/amber99sb-ildn.ff mddir/ambermod.ff

tail -n+30 output.ff/bcl.rtp >> mddir/ambermod.ff/aminoacids.rtp
cat output.ff/atomtypes.atp >> mddir/ambermod.ff/atomtypes.atp
cat output.ff/ffbonded.atp >> mddir/ambermod.ff/ffbonded.atp
cat dat/bchl.hdb >> mddir/ambermod.ff/aminoacids.hdb
cd mddir
pdb2gmx -f ../dat/pdb/4BCL_FIX.pdb -ff ambermod -chainsep id_and_ter -water tip3p
