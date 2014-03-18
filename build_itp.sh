#!/bin/bash

echo "BCL    Ligand" > residuetypes.dat
python prm2gmx.py
cp dat/forcefield.itp output.ff/.
pdb2gmx -ff output -f dat/bchl.gro -o output.ff/bcl.gro -p output.ff/bcl.itp
rm residuetypes.dat
