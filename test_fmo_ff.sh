#!/bin/bash
set -o errexit

FMO_DIR=FMO_conf

OUT_DIR=fmo_test
if [ -e $OUT_DIR ]; then
    if [ "$1" = "INIT" ]; then
        rm -r $OUT_DIR
    fi
fi
if [ ! -e $OUT_DIR ]; then
    mkdir $OUT_DIR
fi

set -o nounset
MDP=dat/mdp
IN_PDB=dat/4BCL/4BCL_FIX.pdb
TOP=4BCL.top
INIT_GRO=4BCL.gro
EM_NAME=em
FF=amber99sb-ildn
cp -r $FMO_DIR/amber_mod.ff fmo_test

cat dat/4BCL/4BCL_FIX.pdb | grep -v BCL | grep -v HOH > $OUT_DIR/temp_config.pdb
cd $OUT_DIR
pdb2gmx -f temp_config.pdb -o $INIT_GRO -p $TOP \
        -ff $FF -water tip3p
grompp -f ../$MDP/em.mdp -c $INIT_GRO -p $TOP -po $EM_NAME -o $EM_NAME
mdrun -v -deffnm em
trjconv -f em.trr -s em.tpr -o em_vid.gro -pbc res -ur compact <<< "System"

