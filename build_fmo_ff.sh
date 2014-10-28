#!/bin/bash

#------------------------------------------------------------
#               REQUIRED SOFTWARE
#------------------------------------------------------------
# GROMACS 4 installation
#   pdb2gmx
# prm2gmx.py            -- Converts .tpg and .prm AMBER94 files to gromacs

#------------------------------------------------------------
#                REQUIRED FILES
#------------------------------------------------------------
# dat/bchl.tpg          -- AMBER94 forcefield params for BChl
# dat/bchl.prm          -- AMBER94 forcefield params for BChl
# dat/bchl.gro          -- One BChl molecule
# dat/forcefield.itp    -- Forcefield header from GROMACS
# dat/pdb/4BCL_FIX.pdb  -- 4BCL pdb file 
# dat/bcl.hdb           -- BChl hydrogen-file for GROMCAS, used for 4BCL pdb
# dat/amber99sb-ildn.ff -- Amber99 forcefield for GROMACS


if [ -e fmo.ff ]; then
    rm -r fmo.ff
fi
mkdir fmo.ff

# ============ADD HYDROGENS TO FMO BCL's, MOVE TO MDDIR=================
cat dat/bchl.hdb > bcl.ff/bcl.hdb

# Rip BCLs out of pdb file
cat dat/pdb/4BCL_FIX.pdb | grep BCL > fmo.ff/4BCL_BCL.pdb 
pdb2gmx -ff bcl -f fmo.ff/4BCL_BCL.pdb -o fmo.ff/4BCL_BCL.pdb -p trash -i trash
rm fmo.ff/\#*
rm trash*


#==============MERGE FF DATA==============================
cp -r dat/gromacs-ff/amber99sb-ildn.ff fmo.ff/amber_mod.ff
cat dat/pdb/4BCL_FIX.pdb | grep ATOM > fmo.ff/4BCL_PROTEIN.pdb 

#tail -n+30 bcl.ff/bcl.rtp >> fmo.ff/amber_mod.ff/aminoacids.rtp
#cat bcl.ff/bcl.rtp        > fmo.ff/amber_mod.ff/bcl.rtp
#cat dat/bchl.hdb             > fmo.ff/amber_mod.ff/bcl.hdb
cp bcl.ff/bcl_cdc.itp         fmo.ff/amber_mod.ff/bcl_cdc.itp
cp bcl.ff/bcl.itp             fmo.ff/amber_mod.ff/bcl.itp
cp bcl.ff/bcl_posre.itp       fmo.ff/amber_mod.ff/bcl_posre.itp
cat bcl.ff/atomtypes.atp  >>  fmo.ff/amber_mod.ff/atomtypes.atp

cat bcl.ff/ffbonded.itp    >> fmo.ff/amber_mod.ff/ffbonded.itp
cat bcl.ff/ffnonbonded.itp >> fmo.ff/amber_mod.ff/ffnonbonded.itp


#================GENERATE PROTEIN .GRO FILE AND MERGE IN BCLs=======================
cd fmo.ff 
pdb2gmx -f 4BCL_PROTEIN.pdb -o conf.pdb -ff amber_mod -chainsep id_and_ter -water tip3p -p 4BCL.top
cat conf.pdb | grep ATOM                     > 4BCL.pdb
cat 4BCL_BCL.pdb | grep ATOM                >> 4BCL.pdb

editconf -f 4BCL.pdb -o 4BCL.gro -d 1 -bt dodecahedron

#================INSERT ITP INCLUDES INTO TOP FILE=======================
head -n-8 4BCL.top                                    > 4BCL_FIX.top
echo ''                                              >> 4BCL_FIX.top
echo '; include BCL forcefield'                      >> 4BCL_FIX.top
echo '#include "./amber_mod.ff/bcl.itp"'            >> 4BCL_FIX.top
echo ''                                              >> 4BCL_FIX.top
echo '; include BCL position restraints'             >> 4BCL_FIX.top
echo '#ifdef POSRES'                                 >> 4BCL_FIX.top
echo '#include "./amber_mod.ff/bcl_posre.itp"'      >> 4BCL_FIX.top
echo '#endif'                                        >> 4BCL_FIX.top
tail -n8 4BCL.top                                    >> 4BCL_FIX.top

echo ''                                              >> 4BCL_FIX.top
echo 'BCL                 7'                         >> 4BCL_FIX.top
echo ''                                              >> 4BCL_FIX.top
mv 4BCL_FIX.top 4BCL.top

rm 4BCL_PROTEIN.pdb
rm 4BCL_BCL.pdb
rm 4BCL.pdb
rm conf.pdb


# ============RUN SIMPLE CODE====================

genbox -cp 4BCL.gro -p 4BCL.top -o 4BCL.gro -cs spc216.gro
grompp -f ../dat/mdp/ions.mdp -c 4BCL.gro -p 4BCL.top -o temp.tpr
genion -s temp.tpr -o 4BCL.gro -p 4BCL.top -pname NA -nname CL -neutral <<< "SOL"
#grompp -f ../dat/mdp/ions.mdp -c 4BCL.gro -p 4BCL.top -o temp.tpr
#trjconv -f 4BCL.gro -s temp.tpr -ur compact -o 4BCL_pbc.gro -pbc res
rm \#*\#
rm temp.tpr

echo ""
echo ""
echo ""
echo ""
echo "Calling make_ndx next!!!"
echo "INSTRUCTIONS: "
echo "1. When prompted, input '[Protein-index] | [BCL-index]' (e.g. '2 | 21')." 
echo "2. Assert that the group is called Protein_BCL"
echo "3. Then, enter 'q' to exit."
echo "Alternatively, if you forget these instructions, you may exit with Ctrl-C safely, rerun the script, and pay more attention"
echo ""
echo ""
echo ""
echo ""
read -n1 -r -p "Press any key to continue... (where's the any key?)" key

make_ndx -f 4BCL.gro -o index.ndx

mkdir em
echo "BCL    Pigment" > residuetypes.dat
grompp -v -p 4BCL.top -c 4BCL.gro -f ../dat/mdp/em.mdp -o em/em -po em/em
cd em
mdrun -v -deffnm em
cd ..

trjconv -f em/em.trr -s em/em.tpr -o em/em_vid.gro -pbc res -ur compact -n index.ndx <<< "Protein_BCL"

grompp -f em/em.mdp -c 4BCL.gro -p 4BCL.top -pp 4BCL_pp.top -o trash -po trash
rm trash*


