#!/bin/bash
set -o errexit


#------------------------------------------------------------
#               REQUIRED SOFTWARE
#------------------------------------------------------------
# GROMACS 4 installation
#   pdb2gmx
PDB2GMX=pdb2gmx
GROMPP=grompp
# prm2gmx.py            -- Converts .tpg and .prm AMBER94 files to gromacs

#------------------------------------------------------------
#                REQUIRED FILES
#------------------------------------------------------------
# dat/bchl.tpg          -- AMBER94 forcefield params for BChl
# dat/bchl.prm          -- AMBER94 forcefield params for BChl
# dat/bchl.gro          -- One BChl molecule atomic positions
# dat/bchl_cdc.txt      -- CDC charges for BChl
# dat/forcefield.itp    -- Forcefield header from GROMACS
# dat/pdb/4BCL_FIX.pdb  -- 4BCL pdb file 
# dat/bcl.hdb           -- BChl hydrogen-file for GROMCAS, used for 4BCL pdb
# dat/amber99sb-ildn.ff -- Amber99 forcefield for GROMACS


# ============GENERATE A FORCEFIELD====================
if [ -e bcl.ff ]; then
    rm -r bcl.ff
fi
mkdir bcl.ff
python prm2gmx.py -AMBER94prm dat/BCHL.prm -AMBER94tpg dat/BCHL.tpg \
    -GMXbonded bcl.ff/ffbonded.itp -GMXnonbonded bcl.ff/ffnonbonded.itp \
    -GMXrtp bcl.ff/bcl.rtp         -GMXatomtypes bcl.ff/atomtypes.atp \
    -suffix QL

# ============GENERATE THE ITP==========================
# Generate the single-molecule ITP for the base pararmeters
cp dat/forcefield.itp bcl.ff/.
echo "BCL    BCL" > residuetypes.dat
# Create the .top file
$PDB2GMX -ff bcl -f dat/bchl.gro -o bcl.ff/bcl.gro -p bcl.ff/bcl.top -i bcl.ff/bcl_posre.itp
# Create the .itp file by cutting the top and bottom off of the .top
tail -n+21 bcl.ff/bcl.top | head -n-12 >> bcl.ff/bcl.itp
cp bcl.ff/bcl.top bcl.ff/bcl_cdc.top

# ============GENERATE THE CDC ITP==========================

# Zero all of the hydrogen charges for the CDC file
cp bcl.ff/bcl.top bcl.ff/bcl_cdc.top
out=""
# Find all lines in the topology file of atomname H*** in BCL (appx regex "BCL \s+ H"):
#   Set their charge to 0.000
sed -i "/^.*BCL\s\+H\S*[A-Z].*/ s/\(^.*BCL\s\+\S\+\s\+[0-9]\+\s\+\)\(\S\+[0-9]\s\+\)\(\S\+[0-9]\).*;.*/\1 0.000 \t\3 ;/" bcl.ff/bcl_cdc.top
# Extra regex for remaining hydrogens, between lines of 110 and 121
#   Set their charge to 0.000
sed -i '110,121 s/\S\+[0-9]\s\+\(\S\+[0-9]\)\s\+;.*/\t0.000  \t\1\t ;/' bcl.ff/bcl_cdc.top
sed -i "/BCL/ s/BCL/BCX/" bcl.ff/bcl_cdc.top

# Insert the CDC charges with sed
while read p ; do
    atomname=$(echo $p | cut -d" " -f2 )
    q_gd=$(echo $p | cut -d" " -f7)
    q_ex=$(echo $p | cut -d" " -f8)
    # Modify all lines with the strings "BCL" "$atomname" on them with the parameters from BCHL_charges.txt
    # i.e. Create the GROMACS forcefield with the CDC charges 
    sed -i "/^.*BCL .* $atomname / s/\(^.*\s\)\(\S\+\)\(.*[0-9]\+.*BCL\s\+\S\+\s\+[0-9]\+\s\+\)\(\S\+[0-9]\s\+\)\(\S\+[0-9]\).*;.*/\1\2\3$q_gd \t\5 \t\2 \t$q_ex\t\5 ;/" bcl.ff/bcl_cdc.top
done < dat/bchl_cdc.txt
# Uncomment to use BCX name instead of BCL
#sed -i 's/BCL/BCX/g' bcl.ff/bcl_cdc.top

tail -n+21 bcl.ff/bcl_cdc.top | head -n-12 >> bcl.ff/bcl_cdc.itp
rm bcl.ff/bcl_cdc.top

# ====================TEST THE TOPOLOGY =============================
cp bcl.ff/bcl.top .
cp bcl.ff/bcl.gro .
$GROMPP -f dat/mdp/em.mdp -c bcl.gro -p bcl.top -o trash -po trash
rm trash*
rm residuetypes.dat
rm bcl.top
rm bcl.gro


