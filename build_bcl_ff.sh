#!/bin/bash
set -o errexit
set -o nounset


#------------------------------------------------------------
#               REQUIRED SOFTWARE
#------------------------------------------------------------
# GROMACS 4 installation
#   pdb2gmx
#   grompp
# prm2gmx.py            -- Converts .tpg and .prm AMBER94 files to gromacs

PDB2GMX=pdb2gmx
GROMPP=grompp
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

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


# ============GENERATE A FORCEFIELD====================
if [ -e $SRCDIR/bcl.ff ]; then
    rm -r $SRCDIR/bcl.ff
fi
mkdir $SRCDIR/bcl.ff
cd $SRCDIR
python $SRCDIR/prm2gmx.py -AMBER94prm $SRCDIR/dat/BCHL.prm \
    -AMBER94tpg $SRCDIR/dat/BCHL.tpg \
    -GMXbonded    $SRCDIR/bcl.ff/ffbonded.itp \
    -GMXnonbonded $SRCDIR/bcl.ff/ffnonbonded.itp \
    -GMXrtp       $SRCDIR/bcl.ff/bcl.rtp \
    -GMXatomtypes $SRCDIR/bcl.ff/atomtypes.atp \
    -suffix QL -bond_scale 2.0 -angle_scale 2.0 -dihedral_scale .5
cd -

# ============GENERATE THE ITP==========================
# Generate the single-molecule ITP for the base pararmeters
cp $SRCDIR/dat/forcefield.itp $SRCDIR/bcl.ff/.
echo "BCL    BCL" > $SRCDIR/residuetypes.dat
# Create the .top file
cd $SRCDIR
$PDB2GMX -ff bcl -f $SRCDIR/dat/bchl.gro -o $SRCDIR/bcl.ff/bcl.gro \
                -p $SRCDIR/bcl.ff/bcl.top -i $SRCDIR/bcl.ff/bcl_posre.itp
cd -
# Create the .itp file by cutting the top and bottom off of the .top
tail -n+21 $SRCDIR/bcl.ff/bcl.top | head -n-12 >> $SRCDIR/bcl.ff/bcl.itp
cp $SRCDIR/bcl.ff/bcl.top $SRCDIR/bcl.ff/bcl_cdc.top

# ============GENERATE THE CDC ITP==========================

# Zero all of the hydrogen charges for the CDC file
cp $SRCDIR/bcl.ff/bcl.top $SRCDIR/bcl.ff/bcl_cdc.top
# Find all lines in the topology file of atomname H*** in BCL (appx regex "BCL \s+ H"):
#   Set their charge to 0.000
sed -i "/^.*BCL\s\+H\S*[A-Z].*/ s/\(^.*BCL\s\+\S\+\s\+[0-9]\+\s\+\)\(\S\+[0-9]\s\+\)\(\S\+[0-9]\).*;.*/\1 0.000 \t\3 ;/" \
        $SRCDIR/bcl.ff/bcl_cdc.top
sed -i "/BCL/ s/BCL/BCX/" $SRCDIR/bcl.ff/bcl_cdc.top
# Extra regex for remaining hydrogens, between lines of 110 and 121
#   Set their charge to 0.000 (count 3 cols past BCX)
sed -i '110,121 s/\(BCX\s\+\S\+\s\+\S\+\s\+\)\S\+/\1 0.000/' \
        $SRCDIR/bcl.ff/bcl_cdc.itp

# Insert the CDC charges with sed
while read p ; do
    atomname=$(echo $p | cut -d" " -f2 )
    q_gd=$(echo $p | cut -d" " -f7)
    q_ex=$(echo $p | cut -d" " -f8)
    # Modify all lines with the strings "BCL" "$atomname" on them with the parameters from BCHL_charges.txt
    # i.e. Create the GROMACS forcefield with the CDC charges 
    sed -i "/^.*BCL .* $atomname / s/\(^.*\s\)\(\S\+\)\(.*[0-9]\+.*BCL\s\+\S\+\s\+[0-9]\+\s\+\)\(\S\+[0-9]\s\+\)\(\S\+[0-9]\).*;.*/\1\2\3$q_gd \t\5 \t\2 \t$q_ex\t\5 ;/" $SRCDIR/bcl.ff/bcl_cdc.top
done < $SRCDIR/dat/bchl_cdc.txt

# Extra regex for remaining hydrogens, between lines of 110 and 121
#   Set their charge to 0.000 (count 3 cols past BCX)
sed -i '102,147 s/\(BCX\s\+\S\+\s\+\S\+\s\+\)\S\+/\1 0.000/' \
        $SRCDIR/bcl.ff/bcl_cdc.itp

tail -n+21 $SRCDIR/bcl.ff/bcl_cdc.top | head -n-12 >> $SRCDIR/bcl.ff/bcl_cdc.itp
rm $SRCDIR/bcl.ff/bcl_cdc.top

# ====================TEST THE TOPOLOGY =============================
cp $SRCDIR/bcl.ff/bcl.gro $SRCDIR/.
cp $SRCDIR/bcl.ff/bcl.top $SRCDIR/.
cd $SRCDIR
$GROMPP -f $SRCDIR/dat/mdp/em.mdp -c $SRCDIR/bcl.gro -p $SRCDIR/bcl.top -o $SRCDIR/trash -po $SRCDIR/trash
cd -
rm $SRCDIR/bcl.gro
rm $SRCDIR/bcl.top
rm $SRCDIR/trash*
rm $SRCDIR/residuetypes.dat
