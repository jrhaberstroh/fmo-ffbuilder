#!/bin/bash
set -o errexit
set -o nounset

SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
FFBASE=${FFBASE=$SRCDIR/dat/gromacs-ff/amber99sb-ildn.ff}
RUNMODE=${1--h}
SOLVATE=${SOLVATE-true}
OUTDIR=${OUTDIR=$SRCDIR/FMO_conf}
OUTDIR=$(realpath $OUTDIR)

if [[ "$RUNMODE" == "-h" ]]; then
    echo "----------------------"
    echo "build_fmo_ff.sh"
    echo "----------------------"
    echo "Pass \$1 as BUILD, EQUIL, or EQUIL2."
    echo "Preprocessed topology will be generated during EQUIL."
    echo ""
    echo "Options: SOLVATE  - ($SOLVATE) generate solvent around the FMO complex"
    echo "Options: OUTDIR - ($OUTDIR) location to write FMO forcefield"
    echo "Options: FFBASE - ($FFBASE) location of base FMO forcefield"
    exit 0
fi



#------------------------------------------------------------
#               REQUIRED SOFTWARE
#------------------------------------------------------------
# GROMACS 4 installation
#   pdb2gmx
#   grompp

#------------------------------------------------------------
#                REQUIRED FILES
#------------------------------------------------------------
# dat/bchl.tpg          -- AMBER94 forcefield params for BChl
# dat/bchl.prm          -- AMBER94 forcefield params for BChl
# dat/bchl.gro          -- One BChl molecule
# dat/pdb/4BCL_FIX.pdb  -- 4BCL pdb file 
# dat/bcl.hdb           -- BChl hydrogen-file for GROMCAS, used for 4BCL pdb
# dat/forcefield.itp    -- Forcefield header from GROMACS matching forcefield used
# dat/gromacs-ff/amber99sb-ildn.ff -- Amber99 forcefield for GROMACS

if [ "$RUNMODE" = "BUILD" ]; then
    
    if [ -e $OUTDIR ]; then
        rm -r $OUTDIR
    fi
    mkdir $OUTDIR
    
    # ============ADD HYDROGENS TO FMO BCL's, MOVE TO MDDIR=================
    cat $SRCDIR/dat/bchl.hdb > $SRCDIR/bcl.ff/bcl.hdb
    
    # Rip BCLs out of pdb file
    cat $SRCDIR/dat/pdb/4BCL_FIX.pdb | grep BCL > $OUTDIR/4BCL_BCL.pdb 
    cd $SRCDIR
    pdb2gmx -ff bcl -f $OUTDIR/4BCL_BCL.pdb -o $OUTDIR/4BCL_BCL.pdb -p trash -i trash
    rm $OUTDIR/\#*
    rm trash*
    cd -
    
    #==============MERGE FF DATA==============================
    FFNEW=$(basename ${FFBASE%.*})
    FFNEW="$FFNEW"_mod.ff
    cp -r $FFBASE $OUTDIR/$FFNEW
    cat $SRCDIR/dat/pdb/4BCL_FIX.pdb | grep ATOM > $OUTDIR/4BCL_PROTEIN.pdb 
    
    #tail -n+30 $SRCDIR/bcl.ff/bcl.rtp >> $OUTDIR/$FFNEW/aminoacids.rtp
    #cat $SRCDIR/bcl.ff/bcl.rtp        > $OUTDIR/$FFNEW/bcl.rtp
    #cat $SRCDIR/dat/bchl.hdb             > $OUTDIR/$FFNEW/bcl.hdb
    cp $SRCDIR/bcl.ff/bcl_cdc.itp         $OUTDIR/$FFNEW/bcl_cdc.itp
    cp $SRCDIR/bcl.ff/bcl.itp             $OUTDIR/$FFNEW/bcl.itp
    cp $SRCDIR/bcl.ff/bcl_posre.itp       $OUTDIR/$FFNEW/bcl_posre.itp
    cat $SRCDIR/bcl.ff/atomtypes.atp  >>  $OUTDIR/$FFNEW/atomtypes.atp
    
    cat $SRCDIR/bcl.ff/ffbonded.itp    >> $OUTDIR/$FFNEW/ffbonded.itp
    cat $SRCDIR/bcl.ff/ffnonbonded.itp >> $OUTDIR/$FFNEW/ffnonbonded.itp
    
    
    #================GENERATE PROTEIN .GRO FILE AND MERGE IN BCLs=======================
    cd $OUTDIR 
    pdb2gmx -f 4BCL_PROTEIN.pdb -o conf.pdb -ff ${FFNEW%.*} -chainsep id_and_ter -water tip3p -p 4BCL.top
    cat conf.pdb | grep ATOM                     > 4BCL.pdb
    cat 4BCL_BCL.pdb | grep ATOM                >> 4BCL.pdb
    
    editconf -f 4BCL.pdb -o 4BCL.gro -d 2 -bt dodecahedron
    
    #================INSERT ITP INCLUDES INTO TOP FILE=======================
    head -n-8 4BCL.top                                    > 4BCL_FIX.top
    echo ''                                              >> 4BCL_FIX.top
    echo '; include BCL forcefield'                      >> 4BCL_FIX.top
    echo "#include \"./$FFNEW/bcl.itp\""                 >> 4BCL_FIX.top
    echo ''                                              >> 4BCL_FIX.top
    echo '; include BCL position restraints'             >> 4BCL_FIX.top
    echo '#ifdef POSRES'                                 >> 4BCL_FIX.top
    echo "#include \"./$FFNEW/bcl_posre.itp\""           >> 4BCL_FIX.top
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
    
    if $SOLVATE; then
        genbox -cp 4BCL.gro -p 4BCL.top -o 4BCL.gro -cs spc216.gro
        grompp -f $SRCDIR/dat/mdp/ions.mdp -c 4BCL.gro -p 4BCL.top -o temp.tpr
        genion -s temp.tpr -o 4BCL.gro -p 4BCL.top \
            -pname NA -nname CL -neutral <<< "SOL"
        # grompp -f $SRCDIR/dat/mdp/ions.mdp -c 4BCL.gro -p 4BCL.top -o temp.tpr
        # trjconv -f 4BCL.gro -s temp.tpr -ur compact -o 4BCL_pbc.gro -pbc res <<< "System"
        rm \#*\#
        rm temp.tpr
    fi

    # Create the index groups for:
    #   1. Protein + BCL    
    #   2. Each individual BCL
    make_ndx -f 4BCL.gro -o temp.ndx 2> /dev/null <<Protein_BCL
q
Protein_BCL
    BCL_IND=$(cat temp.ndx | grep "\[" | grep -n "\[ BCL \]" | head -n1 | cut -f1 -d':')
    echo "Index for BCL group: "$BCL_IND 
    BCL_IND=$(( BCL_IND - 1 ))
    echo "Index for BCL group: "$BCL_IND 
    make_ndx -f 4BCL.gro -n temp.ndx -o index.ndx 2> /dev/null <<BCL_Groups
"Protein" | $BCL_IND
splitres $BCL_IND
q
BCL_Groups
    echo "Index for BCL group: "$BCL_IND 
    rm temp.ndx
    

elif [ "$RUNMODE" = "EQUIL" ]; then
    cd $OUTDIR
    set -o nounset
    if [ -e em ]; then
        rm -r em
    fi
    mkdir em
    echo "BCL    Pigment" > residuetypes.dat
    grompp -v -p 4BCL.top -c 4BCL.gro -f $SRCDIR/dat/mdp/em.mdp -o em/em -po em/em
    cd em
    mdrun -v -deffnm em
    EM_PASS=$?
    echo EM_PASS: $EM_PASS
    cd ..
    trjconv -f em/em.trr -s em/em.tpr -o em/em_vid.gro -pbc res -ur compact -n index.ndx <<< "System"
    grompp -f em/em.mdp -c 4BCL.gro -p 4BCL.top -pp 4BCL_pp.top -o trash -po trash
    rm trash*
elif [ "$RUNMODE" = "EQUIL2" ]; then 
    cd $OUTDIR
    set -o nounset
    if [ -e nvt ]; then
        rm -r nvt
    fi
    mkdir nvt
    grompp -v -p 4BCL.top -c em/em.gro -f $SRCDIR/dat/mdp/nvt.mdp -o nvt/nvt -po nvt/nvt -n index.ndx
    cd nvt
    mdrun -v -deffnm nvt
else
    echo "ERROR: Bad command line argument"
    echo "Pass BUILD or EQUIL[ /2] to run ff building or system equilibration"
    exit 1
fi
