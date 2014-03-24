#!/bin/bash

# REQUIRED:
# dat/bchl.gro (with hydrogens)
# dat/forcefield.itp (basic forcefield control file)
# prm2gmx.py (REQ:  dat/bchl.tpg && dat/bchl.prm)


# ============GENERATE A FORCEFIELD====================
if [ -e output.ff ]; then
    rm -r output.ff
fi
mkdir output.ff
python prm2gmx.py


# ============GENERATE THE SINGLE-MOLECULE ITP====================
cp dat/forcefield.itp output.ff/.
echo "BCL    BCL" > residuetypes.dat
pdb2gmx -ff output -f dat/bchl.gro -o output.ff/bcl.gro -p output.ff/bcl.top -i output.ff/bcl_posre.itp
cp output.ff/bcl.top .
cp output.ff/bcl.gro .

# ====================TEST THE TOPOLOGY =============================
grompp -f dat/mdp/em.mdp -c bcl.gro -p bcl.top -o trash -po trash
tail -n+21 output.ff/bcl.top | head -n-8 >> output.ff/bcl.itp
rm bcl.top
rm bcl.gro
rm trash*
rm residuetypes.dat

TOP_OK=true

if [ "$TOP_OK" = true ]; then


    if [ -e mddir ]; then
        rm -r mddir
    fi
    mkdir mddir
    
    # ============ADD HYDROGENS TO BCL FROM FMO, MOVE TO MDDIR=================
    cat dat/bchl.hdb > output.ff/bcl.hdb
    
    cat dat/pdb/4BCL_FIX.pdb | grep BCL > mddir/4BCL_BCL.pdb 
    pdb2gmx -ff output -f mddir/4BCL_BCL.pdb -o mddir/4BCL_BCL.pdb -p trash -i trash
    rm mddir/\#*
    rm trash*
    
    
    #==============MOVE FF DATA INTO CHARMM_MOD===============================
    cp -r dat/charmm27.ff mddir/charmm_mod.ff
    cat dat/pdb/4BCL_FIX.pdb | grep ATOM > mddir/4BCL_PROTEIN.pdb 
    
    #tail -n+30 output.ff/bcl.rtp >> mddir/charmm_mod.ff/aminoacids.rtp
    #cat output.ff/bcl.rtp        > mddir/charmm_mod.ff/bcl.rtp
    #cat dat/bchl.hdb             > mddir/charmm_mod.ff/bcl.hdb
    head -n-8 output.ff/bcl.itp   > mddir/charmm_mod.ff/bcl.itp
    cat output.ff/bcl_posres.itp  > mddir/charmm_mod.ff/bcl_posres.itp
    cat output.ff/atomtypes.atp  >> mddir/charmm_mod.ff/atomtypes.atp
    
    cat output.ff/ffbonded.itp    >> mddir/charmm_mod.ff/ffbonded.itp
    cat output.ff/ffnonbonded.itp >> mddir/charmm_mod.ff/ffnonbonded.itp
    
    
    #================GENERATE PROTEIN .GRO FILE AND MERGE IN BCLs=======================
    cd mddir 
    pdb2gmx -f 4BCL_PROTEIN.pdb -o conf.pdb -ff charmm_mod -chainsep id_and_ter -water tip3p -p 4BCL.top
    cat conf.pdb | grep ATOM                     > 4BCL.pdb
    cat 4BCL_BCL.pdb | grep ATOM                >> 4BCL.pdb
    
    editconf -f 4BCL.pdb -o 4BCL.gro -d 1
    
    
    
    
    #================INSERT ITP INCLUDES INTO TOP FILE=======================
    head -n-8 4BCL.top                                    > 4BCL_FIX.top
    echo ''                                              >> 4BCL_FIX.top
    echo '; include BCL forcefield'                      >> 4BCL_FIX.top
    echo '#include "./charmm_mod.ff/bcl.itp"'            >> 4BCL_FIX.top
    echo ''                                              >> 4BCL_FIX.top
    echo '; include BCL position restraints'             >> 4BCL_FIX.top
    echo '#ifdef POSRES'                                 >> 4BCL_FIX.top
    echo '#include "./charmm_mod.ff/bcl_posres.itp"'     >> 4BCL_FIX.top
    echo '#endif'                                        >> 4BCL_FIX.top
    tail -n8 4BCL.top                                    >> 4BCL_FIX.top
    echo 'BCL                 7'                         >> 4BCL_FIX.top
    
    mv 4BCL_FIX.top 4BCL.top
    
    
    # ============RUN SIMPLE CODE====================
    
    mkdir em
    echo "BCL    Pigment" > residuetypes.dat
    grompp -v -p 4BCL.top -c 4BCL.gro -f ../dat/mdp/em.mdp -o em/em -po em/em
   
fi 
