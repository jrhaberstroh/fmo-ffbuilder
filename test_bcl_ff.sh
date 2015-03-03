#!/bin/bash

if [ $1 = "INIT" ]; then

    if [ -e bcl_test ]; then
        rm -r bcl_test
    fi
    
    mkdir bcl_test
    cp bcl.ff/bcl.gro bcl_test
    cp bcl.ff/bcl.top bcl_test
    cp -r bcl.ff bcl_test
fi

if [ $1 = "ANNEAL" ]; then
    # Anneal
    grompp_d -f dat/mdp/bcl_anneal.mdp -c bcl_test/bcl.gro -po bcl_test/bcl_anneal -o bcl_test/bcl_anneal -p bcl_test/bcl.top 
    if [ $? -eq 0 ]; then
        cd bcl_test
        mdrun_d -v -deffnm bcl_anneal
        cd ..
    fi
fi

if [ $1 = "MIN" ]; then
    # Minimize
    grompp_d -f dat/mdp/bcl_em.mdp -c bcl_test/bcl_anneal.gro -po bcl_test/bcl_em -o bcl_test/bcl_em -p bcl_test/bcl.top
    grompp=$?
    echo $grompp
    if [ $grompp -eq 0 ]; then
        cd bcl_test
        mdrun_d -v -deffnm bcl_em
        cd ..
    fi
fi

if [ $1 = "EIG" ]; then
    # Diagonalize
    grompp_d -t bcl_test/bcl_em.trr -c bcl_test/bcl_em -f dat/mdp/bcl_nma.mdp -po bcl_test/bcl_nma -o bcl_test/bcl_nma -p bcl_test/bcl.top 
    cd bcl_test
    mdrun_d -v -deffnm bcl_nma
    g_nmeig_d -f bcl_nma.mtx -s bcl_nma.tpr -os -last $(echo "140 * 3" | bc)
    trjconv_d -f eigenvec.trr -o eigenvec.gro -s bcl_nma.tpr
    cd ..
    tail -n420 bcl_test/eigenval.xvg > bcl_test/eigenval.txt
    python eig_test.py
fi
