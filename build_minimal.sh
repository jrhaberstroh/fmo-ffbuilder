#!/bin/bash

cd mddir 
if [ -e em ]; then
    rm -r em; fi
mkdir em
grompp -v -p 4BCL.top -c 4BCL.gro -f ../dat/mdp/em.mdp -o em/em -po em/em
cd em; mdrun -v -deffnm em; cd ..
trjconv -f em/em.trr -s em/em.tpr -o em/em_vid.gro -pbc res -ur compact

