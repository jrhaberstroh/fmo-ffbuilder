#!/bin/bash

if [ -e bcl.ff ]; then
    rm -r bcl.ff
fi
if [ -e fmo.ff ]; then
    rm -r fmo.ff
fi
if [ -e bcl_test ]; then
    rm -r bcl_test
fi

if [ "$1" = "TEST" ]; then
    python test__prm2gmx.py
    if [ $? = 0 ]; then
        read -p 'Press [Enter] key to continue...'
    else
        exit 101
    fi
fi
./build_bcl_ff.sh
if [ "$1" = "TEST" ]; then
    ./test_bcl_ff.sh INIT
    ./test_bcl_ff.sh ANNEAL
    ./test_bcl_ff.sh MIN
    ./test_bcl_ff.sh EIG
    read -p 'Press [Enter] key to continue...'
fi
./build_fmo_ff.sh BUILD
./build_fmo_ff.sh EQUIL



