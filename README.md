# GeCCo-public
public version of our general contraction code project

This code implements the internally contracted multireference theory as published by our group. 
In addition, it also provides arbitrary order coupled-cluster calculations (for closed-shell cases)
and explicitly correlated methods (the latter requires a special interface to a patched version of 
DALTON, please contact the main author, if interested).

Prerequisites:
Linux operating system, Python 3.x, GNU autotools, Fortran compiler (gfortran, ifort), C compiler (gcc, icc)

Quantum chemical software: Interfaces exist for DALTON, Cfour, Molpro; there was also an interface to GAMESS,
but this required a patch of GAMESS (have not looked into that for a long time). Most convenient integration 
is provided for Molpro.

Compilation - quick guide: (all actions should be carried out in the root directory of the project)

(1) Set up configure:

    > autoconf
    
    Ignore warnings, the generated configure script should normally work
    
(2) Configure: Try one of these

    > FC=ifort CC=icc ./configure --with-blas='-mkl=sequential'
    
    > FC=gfortran CC=gcc ./configure --with-blas='-lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm'
    
    I have not tried other blas libraries in a long time. Let me know about any other working solutions.
    
(3) Build:

    > make deps
    
    > make -j 8    # or as many processes as you like to spend
    
   The binary will be <root>/bin/<arch>/gecco.x where <arch> is something like x86_64-linux-gnu-gfortran
   (you can use the same source code for building several versions on different systems or compilers)
    
(4) Set environment variables (e.g. in your .bashrc or similar):
    
    export GECCO_DIR=<path to this project>
    
    export GECCO_BIN=$GECCO_DIR/bin/<arch>/gecco.x   # required for Molpro integration
    
    export GECCO_TMP=/SCR/$USER/GECCOTEST.$$  # for testing; adapt as appropriate; should point to a scratch disk
    
(5) Test:
    
    > cd test
    
    > make -j 8   # to test on 8 threads
    
    a shorter version is called by
    > make essential -j 8

(6) Use it!
    Documentation --> see Wiki
    
