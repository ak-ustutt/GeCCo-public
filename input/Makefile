#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
srcdir:=.

this_lib:=libinput.a

include $(topdir)/make.config

INST_PROGS:= 
INST_LIBS:=$(this_lib)
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

src:=$(wildcard *.f90 *.F *.f)

mod_obj_list:=parse_input.o

# set_orbinf.f is miscompiled by gfortran at -O2 (arm64): corrupts the
# heap during orbital-info setup, causing a SIGSEGV/malloc-heap-corruption
# crash later (e.g. in the very next allocation, inside put_orbinfo.f's
# I/O). Not reproducible with -O0. Same family of bug as split_contr4.f
# and factorize_new.f (see formula/Makefile).
llo_obj_list:=set_orbinf.o

include $(topdir)/make.rules

all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)

