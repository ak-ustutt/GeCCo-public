#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
srcdir:=.

this_lib:=libformula.a

include $(topdir)/make.config

INST_PROGS:= 
INST_LIBS:=$(this_lib)
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

src:=$(wildcard *.f90 *.F *.f)

# split_contr4.f is miscompiled by gfortran (15.2.0, arm64) at -O1/-O2:
# a null-pointer dereference in the inlined pointer-component reallocation
# from resize_contr; not reproducible with -O0. Build it unoptimized until
# this is understood / fixed upstream.
llo_obj_list:=split_contr4.o

include $(topdir)/make.rules

all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)

