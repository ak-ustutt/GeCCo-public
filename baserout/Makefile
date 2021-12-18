#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
srcdir:=.

this_lib:=libbaserout.a

include $(topdir)/make.config

INST_PROGS:= 
INST_LIBS:=$(this_lib)
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

hlo_obj_list:=i8common_entry.o i8elsum.o i8list_cmp.o i8list_ovl.o i8mat_cmp.o i8mltlist.o zero_i8vec.o zero_ivec.o

src:=$(wildcard *.f90 *.F *.f)

include $(topdir)/make.rules

all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)


