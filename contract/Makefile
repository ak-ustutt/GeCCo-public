#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
srcdir:=.

this_lib:=libcontract.a

include $(topdir)/make.config

INST_PROGS:=
INST_LIBS:=$(this_lib)
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

src:=$(wildcard *.f90 *.F *.f)

hlo_obj_list:=contr_blk1blk2_blocked_mm.o contr_blk1blk2_wmaps_c.o contr_op1op2.o contr_op1op2_wmaps_c.o invert.o invsqrt.o

include $(topdir)/make.rules

all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)

