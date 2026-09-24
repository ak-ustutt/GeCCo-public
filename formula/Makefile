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

# factorize_new.f is miscompiled by gfortran (15.2.0, arm64) at -O2: the
# idxsortx calls at the end of the routine crash (wild-pointer SIGSEGV
# inside idxsortx_qs_dsc, called on time_stat/mem_stat/scale_stat); not
# reproducible with -O0, and bisected away from sortrout.f itself.
# (split_contr4.f used to be listed here too, but its crash turned out to
# be a source bug -- an unassociated arc pointer passed to an
# explicit-shape dummy -- now fixed there.)
llo_obj_list:=factorize_new.o

include $(topdir)/make.rules

all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)

