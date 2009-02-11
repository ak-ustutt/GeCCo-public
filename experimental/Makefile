#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
srcdir:=.

this_lib:=libexperimental.a

include $(topdir)/make.config

# switch call to user routines in exp_wrap.F
ifneq ($(strip $(wildcard set_my*.f)),)
  FCFLAGS+=-DEXPRIMENTAL
endif

INST_PROGS:= 
INST_LIBS:=$(this_lib)
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

src:=$(wildcard *.f90 *.F *.f)

include $(topdir)/make.rules

all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)

