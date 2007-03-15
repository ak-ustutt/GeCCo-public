#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
srcdir:=.

include $(topdir)/make.config

INST_PROGS:= 
INST_LIBS:=libmemman.a
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

include $(topdir)/make.rules

#mod_objs:=$(addsuffix .o,$(basename $(INST_MODULES)))
mod_objs:=memman.o
src:=$(wildcard *.f90 *.F *.f)
objs:=$(mod_objs) $(filter-out $(mod_objs),$(addsuffix .o, $(basename $(src))))

libs:=

all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)

$(libdir)/libmemman.a: $(objs)
	$(AR) rv $@ $?
	$(RAN) $@

.PHONY: dep deps
dep deps:
	$(topdir)/config/mkdep90.py $(src) >deps.mk

ifeq (, $(wildcard $(srcdir)/deps.mk))
$(shell $(topdir)/config/mkdep90.py $(src) >deps.mk)
endif
include deps.mk

