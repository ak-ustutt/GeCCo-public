#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
srcdir:=.

include $(topdir)/make.config

INST_PROGS:= gecco.x
INST_LIBS:=
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

include $(topdir)/make.rules

src:=$(filter-out printversion.f,$(wildcard *.f90 *.F *.f))
objs:=$(addsuffix .o, $(basename $(src)))

libs:=-lexperimental -lsolve -lcontract -lformula -lr12_spec -lcc_spec -lformula -loperators -ltargets -linput -ltransform -lstrings -lfileutil -lmemman -lbaserout -lmath -lcrout


all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)

$(bindir)/gecco.x: $(objs) $(libs)
	$(MAKE) printversion.o
	$(FC) $(FLDFLAGS) -o $@ $+ $(libs) printversion.o $(BLAS_LIBS) 

date.h:
	./date.sh "$(FC)" "$(FCFLAGS)"

printversion.o: date.h printversion.f
	$(FC) $(FCFLAGS) -c printversion.f -o printversion.o
	rm date.h

.PHONY: dep deps
dep deps:
	$(topdir)/config/mkdep90.py $(src) >deps.mk

ifeq (, $(wildcard $(srcdir)/deps.mk))
$(shell $(topdir)/config/mkdep90.py $(src) >deps.mk)
endif
include deps.mk

