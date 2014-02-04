#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
srcdir:=.

this_lib:=

include $(topdir)/make.config

INST_PROGS:= gecco.x
INST_LIBS:=
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

src:=$(wildcard *.f90 *.F *.f)

spc_obj_list:=printversion.o

include $(topdir)/make.rules

libs:=-lexperimental -lsolve -lcontract -lformula -lr12_spec -lcc_spec -lformula -loperators -ltargets -linput -ltransform -lstrings -lfileutil -lmemman -lbaserout -lmath -lcrout


#all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)
all: $(INST_PROGS)

$(bindir)/gecco.x: $(objs) $(libs)
	$(MAKE) $(ARCH)/printversion.o
	$(FC) $(FLDFLAGS) -o $@ $+ $(libs) $(ARCH)/printversion.o $(BLAS_LIBS) $(LAPACK_LIBS)

date.h:
	./date.sh "$(FC)" "$(FCFLAGS)"

$(ARCH)/printversion.o: date.h printversion.f
	$(FC) $(FCFLAGS) -c printversion.f -o $(ARCH)/printversion.o
	rm date.h

