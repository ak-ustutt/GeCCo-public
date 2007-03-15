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
INST_LIBS:= libcrout.a
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

include $(topdir)/make.rules

src:=$(wildcard *.c)
objs:=$(addsuffix .o, $(basename $(src)))

libs:=

all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)

$(libdir)/libcrout.a: $(objs)
	$(AR) rv $@ $?
	$(RAN) $@


