#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
srcdir:=.

this_lib:=libcrout.a

include $(topdir)/make.config

INST_PROGS:= 
INST_LIBS:= $(this_lib)
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

#src:=$(wildcard *.c)

include $(topdir)/make.rules

#src:=$(wildcard *.c)
# objs:=$(addprefix $(ARCH)/,$(addsuffix .o, $(basename $(src))))
#objs:=$(addprefix $(ARCH)/,$(patsubst %.c,%.o,$(src)))

libs:=

all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)

#$(libdir)/libcrout.a: objdir $(all_objs)
#	$(AR) rv $@ $(all_objs)
#	$(RAN) $@


