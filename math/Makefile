#
# $Id: Makefile.in 1.3 05/05/03 19:18:55+02:00 jonas@radon $
#

#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
srcdir:=.

this_lib:=libmath.a

include $(topdir)/make.config

INST_PROGS:=
INST_LIBS:= $(this_lib)
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

src:=$(wildcard *.f90 *.F *.f)

mod_obj_list:=
spc_obj_list:=
hlo_obj_list:=check_diag.o check_inv.o dmdiagm.o dvdsbsp_assemble_residual.o dvdsbsp_contract.o dvdsbsp_del.o dvdsbsp_diagonalize.o dvdsbsp_init.o dvdsbsp_orthvec.o dvdsbsp_update.o dvdsbsp_util.o eigen_asym.o eigen_asym_met.o eispack.o gaussj.o invsqrt_mat.o jacobi.o linpack.o mat_svd_traf.o mgs.o regular_tikhonov.o spinsym_traf.o svd_drv.o svd_get_left.o vecsp_assemble_vector.o vecsp_del.o vecsp_get_list.o vecsp_init.o vecsp_orthvec.o
llo_obj_list:=

libs:=

include $(topdir)/make.rules

all: $(INST_PROGS) $(INST_LIBS) $(INST_MODULES)


