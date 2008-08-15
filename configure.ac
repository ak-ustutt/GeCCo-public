#
# Process this file with autoconf to produce a configure script.
#
# $Id: configure.ac 1.7 05/05/05 22:13:58+02:00 jonas@radon $

AC_INIT([GeCCo],[0.5],[andreas.koehn@uni-mainz.de])
AC_CONFIG_SRCDIR(gecco/gecco.f)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR(config)
AC_PREFIX_DEFAULT([/usr/local])
AC_PROG_MAKE_SET
AC_CANONICAL_HOST

# Check for programs.

CHECK_GNU_MAKE
if test "$ifGNUmake" = "#"; then
	AC_MSG_ERROR([No working version of GNU make found!])
fi
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CHECK_PROGS(LD,ld,false)
AC_CHECK_PROGS(AR,ar,false)
AC_CHECK_PROGS(STRIP,strip,true)
AC_CHECK_PROGS(ETAGS,etags,true)
AC_CHECK_PROGS(CTAGS,ctags,true)
AC_PROG_LN_S

AC_CACHE_SAVE

# ACX_ macros are defined in aclocal.m4
# Setup compilers for different architechtures
case $host in 
	*i*86*-linux*)
		try_cc="gcc icc pgcc cc"
		try_fc="ifort ifc pgf90 f90"
		build_type="i386_linux"
		;;
	*ia*64*-linux*)
		try_cc="gcc icc cc"
		try_fc="ifort ifc pgf90 f90"
		build_type="i64_linux"
		;;
	*x86*64*-linux*)
		try_cc="gcc icc pgcc cc"
		try_fc="pgf90 ifort ifc f90"
		build_type="x86_64-linux"
		;;
	*alpha*-linux*)
		try_cc="ccc cc gcc"
		try_fc="fort f90"
		build_type="alpha_linux"
		;;
	*dec*-osf*)
		try_cc="ccc cc gcc"
		try_fc="fort f90"
		build_type="dec_osf"
		;;
	*ibm*-aix*)
		try_cc="cc gcc"
		try_fc="xlf90 f90"
		build_type="ibm_aix"
		;;
	*sgi*-irix*)
		try_cc="cc gcc"
		try_fc="f90"
		build_type="sgi_irix"
		;;
	*hppa*-hpux*)
		try_cc="cc gcc"
		try_fc="f90"
		build_type="hp_hpux"
		;;
	*sparc*-solaris*)
		try_cc="cc gcc"
		try_fc="f90"
		build_type="sun_solaris"
		;;
	*) 	
		AC_MSG_WARN([Unknown architecture! This might not work...])
		try_cc="cc gcc ccc cl"
		try_fc="f90 xlf90 ifc pgf90 epcf90 f95 xlf95 lf95 fort"
		build_type="unknown"
		;;
esac

AC_PROG_CC([$try_cc])
AC_PROG_FC([$try_fc])

ACX_BUILD_FLAGS([$build_type])
AC_CACHE_SAVE

AX_FC_SEARCH_PATH

AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS

#AC_PROG_YACC
#AC_PROG_LEX
#test "$YACC" = "bison -y" && YACC="$YACC -d"

dnl ACX_GETKW([], AC_MSG_ERROR([libgetkw.a not found!]))

AC_LANG_PUSH([Fortran])
ACX_BLAS(,AC_MSG_ERROR([No working blas library found!]))
ACX_LAPACK(,AC_MSG_ERROR([No working lapack library found!]))
AC_LANG_POP([Fortran])

AC_LANG_PUSH([Fortran])
AC_FC_SRCEXT([f90])
#AC_FC_FREEFORM
AC_LANG_POP([Fortran])

AX_F90_MODULE_EXTENSION
if test "x$ax_f90_modext" = "xunknown" ; then
     AC_MSG_ERROR([unable to find f90 modules extension])
fi

AC_CACHE_SAVE
#AC_MSG_NOTICE([Checking if MPI is available])
#AC_LANG_PUSH([Fortran 90])
#ACX_MPI([FC=$MPIFC])
#AC_LANG_POP
#fdefs="$fdefs -DHAVE_MPI"
#AC_CACHE_SAVE

ACX_SUBST_BUILD_FLAGS 
AC_SUBST(top_bindir,['$(topdir)/bin/$(ARCH)'])
AC_SUBST(top_includedir,['$(topdir)/include'])
AC_SUBST(top_libdir,['$(topdir)/lib/$(ARCH)'])

AC_CONFIG_FILES([Makefile make.config make.rules])
AC_OUTPUT
