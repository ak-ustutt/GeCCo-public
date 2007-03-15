dnl @synopsis AX_F90_INTERNAL_HEADMOD(MESSAGE, FILE-REGEXP, FLAG, FUNCTION-BODY, OUTPUT-VAR[, SEARCH-PATH [, ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl
dnl Internal macro used by AX_F90_HEADER and AC_F90_MODULE.
dnl
dnl @category Fortran
dnl @author Luc Maisonobe <luc@spaceroots.org>
dnl @version 2005-01-14
dnl @license AllPermissive

AC_DEFUN([AX_F90_INTERNAL_HEADMOD],[
AS_VAR_PUSHDEF([ax_include],[ax_f90_headmod_$2])
AC_MSG_CHECKING([$1])
AC_LANG_PUSH(Fortran)
AS_VAR_SET(ax_include,"not found")

if test "x$6" = x ; then
  ax_search=""
else
  ax_search="$6"
fi

if test "$prefix" = "NONE"; then
  ax_search="$ax_search:$ac_default_prefix"
else
  ax_search="$ax_search:$prefix:$ac_default_prefix"
fi

for ax_base in "" `echo $LD_LIBRARY_PATH | tr ':' '\012'` ; do
  if test "x$ax_base" != x ; then
    changequote(,)dnl
    ax_base=`echo $ax_base | sed 's,/[^/]*$,,'`
    changequote([,])dnl
    ax_search="${ax_search}:${ax_base}"
  fi
done
for ax_base in `echo $ax_search | tr ':' '\012'` ; do
 if test "AS_VAR_GET(ax_include)" = "not found" ; then
   for ax_mod in "" `find $ax_base -follow -name $2 -print` ; do
     if test "x$ax_mod" != x ; then
       changequote(,)dnl
       ax_dir=`echo $ax_mod | sed 's,/[^/]*$,,'`
       changequote([,])dnl
       ax_save_FCFLAGS="$FCFLAGS"
       FCFLAGS="$ax_save_FCFLAGS $3$ax_dir"
       AC_COMPILE_IFELSE([subroutine conftest_routine
$4
          end subroutine conftest_routine
         ],AS_VAR_SET(ax_include,"$3$ax_dir"),[])
       FCFLAGS="$ax_save_FCFLAGS"
     fi
   done
 fi
done
AC_LANG_POP(Fortran)
AC_MSG_RESULT([AS_VAR_GET(ax_include)])
if test "AS_VAR_GET(ax_include)" = "not found"; then
 $5=""
 $8
else
 $5="AS_VAR_GET(ax_include)"
 $7
fi
AC_SUBST($5)
AS_VAR_POPDEF([ax_include])
])
