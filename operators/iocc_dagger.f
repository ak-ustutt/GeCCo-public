*----------------------------------------------------------------------*
*     functions to calculate with hpv/ca occupations
*----------------------------------------------------------------------*
      function iocc_dagger(iocc_in)
*----------------------------------------------------------------------*
*     return daggered occupation
*     include interface file ifc_ioccfunc.inc in calling routines!
*----------------------------------------------------------------------*
      
      implicit none
      include 'opdim.h'

      integer :: iocc_dagger(ngastp,2)

      integer, intent(in) :: iocc_in(ngastp,2)

      ! function result and argument may be the same:
      integer :: iscr(ngastp,2)

      iscr(1:ngastp,1) = iocc_in(1:ngastp,2)
      iscr(1:ngastp,2) = iocc_in(1:ngastp,1)

      iocc_dagger = iscr

      return
      end
