      function iocc_dagger_n(iocc_in,njoined)
*----------------------------------------------------------------------*
*     return daggered occupation (for multi-vertex case)
*     include interface file ifc_ioccfunc.inc in calling routines!
*----------------------------------------------------------------------*
      
      implicit none
      include 'opdim.h'

      integer, intent(in) :: njoined,iocc_in(ngastp,2,njoined)

      integer :: iocc_dagger_n(ngastp,2,njoined)

      ! function result and argument may be the same:
      integer :: iscr(ngastp,2,njoined), ijoin, ijoin_d

      do ijoin = 1, njoined
        ijoin_d = njoined+1-ijoin
        iscr(1:ngastp,1,ijoin_d) = iocc_in(1:ngastp,2,ijoin)
        iscr(1:ngastp,2,ijoin_d) = iocc_in(1:ngastp,1,ijoin)
      end do

      iocc_dagger_n = iscr

      return
      end
