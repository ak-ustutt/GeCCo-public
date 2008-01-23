      function irest_dagger_n(irest_in,njoined,ngas)
*----------------------------------------------------------------------*
*     return daggered restriction (for multi-vertex case)
*     include interface file ifc_ioccfunc.inc in calling routines!
*----------------------------------------------------------------------*
      
      implicit none
      include 'opdim.h'

      integer, intent(in) :: njoined,ngas,irest_in(2,ngas,2,2,njoined)

      integer :: irest_dagger_n(2,ngas,2,2,njoined)

      ! function result and argument may be the same:
      integer :: iscr(2,ngas,2,2,njoined), ijoin, ijoin_d

      do ijoin = 1, njoined
        ijoin_d = njoined+1-ijoin
        iscr(1:2,1:ngas,1,1:2,ijoin_d)=irest_in(1:2,1:ngas,2,1:2,ijoin)
        iscr(1:2,1:ngas,2,1:2,ijoin_d)=irest_in(1:2,1:ngas,1,1:2,ijoin)
      end do

      irest_dagger_n = iscr

      return
      end
