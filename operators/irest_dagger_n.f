      function irest_dagger_n(irest_in,njoined,ngas,nspin)
*----------------------------------------------------------------------*
*     return daggered restriction (for multi-vertex case)
*     include interface file ifc_ioccfunc.inc in calling routines!
*----------------------------------------------------------------------*
      
      implicit none
      include 'opdim.h'

      integer, intent(in) :: njoined,ngas,nspin,
     &     irest_in(2,ngas,2,2,nspin,njoined)

      integer :: irest_dagger_n(2,ngas,2,2,nspin,njoined)

      ! function result and argument may be the same:
      integer :: iscr(2,ngas,2,2,nspin,njoined), ijoin, ijoin_d, ispin

      do ijoin = 1, njoined
        ijoin_d = njoined+1-ijoin
        do ispin = 1, nspin
          iscr(1:2,1:ngas,1,1:2,ispin,ijoin_d)=
     &         irest_in(1:2,1:ngas,2,1:2,ispin,ijoin)
          iscr(1:2,1:ngas,2,1:2,ispin,ijoin_d)=
     &         irest_in(1:2,1:ngas,1,1:2,ispin,ijoin)
        end do
      end do

      irest_dagger_n = iscr

      return
      end
