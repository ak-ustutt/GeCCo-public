*----------------------------------------------------------------------*
      integer function nca_i8occ(i8occ)
*----------------------------------------------------------------------*
*     return the number of C and A of i8-stored occupation
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'

      integer(8), intent(in) ::
     &     i8occ
      integer(8) ::
     &     intwk
      
      nca_i8occ = 0
      
      intwk = i8occ
      do while (intwk.gt.0)
        nca_i8occ = nca_i8occ + mod(intwk,pack_base)
        intwk = intwk/pack_base
      end do

      return
      end
