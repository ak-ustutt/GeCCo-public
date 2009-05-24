*----------------------------------------------------------------------*
      subroutine init_reorder(reo)
*----------------------------------------------------------------------*
*     initialize all pointers 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(reorder), intent(inout) ::
     &     reo

      reo%nreo = 0
      
      reo%occ_opin => null()
      reo%occ_opout => null()
      reo%rst_opin => null()
      reo%rst_opout => null()
      reo%occ_shift => null()
      reo%occ_op0 => null()
      reo%from_to => null()
      reo%merge_stp1    => null()
      reo%merge_stp1inv => null()
      reo%merge_stp2    => null()
      reo%merge_stp2inv => null()
      
      return
      end
