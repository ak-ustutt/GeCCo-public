*----------------------------------------------------------------------*
      subroutine dealloc_reorder(reo)
*----------------------------------------------------------------------*
*     deallocate all subfields in binary contraction bcontr
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'def_contraction.h'

      type(reorder), intent(inout) ::
     &     reo

      if (associated(reo%occ_opin))  deallocate(reo%occ_opin)
      if (associated(reo%occ_opout)) deallocate(reo%occ_opout)
      if (associated(reo%rst_opin))  deallocate(reo%rst_opin)
      if (associated(reo%rst_opout)) deallocate(reo%rst_opout)
      if (associated(reo%occ_shift)) deallocate(reo%occ_shift)
      if (associated(reo%occ_op0))   deallocate(reo%occ_op0)
      if (associated(reo%from_to))   deallocate(reo%from_to)
      if (associated(reo%merge_stp1))   deallocate(reo%merge_stp1)
      if (associated(reo%merge_stp1inv))deallocate(reo%merge_stp1inv)
      if (associated(reo%merge_stp2))   deallocate(reo%merge_stp2)
      if (associated(reo%merge_stp2inv))deallocate(reo%merge_stp2inv)

      reo%occ_opin => null()
      reo%occ_opout => null()
      reo%rst_opin => null()
      reo%rst_opout => null()
      reo%occ_shift => null()
      reo%occ_op0 => null()
      reo%from_to => null()
      reo%merge_stp1 => null()
      reo%merge_stp1inv => null()
      reo%merge_stp2 => null()
      reo%merge_stp2inv => null()

      return
      end
