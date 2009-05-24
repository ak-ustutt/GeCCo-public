*----------------------------------------------------------------------*
      subroutine dealloc_bcontr(bcontr)
*----------------------------------------------------------------------*
*     deallocate all subfields in binary contraction bcontr
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'def_contraction.h'

      type(binary_contr), intent(inout) ::
     &     bcontr

      if (associated(bcontr%occ_res)) deallocate(bcontr%occ_res)
      if (associated(bcontr%occ_op1)) deallocate(bcontr%occ_op1)
      if (associated(bcontr%occ_op2)) deallocate(bcontr%occ_op2)
      if (associated(bcontr%rst_res)) deallocate(bcontr%rst_res)
      if (associated(bcontr%rst_op1)) deallocate(bcontr%rst_op1)
      if (associated(bcontr%rst_op2)) deallocate(bcontr%rst_op2)
      if (associated(bcontr%occ_ex1)) deallocate(bcontr%occ_ex1)
      if (associated(bcontr%occ_ex2)) deallocate(bcontr%occ_ex2)
      if (associated(bcontr%occ_cnt)) deallocate(bcontr%occ_cnt)
      if (associated(bcontr%merge_op1)) deallocate(bcontr%merge_op1)
      if (associated(bcontr%merge_op2)) deallocate(bcontr%merge_op2)
      if (associated(bcontr%merge_op1op2))
     &     deallocate(bcontr%merge_op1op2)
      if (associated(bcontr%merge_op2op1))
     &     deallocate(bcontr%merge_op2op1)

      return
      end
