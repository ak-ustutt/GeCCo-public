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
      if (associated(bcontr%rst_cnt)) deallocate(bcontr%rst_cnt)
      if (associated(bcontr%rst_ex1)) deallocate(bcontr%rst_ex1)
      if (associated(bcontr%rst_ex2)) deallocate(bcontr%rst_ex2)
      if (associated(bcontr%merge_op1)) deallocate(bcontr%merge_op1)
      if (associated(bcontr%merge_op2)) deallocate(bcontr%merge_op2)
      if (associated(bcontr%merge_op1op2))
     &     deallocate(bcontr%merge_op1op2)
      if (associated(bcontr%merge_op2op1))
     &     deallocate(bcontr%merge_op2op1)
      if (associated(bcontr%itf_index_info))
     &     deallocate(bcontr%itf_index_info)
      if (associated(bcontr%idx_op1)) deallocate(bcontr%idx_op1)
      if (associated(bcontr%idx_op2)) deallocate(bcontr%idx_op2)
      if (associated(bcontr%idx_op1op2)) deallocate(bcontr%idx_op1op2)

      bcontr%occ_res => null()
      bcontr%occ_op1 => null()
      bcontr%occ_op2 => null()
      bcontr%rst_res => null()
      bcontr%rst_op1 => null()
      bcontr%rst_op2 => null()
      bcontr%occ_ex1 => null()
      bcontr%occ_ex2 => null()
      bcontr%occ_cnt => null()
      bcontr%rst_cnt => null()
      bcontr%rst_ex1 => null()
      bcontr%rst_ex2 => null()
      bcontr%merge_op1 => null()
      bcontr%merge_op2 => null()
      bcontr%merge_op1op2 => null()
      bcontr%merge_op2op1 => null()
      bcontr%itf_index_info => null()
      bcontr%idx_op1 => null()
      bcontr%idx_op2 => null()
      bcontr%idx_op1op2 => null()

      return
      end
