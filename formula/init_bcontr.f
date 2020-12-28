*----------------------------------------------------------------------*
      subroutine init_bcontr(bcontr)
*----------------------------------------------------------------------*
*     initialize all pointers
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(binary_contr), intent(inout) ::
     &     bcontr

      bcontr%n_operands = 0
      bcontr%n_cnt = 0
      bcontr%perm = .false.

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
      bcontr%index_info_available = .false.
      bcontr%idx_op1 => null()
      bcontr%idx_op2 => null()
      bcontr%idx_op1op2 => null()

      return
      end
