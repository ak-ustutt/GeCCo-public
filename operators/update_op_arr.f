*----------------------------------------------------------------------*
      subroutine update_op_arr(op_info)
*----------------------------------------------------------------------*
*     setup or update direct access pointer array for linked list
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_operator_info.h'

      type(operator_info) ::
     &     op_info
      integer ::
     &     iop, idx
      
      integer, external ::
     &     idx_mel_list

      if (associated(op_info%op_arr)) deallocate(op_info%op_arr)

      if (op_info%nops.gt.0) then
        allocate(op_info%op_arr(op_info%nops))
        call op_list2arr2(op_info%op_list,op_info%op_arr,op_info%nops)
      end if

      ! update op-index -> list-index lookup table
      if (associated(op_info%op2list)) deallocate(op_info%op2list)
      allocate(op_info%op2list(op_info%nops))
      do iop = 1, op_info%nops
        idx =
!     &     idx_mel_list(trim(op_info%op_arr(iop)%op%assoc_list),op_info)
     &     idx_mel_list(op_info%op_arr(iop)%op%assoc_list,op_info)
! note for the above: trim is used inside the routine anyway
        op_info%op2list(iop) = idx
      end do

      return
      end

