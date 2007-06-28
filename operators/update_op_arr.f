*----------------------------------------------------------------------*
      subroutine update_op_arr(op_info)
*----------------------------------------------------------------------*
*     setup or update direct acces pointer array for linked list
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_operator_info.h'

      type(operator_info) ::
     &     op_info
      integer ::
     &     iop

      if (associated(op_info%op_arr)) deallocate(op_info%op_arr)

      if (op_info%nops.gt.0) then
        allocate(op_info%op_arr(op_info%nops))
        call op_list2arr2(op_info%op_list,op_info%op_arr,op_info%nops)
      end if

      ! update index -> ID lookup table
      if (associated(op_info%idx2id)) deallocate(op_info%idx2id)
      allocate(op_info%idx2id(op_info%nops))
      do iop = 1, op_info%nops
        op_info%idx2id(iop) = op_info%op_arr(iop)%op%id
      end do


      return
      end

