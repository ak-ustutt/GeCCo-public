*----------------------------------------------------------------------*
!>     setup or update direct acces pointer array for linked list
!>     here: ME-lists
!>
!>    destroy array and recreate from linked list
!>    @param op_info info struct of the operator whose array is updated.
*----------------------------------------------------------------------*
      subroutine update_mel_arr(op_info)
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_operator_info.h'

      type(operator_info) ::
     &     op_info
      integer ::
     &     imel

      if (associated(op_info%mel_arr)) deallocate(op_info%mel_arr)

      if (op_info%nmels.gt.0) then
        allocate(op_info%mel_arr(op_info%nmels))
        call mel_list2arr(op_info%mel_list,
     &                    op_info%mel_arr,op_info%nmels)
      end if

      return
      end

