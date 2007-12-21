*----------------------------------------------------------------------*
      subroutine update_mel_arr(op_info)
*----------------------------------------------------------------------*
*     setup or update direct acces pointer array for linked list
*     here: ME-lists
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

