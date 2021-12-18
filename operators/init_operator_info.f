      subroutine init_operator_info(op_info)

      implicit none

      include 'mdef_operator_info.h'

      type(operator_info), intent(inout) ::
     &     op_info

      ! no operators yet
      op_info%nops = 0
      ! no ME-lists yet
      op_info%nmels = 0
      ! initialize operator list
      allocate(op_info%op_list)
      nullify(op_info%op_list%op)
      nullify(op_info%op_list%prev)
      nullify(op_info%op_list%next)
      ! initialize pointer array
      nullify(op_info%op_arr)
      ! initialize ME-list
      allocate(op_info%mel_list)
      nullify(op_info%mel_list%mel)
      nullify(op_info%mel_list%prev)
      nullify(op_info%mel_list%next)
      ! initialize pointer array
      nullify(op_info%mel_arr)
      ! initialize op -> list lookup table
      nullify(op_info%op2list)

      return
      end
