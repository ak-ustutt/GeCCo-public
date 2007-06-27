      subroutine init_operator_info(op_info)

      implicit none

      include 'mdef_operator_info.h'

      type(operator_info), intent(inout) ::
     &     op_info

      ! no operators yet
      op_info%nops = 0
      ! initialize operator list
      allocate(op_info%op_list)
      nullify(op_info%op_list%op)
      nullify(op_info%op_list%prev)
      nullify(op_info%op_list%next)
      ! initialize pointer array
      nullify(op_info%op_arr)
      ! initialize operator file list
      allocate(op_info%opfil_list)
      nullify(op_info%opfil_list%fhand)
      nullify(op_info%opfil_list%prev)
      nullify(op_info%opfil_list%next)
      ! initialize pointer array
      nullify(op_info%opfil_arr)
      ! initialize ID-counter
      op_info%id_cnt = 0
      nullify(op_info%idx2id)

      return
      end
