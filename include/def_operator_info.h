      type operator_info
        integer ::
     &       nops, nmels, id_cnt
	integer, pointer ::
     &       op2list(:)
        type(operator_list), pointer ::
     &       op_list
        type(operator_array), pointer ::
     &       op_arr(:)
        type(me_list_list), pointer ::
     &       mel_list
        type(me_list_array), pointer ::
     &       mel_arr(:)
      end type operator_info
