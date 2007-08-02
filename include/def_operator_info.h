      type operator_info
        integer ::
     &       nops, id_cnt
	integer, pointer ::
     &       idx2id(:)
        type(operator_list), pointer ::
     &       op_list
        type(operator_array), pointer ::
     &       op_arr(:)
        type(file_list), pointer ::  ! obsolete
     &       opfil_list              ! obsolete
        type(file_array), pointer ::
     &       opfil_arr(:)
      end type operator_info
