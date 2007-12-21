      type target_info

      integer ::
     &     ntargets
      integer(8), pointer ::
     &     last_mod(:)
      type(target_list), pointer ::
     &     list
      type(target_array), pointer ::
     &     array(:)

      end type target_info
