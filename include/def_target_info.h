      type target_info

      integer ::
     &     ntargets, nactions
      integer(8), pointer ::
     &     last_mod(:)
      type(target_list), pointer ::
     &     list
      type(target_array), pointer ::
     &     array(:)

      type(action_list), pointer ::
     &     act_list
      type(action_array), pointer ::
     &     act_array(:)

      end type target_info
