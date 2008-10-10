*----------------------------------------------------------------------*
*     Define the properties of the operator to search for when trying 
*     to delete formula nodes.
*----------------------------------------------------------------------*
      type del_cond
        character*64 ::
     &     op_name
        logical ::
     &     transposed
        integer ::
     &       num_op_restr(2),
     &       part_num_restr(2)
      end type del_cond
*----------------------------------------------------------------------*
*     Define the array for the deletion conditions.
*     The members of the array are linked by the logical AND operation.
*----------------------------------------------------------------------*
      type del_cond_array
        integer ::
     &       and_dim
        type(del_cond), pointer ::
     &     del_cond_arr(:)
      end type del_cond_array
*----------------------------------------------------------------------*
*     Define the list constaining the above arrays.
*     The members of the list are linked by the logical OR operation.
*----------------------------------------------------------------------*
      type del_cond_list
        integer ::
     &     or_dim
        type(del_cond_array), pointer ::
     &     del_cond_item(:)
      end type del_cond_list
