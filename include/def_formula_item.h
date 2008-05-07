      type formula_item
        ! command: see below
        ! target: index of operator
        integer ::
     &       command, target
        ! definition of temporary intermediate
        type(operator), pointer ::
     &       interm
        ! definition of contraction
        type(contraction), pointer ::
     &       contr
        type(formula_item), pointer ::
     &       prev, next
      end type formula_item

      integer, parameter ::
     &     command_end_of_formula = -1,
     &     command_set_target_init = 0,   
     &     command_set_target_update = 1,
     &     command_new_intermediate = 2,
     &     command_del_intermediate = 3,
     &     command_add_contribution = 4,
     &     command_symmetrise = 5,
     &     command_internal = 6


      ! 0: target operator is op(target), initialize to 0
      ! 1: dto. but update only
      ! 2: define new intermediate, target is new index
      !    sets target to that operator
      ! 3: delete intermediate associated with op(target)
      ! 4: add the contribution described by contr
