      type formula_item
        ! command: see below
        ! target: index of operator
        integer ::
     &       command, target
        ! definition of temporary intermediate
        type(operator), pointer ::
     &       interm
	character(len=len_opname), pointer ::
     &       parent1, parent2
        logical ::
     &       tra, tra1, tra2 ! interm./parents transposed?
        ! definition of contraction
        type(contraction), pointer ::
     &       contr
        ! pre-processed binary contraction
	type(binary_contr), pointer ::
     &       bcontr 
        ! pre-processed reordering info
	type(reorder), pointer ::
     &       reo
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
     &     command_reorder = 6,
     &     command_add_intm = 7,
     &     command_add_bc = 8,
     &     command_add_bc_reo = 9,
     &     command_bc = 10,
     &     command_bc_reo = 11,
     &     command_internal = 12,
     &     command_add_reo = 13


      ! 0: target operator is op(target), initialize to 0
      ! 1: dto. but update only
      ! 2: define new intermediate, target is new index
      !    sets target to that operator
      ! 3: delete intermediate associated with op(target)
      ! 4: add the contribution described by contr
