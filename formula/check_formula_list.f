*----------------------------------------------------------------------*
      subroutine check_formula_list(formula_list,op_info)
*----------------------------------------------------------------------*
*     check formula list, in pt. whether the chain is linked
*     consistently, and whether all necessary subarrays are allocated
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'

      integer, parameter :: maxlen = 10 000 000

      type(formula_item), intent(in), target ::
     &     formula_list
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     current

      integer ::
     &     idx, isave, warning, error
      
      current => formula_list

      warning = 0
      error = 0
      ! test initial node
      if (current%command.ne.command_set_target_init) then
        write(luout,*) 'formula does not start with [INIT]'
        warning = warning+1
      end if
      if (associated(current%prev)) then
        write(luout,*) 'there is a link to a previous item'// 
     &                  ' (not at beginning of formula?)'
        warning = warning+1
      end if

      isave = 0
      do 
        select case(current%command)
        case(command_set_target_init,
     &     command_set_target_update,
     &     command_symmetrise,
     &     command_internal,command_end_of_formula)
        ! nothing to do
        case(command_add_contribution)
        if (.not.associated(current%contr)) then
          error = error + 1
          write(luout,*) '[CONTR] without allocated data structure' 
        end if
        case(command_del_intermediate)
        if (.not.associated(current%label)) then
          error = error + 1
          write(luout,*) '[DEL] without allocated label'
        end if
        case(command_new_intermediate)
        if (.not.associated(current%interm)) then
          error = error + 1
          write(luout,*) 
     &       '[NEW INTERMEDIATE] without allocated data structure'
        end if
        if (.not.associated(current%parent1)) then
          error = error + 1
          write(luout,*) '[NEW INTERMEDIATE] without parent1 label'
        end if
        if (.not.associated(current%parent2)) then
          error = error + 1
          write(luout,*) '[NEW INTERMEDIATE] without parent2 label'
        end if
        case(command_reorder)
        if (.not.associated(current%reo)) then
          error = error + 1
          write(luout,*) 
     &       '[REORDER] without allocated data structure'
        end if
        case(command_add_bc_reo,command_bc_reo,command_add_reo)
        if (.not.associated(current%reo)) then
          error = error + 1
          write(luout,*) 
     &       '[CONTR][REORDER] without allocated data structure "reo"'
        end if
        if (.not.associated(current%bcontr)) then
          error = error + 1
          write(luout,*) 
     &      '[CONTR][REORDER] without allocated data structure "bcontr"'
        end if
        case(command_add_intm,command_cp_intm,command_add_bc,command_bc)
        if (.not.associated(current%bcontr)) then
          error = error + 1
          write(luout,*) 
     &      '[CONTR] without allocated data structure "bcontr"'
        end if
        case default
        error = error + 1
        write (luout,*) 'unknown command ', current%command
        end select

        if (current%command.eq.command_end_of_formula) then
          if (associated(current%next)) then
            warning = warning+1
            write(luout,*) 
     &          'there is still an allocated item after [END]'
          end if
          exit
        end if

        if (.not.associated(current%next)) then
          error = error + 1
          write(luout,*) 'chain ends before [END]'
          exit
        end if
        if (.not.associated(current%next%prev,current)) then
          error = error + 1
          write(luout,*) 'forward-backward inconsistency detected'
          idx = -1
          call print_form_item(luout,idx,current,op_info)
          write(luout,*) '... but current%next%prev is ...'
          idx = -1
          call print_form_item(luout,idx,current%next%prev,op_info)
        end if
        current => current%next
        if (.not.associated(current%prev%next,current)) then
          error = error + 1
          write(luout,*) 'backward-forward inconsistency detected'
          idx = -1
          call print_form_item(luout,idx,current,op_info)
          write(luout,*) '... but current%prev%next is ...'
          idx = -1
          call print_form_item(luout,idx,current%prev%next,op_info)
        end if

        isave = isave+1
        if (isave.gt.maxlen) then
          error = error + 1
          write(luout,*) 'Infinite loop? Or formula longer than ',maxlen
          exit
        end if
      end do
     
      if (warning.gt.0) then
        call warn('check_formula_list','unclean formula')
      end if
      if (error.gt.0) then
        call quit(0,'check_formula_list','buggy formula')
      end if

      return
      end
