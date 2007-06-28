*----------------------------------------------------------------------*
      subroutine new_formula_item(form_pnt,command,target)
*----------------------------------------------------------------------*
*     allocate new entry for linked list
*     on entry, form_pnt should point to last entry (containing [END])
*     on exit, form_pnt is ready to take the information coming
*       along with command "command" (i.e. contr or interm is alloc.d)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'
      
      type(formula_item), intent(inout), target ::
     &     form_pnt
      integer, intent(in) ::
     &     command, target

      select case(command)
      case(command_set_target_init,
     &     command_set_target_update,
     &     command_del_intermediate)
      ! nothing to do
      case(command_add_contribution)
        allocate(form_pnt%contr)
        call init_contr(form_pnt%contr)
      case(command_new_intermediate)
        allocate(form_pnt%interm)
      case default
        call quit(1,'new_formula_item','unknown mode')
      end select

      form_pnt%command = command
      form_pnt%target = target
      
      allocate(form_pnt%next)
      form_pnt%next%prev => form_pnt
      form_pnt%next%next => null()
      form_pnt%next%contr => null()
      form_pnt%next%interm => null()
      form_pnt%next%command = command_end_of_formula
      
      return
      end
