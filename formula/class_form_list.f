*----------------------------------------------------------------------*
      subroutine class_form_list(lulog,form_head,op_info,mode)
*----------------------------------------------------------------------*
*     print formula on linked list to unit lulog
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     lulog, mode
      type(formula_item), intent(in), target ::
     &     form_head
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idx

      type(formula_item), pointer ::
     &     form_ptr

      idx = 0
      form_ptr => form_head
      do
        select case(form_ptr%command)
        case(command_end_of_formula)
          write(lulog,*) '[END]'
        case(command_set_target_init)
c          write(lulog,*) '[INIT TARGET]',form_ptr%target
          write(lulog,*) '   term #    class'
        case(command_set_target_update)
          write(lulog,*) '[SET TARGET]',form_ptr%target
        case(command_new_intermediate)
          write(lulog,*) '[NEW INTERMEDIATE]',form_ptr%target
        case(command_del_intermediate)
          write(lulog,*) '[DELETE INTERMEDIATE]',form_ptr%target
        case(command_add_contribution)
          idx = idx+1
c          write(lulog,*) '[ADD]',form_ptr%target,'( term #',idx,')'
          select case(mode)
          case(1)
            call class_contr1(lulog,form_ptr%contr,op_info,idx)
          case default
            write(lulog,*) 'unknown mode ',mode
          end select  
        case(command_symmetrise)
          write(lulog,*) '[SYMMETRISE TARGET]',form_ptr%target
        case default
          write(lulog,*) 'unknown command ',form_ptr%command,
     &                   form_ptr%target
        end select

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      return
      end

