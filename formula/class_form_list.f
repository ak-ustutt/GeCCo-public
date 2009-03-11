*----------------------------------------------------------------------*
      subroutine class_form_list(luout,form_head,op_info,mode)
*----------------------------------------------------------------------*
*     print formula on linked list to unit luout
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     luout, mode
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
          write(luout,*) '[END]'
        case(command_set_target_init)
c          write(luout,*) '[INIT TARGET]',form_ptr%target
          write(luout,*) '   term #    class'
        case(command_set_target_update)
          write(luout,*) '[SET TARGET]',form_ptr%target
        case(command_new_intermediate)
          write(luout,*) '[NEW INTERMEDIATE]',form_ptr%target
        case(command_del_intermediate)
          write(luout,*) '[DELETE INTERMEDIATE]',form_ptr%target
        case(command_add_contribution)
          idx = idx+1
c          write(luout,*) '[ADD]',form_ptr%target,'( term #',idx,')'
          select case(mode)
          case(1)
            call class_contr1(luout,form_ptr%contr,op_info,idx)
          case default
            write(luout,*) 'unknown mode ',mode
          end select  
        case(command_symmetrise)
          write(luout,*) '[SYMMETRISE TARGET]',form_ptr%target
        case default
          write(luout,*) 'unknown command ',form_ptr%command,
     &                   form_ptr%target
        end select

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      return
      end

