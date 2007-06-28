*----------------------------------------------------------------------*
      subroutine tex_form_list(luout,form_head,op_info)
*----------------------------------------------------------------------*
*     print formula on linked list in TeX style to unit luout
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     luout
      type(formula_item), intent(in), target ::
     &     form_head
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     form_ptr
      integer ::
     &     iblk_cur
      logical ::
     &     first

      form_ptr => form_head
      do
        select case(form_ptr%command)
        case(command_end_of_formula)
          exit
        case(command_set_target_init)
          first = .true.
        case(command_set_target_update)
          first = .true.
        case(command_new_intermediate)
          first = .true.
        case(command_del_intermediate)
        case(command_add_contribution)
          if (.not.first) then
            if (iblk_cur.ne.form_ptr%contr%iblk_res) first = .true.
          end if
          call tex_contr(luout,first,form_ptr%contr,op_info)
          iblk_cur = form_ptr%contr%iblk_res
          first = .false.
        end select

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      return
      end

