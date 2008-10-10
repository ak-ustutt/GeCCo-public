*----------------------------------------------------------------------*
      subroutine tex_form_list(lutex,form_head,op_info)
*----------------------------------------------------------------------*
*     print formula on linked list in TeX style to unit lutex
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, parameter ::
     &     max_count = 35

      integer, intent(in) ::
     &     lutex
      type(formula_item), intent(in), target ::
     &     form_head
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     form_ptr
      integer ::
     &     iblk_cur
      logical ::
     &     first, newline
      integer ::
     &     length_count, length_term

      length_count = 0
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
          if (length_count.gt.max_count) then
            newline = .true.
            length_count = 0
          end if
          call tex_contr(lutex,length_term,
     &                   first,newline,form_ptr%contr,op_info)
          length_count = length_count + length_term
          newline = .false.
          iblk_cur = form_ptr%contr%iblk_res
          first = .false.
        end select

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      return
      end

