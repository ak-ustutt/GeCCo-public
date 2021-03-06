*----------------------------------------------------------------------*
      integer function form_count(form_head)
*----------------------------------------------------------------------*
*     count and print number of formula items on linked list
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      type(formula_item), intent(in), target ::
     &     form_head

      integer ::
     &     idx

      type(formula_item), pointer ::
     &     form_ptr

      idx = 0
      form_ptr => form_head
      do
        idx = idx + 1
        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do
      form_count = idx

      return
      end

