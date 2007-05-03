*----------------------------------------------------------------------*
      subroutine write_form_list(ffform,form_head)
*----------------------------------------------------------------------*
*     write formula on linked list to file ffform
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_formula.h'

      type(filinf), intent(inout) ::
     &     ffform
      type(formula), intent(in), target ::
     &     form_head

      type(formula), pointer ::
     &     form_ptr
      logical ::
     &     closeit

      if (ffform%unit.le.0) then
        call file_open(ffform)
        closeit = .true.
      else
        closeit = .false.
      end if
      write(ffform%unit) 10,'test title'
      form_ptr => form_head
      do
        call wr_formula(ffform,form_ptr)

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      if (closeit)
     &     call file_close_keep(ffform)

      return
      end

