*----------------------------------------------------------------------*
      subroutine write_form_list(ffform,form_head,title)
*----------------------------------------------------------------------*
*     write formula on linked list to file ffform
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_formula_item.h'

      character, intent(in) ::
     &     title*(*)
      type(filinf), intent(inout) ::
     &     ffform
      type(formula_item), intent(in), target ::
     &     form_head

      type(formula_item), pointer ::
     &     form_ptr
      logical ::
     &     closeit
      integer ::
     &     nterms

      if (ffform%unit.le.0) then
        call file_open(ffform)
        closeit = .true.
      else
        closeit = .false.
      end if
      write(ffform%unit) len_trim(title),title
      form_ptr => form_head
      nterms = 0
      do
        nterms = nterms+1
        call wr_formula(ffform,form_ptr)

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do
      write(lulog,*) 'wrote ',nterms,' entries'

      if (closeit)
     &     call file_close_keep(ffform)

      return
      end

