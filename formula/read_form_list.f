*----------------------------------------------------------------------*
      subroutine read_form_list(ffform,form_head,init)
*----------------------------------------------------------------------*
*     read formula from file ffform to linked list 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_formula_item.h'

      type(filinf), intent(inout) ::
     &     ffform
      type(formula_item), intent(in), target ::
     &     form_head
      logical, intent(in) ::
     &     init

      type(formula_item), pointer ::
     &     form_ptr
      logical ::
     &     closeit

      logical, external ::
     &     rd_formula
c dbg
      integer ::
     &     nterms
c dbg
      if (ffform%unit.le.0) then
        call file_open(ffform)
        closeit = .true.
      else
        closeit = .false.
      end if

      rewind ffform%unit
      ! preliminary: overread title record
      read(ffform%unit)

      form_ptr => form_head
      ! we need init_formula here??? -->yes!
      if (init) call init_formula(form_ptr)
      nterms = 1 ! there is at least the [END] entry (and we exit befor incr. nterms)
      do while(rd_formula(ffform,form_ptr))
        if (form_ptr%command.eq.command_end_of_formula) exit
        nterms = nterms+1
c        allocate(form_ptr%next)
c        form_ptr%next%prev => form_ptr
        form_ptr => form_ptr%next
c        nullify(form_ptr%next)
c        nullify(form_ptr%contr)
c        nullify(form_ptr%interm)
c        form_ptr%command = command_end_of_formula
      end do
      write(lulog,*) 'read ',nterms,' entries'

      if (closeit)
     &     call file_close_keep(ffform)

      return
      end

