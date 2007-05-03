*----------------------------------------------------------------------*
      subroutine factorize(form_head,op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     read formula from file ffform to linked list 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_contraction.h'
      include 'def_formula.h'

      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(formula), intent(in), target ::
     &     form_head

      type(formula), pointer ::
     &     form_ptr
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      call atim(cpu0,sys0,wall0)

      form_ptr => form_head
      do while(form_ptr%command.ne.command_end_of_formula)
        if (form_ptr%command.eq.command_add_contribution) then
          call form_fact(form_ptr%contr,
     &         op_info,str_info,orb_info)
        end if
        if (.not.associated(form_ptr%next))
     &       call quit(1,'form_opt','buggy formula list')
        form_ptr => form_ptr%next
      end do

      call atim(cpu,sys,wall)
      if (iprlvl.ge.5) 
     &     call prtim(luout,'factorization',
     &     cpu-cpu0,sys-sys0,wall-wall0)


      return
      end

