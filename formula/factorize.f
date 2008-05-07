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
      include 'def_formula_item.h'

      integer, parameter ::
     &     ntest = 00

      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(formula_item), intent(in), target ::
     &     form_head

      type(formula_item), pointer ::
     &     form_ptr
      integer ::
     &     iterm,
     &     iscale_stat(ngastp,2)
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      call atim_csw(cpu0,sys0,wall0)

      call write_title(luout,wst_subsection,
     &     'Formula factorization')

      iscale_stat = 0

      form_ptr => form_head
      iterm = 0
      do while(form_ptr%command.ne.command_end_of_formula)
        if (form_ptr%command.eq.command_add_contribution) then
          iterm = iterm + 1
          if (iprlvl.ge.10)
     &         write(luout,*) 'factorizing term # ',iterm
          call form_fact2(form_ptr%contr,
     &         op_info,str_info,orb_info,iscale_stat)
        end if
        if (.not.associated(form_ptr%next))
     &       call quit(1,'form_opt','buggy formula list')
        form_ptr => form_ptr%next
      end do

      call write_title(luout,wst_subsection,
     &     'Summary')
      
      write(luout,'(x,a,"H^",i2," P^",i2," V^",i2," X^",i2)')
     &       'Most expensive contraction:  ',iscale_stat(1:4,1)
      write(luout,'(x,a,"H^",i2," P^",i2," V^",i2," X^",i2)')
     &       'Largest intermediate      :  ',iscale_stat(1:4,2)

      call atim_csw(cpu,sys,wall)
      if (iprlvl.ge.5) 
     &     call prtim(luout,'factorization',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

