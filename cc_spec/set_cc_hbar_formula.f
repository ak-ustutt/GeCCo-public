*----------------------------------------------------------------------*
      subroutine set_cc_hbar_formula(formula_hbar,
     &     title,name_hbar,name_h,name_t,
     &     op_info)
*----------------------------------------------------------------------*
*     generate the formula for Hbar = e^{-T}He^{T}
*     restricted to the shape given for the Hbar operator
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
c      include 'par_formnames_gen.h'

      type(formula), intent(inout), target ::
     &     formula_hbar

      character*(*), intent(in) ::
     &     name_hbar, name_t, name_h, title

      type(operator_info), intent(in) ::
     &     op_info

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_hbar
      type(formula_item), pointer ::
     &     form_pnt

      integer ::
     &     nterms,
     &     idxhbar,idxham,idxtop

      integer, external ::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(lulog,*) '======================'
        write(lulog,*) ' output from set_hbar'
        write(lulog,*) '======================'
      end if

      call atim_csw(cpu0,sys0,wall0)

      idxhbar = idx_oplist2(name_hbar,op_info)
      idxham  = idx_oplist2(name_h,op_info)
      idxtop  = idx_oplist2(name_t,op_info)
      if (idxham.lt.0.or.idxhbar.lt.0.or.idxtop.lt.0)
     &     call quit(1,'set_cc_hbar_formula',
     &     'required operators are not yet defined')

      ! initialize formula
      call init_formula(form_hbar)
      form_pnt => form_hbar
      ! put [INIT] at the beginning
      call new_formula_item(form_pnt,command_set_target_init,idxhbar)
      form_pnt => form_pnt%next

      ! expand e^{-T} H e^T 
      call expand_op_bch(form_pnt,4,idxhbar,
     &     1d0,-1,idxham,1d0,idxtop,-1,-1,op_info)

      ! reorder
      call reorder_formula(form_hbar,op_info)

      ! write to disc
      formula_hbar%comment = title
      write(name,'(a,".fml")') trim(formula_hbar%label)
      call file_init(formula_hbar%fhand,name,ftyp_sq_unf,0)
      call write_form_list(formula_hbar%fhand,form_hbar,title)

      if (ntest.ge.100) then
        write(lulog,*) 'generated:'
        call print_form_list(lulog,form_hbar,op_info)
      end if

      call dealloc_formula_list(form_hbar)

      call atim_csw(cpu,sys,wall)
      call prtim(lulog,'CC Hbar',cpu-cpu0,sys-sys0,wall-wall0)

      end
