*----------------------------------------------------------------------*
      subroutine set_hhat2(formula_hhat,op_info,
     &     idxhhat,idxham,idxtop)
*----------------------------------------------------------------------*
*     generate the formula for Hhat = e^{-T1}He^{T1}
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'par_formnames_gen.h'

      type(formula), intent(inout), target ::
     &     formula_hhat

      integer, intent(in) ::
     &     idxhhat,idxham,idxtop

      type(operator_info), intent(in) ::
     &     op_info

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_hhat
      type(formula_item), pointer ::
     &     form_pnt

      integer ::
     &     nterms 

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(luout,*) '======================'
        write(luout,*) ' output from set_hhat'
        write(luout,*) '======================'
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! initialize formula
      call init_formula(form_hhat)
      form_pnt => form_hhat
      ! put [INIT] at the beginning
      call new_formula_item(form_pnt,command_set_target_init,idxhhat)
      form_pnt => form_pnt%next

      ! expand e^{-T1} H e^T1 
      call expand_op_bch(form_pnt,4,idxhhat,
     &     1d0,-1,idxham,1d0,idxtop,1,1,op_info)

      ! reorder
      call reorder_formula(form_hhat,op_info)
c dbg
      print *,'generated (1):'
      call print_form_list(luout,form_hhat,op_info)
c dbg      
      ! write to disc
      formula_hhat%label = label_cchhat
      formula_hhat%comment = title_cchhat
      write(name,'(a,".fml")') label_cchhat
      call file_init(formula_hhat%fhand,name,ftyp_sq_unf,0)
      call write_form_list(formula_hhat%fhand,form_hhat,title_cchhat)

      call dealloc_formula_list(form_hhat)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'CC Hhat',cpu-cpu0,sys-sys0,wall-wall0)

      end
