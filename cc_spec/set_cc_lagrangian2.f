*----------------------------------------------------------------------*
      subroutine set_cc_lagrangian2(form_cclag,op_info,
     &     idxham,idxtbar,idxtop,idxecc)
*----------------------------------------------------------------------*
*
*     set up sequence of operators, integrals and contractions that
*     defines a CC-Lagrangian within the chosen operator space 
*
*     new version using gen_contr()
*     written by andreas, june 2007
*
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
      include 'par_formnames_gen.h'

      type(formula), intent(inout), target ::
     &     form_cclag

      integer, intent(in) ::
     &     idxham,idxtbar,idxtop,idxecc

      type(operator_info), intent(in) ::
     &     op_info

      ! local variables

      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_lag
      type(formula_item), pointer ::
     &     form_pnt

      integer ::
     &     nterms 

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(luout,*) '==============================='
        write(luout,*) ' output from set_cc_lagrangian'
        write(luout,*) '==============================='
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! initialize formula
      call init_formula(form_lag)
      form_pnt => form_lag
      ! put [INIT] at the beginning
      call new_formula_item(form_pnt,command_set_target_init,idxecc)
      form_pnt => form_pnt%next

      ! expand <0|(1+Tbar) e^{-T} H e^T|0> =
      ! <0| e^{-T} H e^T|0> +
      call expand_op_bch(form_pnt,2,idxecc,
     &     1d0,-1,idxham,1d0,idxtop,1,2,op_info)

      ! advance pointer
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      end do
      ! <0|Tbar e^{-T} H e^T|0>
      call expand_op_bch(form_pnt,4,idxecc,
     &     1d0,idxtbar,idxham,1d0,idxtop,1,-1,op_info)

      ! insert here procedure to produce approx. expansions      

      ! post_processing and term counting:
      call cc_form_post(form_lag,nterms,idxtbar,idxham,idxtop,op_info)

      ! assign canonical name and comment
      form_cclag%label = label_cclg0
      form_cclag%comment = title_cclg0
      ! write to disc
      write(name,'(a,".fml")') label_cclg0
      call file_init(form_cclag%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_cclag%fhand,form_lag,title_cclg0)

      call dealloc_formula_list(form_lag)

      call atim_csw(cpu,sys,wall)
      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'CC Lagrangian',cpu-cpu0,sys-sys0,wall-wall0)

      end
