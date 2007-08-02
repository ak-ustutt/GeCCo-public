*----------------------------------------------------------------------*
      subroutine set_r12_lagrangian(form_cclag,op_info,
     &     idxham,idxtbar,idxrba,idxcba,idxtop,idxr12,idxc12,idxecc,
     &     idxsop,idxsba)
*----------------------------------------------------------------------*
*
*     set up sequence of operators, integrals and contractions that
*     defines an R12-Lagrangian within the chosen operator space 
*
*     modified from set_cc_lagrangian2 by GWR June 2007
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
     &     idxham,idxtbar,idxtop,idxecc,idxrba,idxcba,idxr12,idxc12,
     &     idxsop,idxsba

      type(operator_info), intent(inout) ::
     &     op_info

      ! local variables

      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_lag, form_t_cr
      type(formula_item), pointer ::
     &     form_pnt, fl_t_cr_pnt

      integer ::
     &     nterms

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(luout,*) '==============================='
        write(luout,*) ' output from set_r12_lagrangian'
        write(luout,*) '==============================='
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! initialize formula
      call init_formula(form_lag)
      form_pnt => form_lag
      ! put [INIT] at the beginning
      call new_formula_item(form_pnt,command_set_target_init,idxecc)
      form_pnt => form_pnt%next

      ! When doing R12 calculation, must first combine the T and R 
      ! operators (S=T+CR).
      call init_formula(form_t_cr)
      fl_t_cr_pnt => form_t_cr
      call new_formula_item(fl_t_cr_pnt,command_set_target_init,idxsop)
      fl_t_cr_pnt => fl_t_cr_pnt%next
      call expand_op_product(fl_t_cr_pnt,idxsop,
     &     1d0,1,idxtop,-1,-1,
     &     0,0,op_info)
      do while(associated(fl_t_cr_pnt%next))
        fl_t_cr_pnt => fl_t_cr_pnt%next
      end do
      call expand_op_product(fl_t_cr_pnt,idxsop,
     &     1d0,2,(/idxc12,idxr12/),-1,-1,
     &     (/1,2/),1,op_info)

      call print_form_list(luout,form_t_cr,op_info)

      ! Must also form S+.
      fl_t_cr_pnt => fl_t_cr_pnt%next      
      call new_formula_item(fl_t_cr_pnt,command_set_target_init,idxsba)
      fl_t_cr_pnt => fl_t_cr_pnt%next
      call expand_op_product(fl_t_cr_pnt,idxsba,
     &     1d0,1,idxtbar,-1,-1,
     &     0,0,op_info)
      do while(associated(fl_t_cr_pnt%next))
        fl_t_cr_pnt => fl_t_cr_pnt%next
      end do
      call expand_op_product(fl_t_cr_pnt,idxsba,
     &     1d0,2,(/idxrba,idxcba/),-1,-1,
     &     (/1,2/),1,op_info)

      call print_form_list(luout,form_t_cr,op_info)

      ! expand <0|(1+Tbar) e^{-T} H e^T|0> =
      ! <0| e^{-T} H e^T|0> +
      call expand_op_bch(form_pnt,2,idxecc,
     &     1d0,-1,idxham,1d0,idxsop,1,2,op_info)

      ! advance pointer
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      end do
      ! <0|Tbar e^{-T} H e^T|0>
      call expand_op_bch(form_pnt,4,idxecc,
     &     1d0,idxsba,idxham,1d0,idxsop,1,-1,op_info)
      ! insert here procedure to produce approx. expansions      

      ! post_processing and term counting:
      call cc_form_post(form_lag,nterms,idxsba,idxham,idxsop,op_info)

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

      return
      end
