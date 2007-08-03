*----------------------------------------------------------------------*
      subroutine set_r12_lagrangian(form_cclag,op_info,orb_info,
     &     idxham,idxtbar,idxrba,idxcba,idxtop,idxr12,idxc12,idxecc)
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
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'par_formnames_gen.h'
      include 'par_opnames_gen.h'
      include 'explicit.h'

      type(formula), intent(inout), target ::
     &     form_cclag

      integer, intent(in) ::
     &     idxham,idxtbar,idxtop,idxecc,idxrba,idxcba,idxr12,idxc12

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf) ::
     &     orb_info

      ! local variables

      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_lag, form_t_cr, form_tbar_cbarr
      type(formula_item), pointer ::
     &     form_pnt, fl_t_cr_pnt

      integer ::
     &     nterms, idx_sop, idx_sbar, idx_r12, idx_top, ndef
      integer ::
     &     occ_def(ngastp,2,2)

      type(operator), pointer::
     &     sop_pnt, sbar_pnt

      integer, external::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(luout,*) '==============================='
        write(luout,*) ' output from set_r12_lagrangian'
        write(luout,*) '==============================='
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! Definition of the S=T+CR operator.
      call add_operator(op_sop,op_info)
      idx_sop = idx_oplist2(op_sop,op_info)
      sop_pnt => op_info%op_arr(idx_sop)%op

      ! set CR part:
      idx_r12 = idx_oplist2(op_r12,op_info)
      if (trir12.eq.0) then
        call clone_operator(sop_pnt,op_info%op_arr(idx_r12)%op,orb_info)
      else
        ndef = 1
        occ_def(1:ngastp,1,1) = (/0,1,0,2/)
        occ_def(1:ngastp,2,1) = (/3,0,0,0/)
        if (ansatze.gt.1) then
          ndef = 2
          occ_def(1:ngastp,1,2) = (/0,2,0,1/)
          occ_def(1:ngastp,2,2) = (/3,0,0,0/)
        end if
        call set_uop(sop_pnt,op_sop,.false.,0,0,1,1,0,
     &       occ_def,ndef,orb_info)
        call join_operator(sop_pnt,op_info%op_arr(idx_r12)%op,orb_info)
      end if

      ! join with T
      idx_top = idx_oplist2(op_top,op_info)
      call join_operator(sop_pnt,op_info%op_arr(idx_top)%op,orb_info)

      ! define Sbar for the projection
      call add_operator(op_sba,op_info)
      idx_sbar = idx_oplist2(op_sba,op_info)
      sbar_pnt => op_info%op_arr(idx_sbar)%op
      call clone_operator(sbar_pnt,sop_pnt,orb_info)
      sbar_pnt%dagger = .true.

      ! When doing R12 calculation, must first combine the C and R 
      ! operators (S=T+CR).
      call init_formula(form_t_cr)
      fl_t_cr_pnt => form_t_cr
      call new_formula_item(fl_t_cr_pnt,command_set_target_init,idx_sop)
      fl_t_cr_pnt => fl_t_cr_pnt%next
      call expand_op_product(fl_t_cr_pnt,idx_sop,
     &     1d0,1,idxtop,-1,-1,
     &     0,0,op_info)
      do while(associated(fl_t_cr_pnt%next))
        fl_t_cr_pnt => fl_t_cr_pnt%next
      end do
      call expand_op_product(fl_t_cr_pnt,idx_sop,
     &     1d0,2,(/idxc12,idxr12/),-1,-1,
     &     (/1,2/),1,op_info)

      call write_title(luout,wst_title,'T + CR')
      call print_form_list(luout,form_t_cr,op_info)

      ! Must also form SBAR.
      call init_formula(form_tbar_cbarr)
      fl_t_cr_pnt => form_tbar_cbarr
      call new_formula_item(fl_t_cr_pnt,
     &                      command_set_target_init,idx_sbar)
      fl_t_cr_pnt => fl_t_cr_pnt%next
      call expand_op_product(fl_t_cr_pnt,idx_sbar,
     &     1d0,1,idxtbar,-1,-1,
     &     0,0,op_info)
      do while(associated(fl_t_cr_pnt%next))
        fl_t_cr_pnt => fl_t_cr_pnt%next
      end do
      call expand_op_product(fl_t_cr_pnt,idx_sbar,
     &     1d0,2,(/idxrba,idxcba/),-1,-1,
     &     (/1,2/),1,op_info)

      call write_title(luout,wst_title,'TBAR + R CBAR')
      call print_form_list(luout,form_tbar_cbarr,op_info)

      ! and now: the actual formula
      ! initialize formula
      call init_formula(form_lag)
      form_pnt => form_lag
      ! put [INIT] at the beginning
      call new_formula_item(form_pnt,command_set_target_init,idxecc)
      form_pnt => form_pnt%next

      ! expand <0|(1+Sbar) e^{-S} H e^S|0> =
      ! <0| e^{-S} H e^S|0> +
      call expand_op_bch(form_pnt,2,idxecc,
     &     1d0,-1,idxham,1d0,idx_sop,1,-1,op_info)

      ! advance pointer
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      end do
      ! <0|Sbar e^{-S} H e^S|0>
      call expand_op_bch(form_pnt,4,idxecc,
     &     1d0,idx_sbar,idxham,1d0,idx_sop,1,-1,op_info)
      ! insert here procedure to produce approx. expansions      

      call write_title(luout,wst_title,'raw formula')
      call print_form_list(luout,form_lag,op_info)

      ! replace S by T+CR
      call expand_subexpr(form_lag,form_t_cr,op_info)
      call sum_terms(form_lag,op_info)

      call write_title(luout,wst_title,'after replacing S')
      call print_form_list(luout,form_lag,op_info)

      ! replace Sbar by Tbar + R^t CBAR
      call expand_subexpr(form_lag,form_tbar_cbarr,op_info)

      call write_title(luout,wst_title,'Final formula')
      call print_form_list(luout,form_lag,op_info)

      ! to come
      ! define X, B, ... in terms of R12 
      ! use factor_out_subexpr to express Lagrangian through X, B
      ! i.e. eliminate all formal operators

      ! define X, B, ... in terms of actually available integrals
      ! 

c      ! post_processing and term counting:
c      call cc_form_post(form_lag,nterms,idxsba,idxham,idxsop,op_info)

      ! assign canonical name and comment
      form_cclag%label = label_cclg0
      form_cclag%comment = title_cclg0
      ! write to disc
      write(name,'(a,".fml")') label_cclg0
      call file_init(form_cclag%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_cclag%fhand,form_lag,title_cclg0)
      call dealloc_formula_list(form_lag)

      call mem_map(.true.)
      ! remove the formal operators
      call del_operator(idx_sbar,op_info)
      call del_operator(idx_sop,op_info)

      call atim_csw(cpu,sys,wall)
      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'CC Lagrangian',cpu-cpu0,sys-sys0,wall-wall0)

      call quit(1,'test','exit')

      return
      end
