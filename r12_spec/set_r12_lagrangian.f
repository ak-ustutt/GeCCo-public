*----------------------------------------------------------------------*
      subroutine set_r12_lagrangian(form_cclag,op_info,orb_info,
     &     idxham,idxtbar,idxrba,idxcba,idxtop,idxr12,idxc12,idxlcc)
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
     &     ntest = 100

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
     &     idxham,idxtbar,idxtop,idxlcc,idxrba,idxcba,idxr12,idxc12

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf), intent(inout) ::
     &     orb_info

      ! local variables

      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_lag, form_t_cr, form_tbar_cbarr, form_v, form_gr,
     &     form_b, form_b_temp,form_rfr, form_h, form_v_bar,
     &     form_gr_bar, form_gr_temp_3, form_r
      type(formula_item), pointer ::
     &     form_pnt, fl_t_cr_pnt, form_gr_pnt, form_rfr_pnt, form_h_pnt,
     &     form_gr_bar_pnt, form_gr_temp_3_pnt, form_r_pnt

      integer ::
     &     nterms, idx_sop, idx_sbar, idx_r12, idx_top, ndef, idxrint,
     &     idxopv, idxopb, idx_ctemp, idx_ftemp, idxopvba, idx_v,
     &     idx_v1, idx_v2, idx_v3, idx_gtemp, idx_rtemp
      integer, allocatable ::
     &     occ_def(:,:,:)
      type(operator_array), pointer ::
     &     ops_array(:)

      type(operator), pointer::
     &     sop_pnt, sbar_pnt, ctemp_pnt, ftemp_pnt, v_pnt, v1_pnt,
     &     v2_pnt, v3_pnt, g_temp_pnt, r_temp_pnt

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
        allocate (occ_def(ngastp,2,2))
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
        deallocate(occ_def)
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
     &     0,0,.false.,op_info)
      do while(associated(fl_t_cr_pnt%next))
        fl_t_cr_pnt => fl_t_cr_pnt%next
      end do
      call expand_op_product(fl_t_cr_pnt,idx_sop,
     &     1d0,2,(/idxc12,idxr12/),-1,-1,
     &     (/1,2/),1,.false.,op_info)

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
     &     0,0,.false.,op_info)
      do while(associated(fl_t_cr_pnt%next))
        fl_t_cr_pnt => fl_t_cr_pnt%next
      end do
      call expand_op_product(fl_t_cr_pnt,idx_sbar,
     &     1d0,2,(/idxrba,idxcba/),-1,-1,
     &     (/1,2/),1,.false.,op_info)

      call write_title(luout,wst_title,'TBAR + R CBAR')
      call print_form_list(luout,form_tbar_cbarr,op_info)

      ! and now: the actual formula
      ! initialize formula
      call init_formula(form_lag)
      form_pnt => form_lag
      ! put [INIT] at the beginning
      call new_formula_item(form_pnt,command_set_target_init,idxlcc)
      form_pnt => form_pnt%next

      ! expand <0|(1+Sbar) e^{-S} H e^S|0> =
      ! <0| e^{-S} H e^S|0> +
      call expand_op_bch(form_pnt,2,idxlcc,
     &     1d0,-1,idxham,1d0,idx_sop,1,-1,op_info)

      ! advance pointer
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      end do
      ! <0|Sbar e^{-S} H e^S|0>
      call expand_op_bch(form_pnt,4,idxlcc,
     &     1d0,idx_sbar,idxham,1d0,idx_sop,1,-1,op_info)
      call sum_terms(form_lag,op_info)
      ! insert here procedure to produce approx. expansions      

      call write_title(luout,wst_title,'raw formula')
      call print_form_list(luout,form_lag,op_info)

      ! post_processing and term counting:
      call sum_terms(form_lag,op_info)
      call cc_form_post(form_lag,nterms,
     &     idx_sbar,idxham,idx_sop,op_info)

      ! replace S by T+CR
      call expand_subexpr(form_lag,form_t_cr,op_info)

      call write_title(luout,wst_title,'after replacing S')
      call print_form_list(luout,form_lag,op_info)

      ! replace Sbar by Tbar + R^t CBAR
      call expand_subexpr(form_lag,form_tbar_cbarr,op_info)

      call write_title(luout,wst_title,'Final formula')
      call print_form_list(luout,form_lag,op_info)

      ! Locate indices of the R12-method's operators.
      idxrint=idx_oplist2(op_rint,op_info)
      if(idxrint.le.0)
     &     call quit(1,'set_r12_lagrangian','No R12 integrals.')



      ! define X, B, ... in terms of R12 
      ! The form of these operators, in terms of the coefficients to
      ! which they must couple, has been setup in the operator 
      ! definition routine.
      if(ntest.ge.100)then
        write(luout,*)'Formation of R12 Intermediates.'
      endif  



      ! Define V= contraction between 2-electron part of H and R.
      idxopv=idx_oplist2(op_v_inter,op_info)
      call init_formula(form_gr)
      form_gr_pnt => form_gr
      call new_formula_item(form_gr_pnt,command_set_target_init,idxlcc)
      form_gr_pnt => form_gr_pnt%next

      call expand_op_product(form_gr_pnt,idxlcc,
     &     1d0,3,(/idxham,idxc12,idxr12/),(/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3/),1,.false.,op_info)

c dbg
c      call print_form_list(luout,form_gr,op_info)
c dbg

      call form_deriv3(form_v,form_gr,
     &     1,idxc12,0,idxopv,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'V Formula')
        call print_form_list(luout,form_v,op_info)
      endif  



      ! Define the adjoint of V.
      ! NB. This operator is not defined to be daggered when setup. This 
      ! needs to be remedied. Perhaps. GWR 14/08/2007
      idxopvba=idx_oplist2(op_vbar_inter,op_info)
      call init_formula(form_gr_bar)
      form_gr_bar_pnt => form_gr_bar
      call new_formula_item(form_gr_bar_pnt,command_set_target_init,
     &     idxlcc)
      form_gr_bar_pnt => form_gr_bar_pnt%next

      call expand_op_product(form_gr_bar_pnt,idxlcc,
     &     1d0,3,(/idxrba,idxcba,idxham/),(/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3/),1,.false.,op_info)
c dbg
c      call print_form_list(luout,form_gr_bar,op_info)
c dbg       

      call form_deriv3(form_v_bar,form_gr_bar,
     &     1,idxcba,0,idxopvba,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'V+ Formula')
        call print_form_list(luout,form_v_bar,op_info)
      endif  



      ! Define a temporary intermediate, cloned from C12. This used as
      ! a dummy to allow the derivative of B to be formed correctly.
      call add_operator(op_c_temp,op_info)
      idx_ctemp=idx_oplist2(op_c_temp,op_info)
      ctemp_pnt => op_info%op_arr(idx_ctemp)%op
      call clone_operator(ctemp_pnt,op_info%op_arr(idxc12)%op,orb_info)
      ! Also need to define a temporary dummy for the particle/external
      ! only Fock operators i.e. those which are required for the 
      ! intermediate, B.
      call add_operator(op_f_temp,op_info)
      idx_ftemp=idx_oplist2(op_f_temp,op_info)
      ftemp_pnt => op_info%op_arr(idx_ftemp)%op
      ndef = 4
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/0,1,0,0/)
      occ_def(1:ngastp,2,1) = (/0,1,0,0/)
      occ_def(1:ngastp,1,2) = (/0,1,0,0/)
      occ_def(1:ngastp,2,2) = (/0,0,0,1/)
      occ_def(1:ngastp,1,3) = (/0,0,0,1/)
      occ_def(1:ngastp,2,3) = (/0,1,0,0/)
      occ_def(1:ngastp,1,4) = (/0,0,0,1/)
      occ_def(1:ngastp,2,4) = (/0,0,0,1/)
      call set_uop(ftemp_pnt,op_f_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)
      ! Expand F(Temp) in terms of the actual Hamiltonian operator.
      call init_formula(form_h)
      form_h_pnt => form_h
      call new_formula_item(form_h_pnt,command_set_target_init,
     &     idx_ftemp)
      form_h_pnt => form_h_pnt%next
      call expand_op_product(form_h_pnt,idx_ftemp,
     &     1d0,1,idxham,-1,-1,
     &     0,0,.false.,op_info)

c dbg
c      call print_form_list(luout,form_rfr,op_info)
c dbg

      ! Define B= contraction between R+, F and R.
      idxopb=idx_oplist2(op_b_inter,op_info)
      call init_formula(form_rfr)
      form_rfr_pnt => form_rfr
      call new_formula_item(form_rfr_pnt,command_set_target_init,idxlcc)
      form_rfr_pnt => form_rfr_pnt%next
      call expand_op_product(form_rfr_pnt,idxlcc,
     &     1d0,4,(/idxrba,idx_ftemp,idx_ctemp,idxr12/),
     &     (/-1,-1,-1,-1/),(/-1,-1,-1,-1/),
     &     (/1,3,1,4,2,4,3,4/),4,.false.,op_info)

      ! Replace F(Temp) with H. Then form the derivative.
      call expand_subexpr(form_rfr,form_h,op_info)
      
      call form_deriv3(form_b,form_rfr,
     &     1,idx_ctemp,0,idxopb,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'B Formula')
        call print_form_list(luout,form_b,op_info)
      endif  

      ! use factor_out_subexpr to express Lagrangian through X, B
      ! i.e. eliminate all formal operators
      call factor_out_subexpr(form_lag,form_v,op_info)
      call factor_out_subexpr(form_lag,form_v_bar,op_info)
      call factor_out_subexpr(form_lag,form_b,op_info)

      call dealloc_formula_list(form_h)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'Factored R12 Lagrangian')
        call print_form_list(luout,form_lag,op_info)
      endif  



      ! Define V, B, etc. in terms of actually available integrals. Do 
      ! this by defining temporary intermediates of the same shape, then
      ! expand these in terms of the required integrals. Subsequently sum 
      ! the individual contributions.

      ! V.
      ! Set up temporary operator intermediates.
      call add_operator(op_v_temp,op_info)
      idx_v=idx_oplist2(op_v_temp,op_info)
      v_pnt => op_info%op_arr(idx_v)%op
      allocate(ops_array(2))
      ops_array(1)%op => op_info%op_arr(idxlcc)%op
      ops_array(2)%op => ctemp_pnt
      call set_gen_intermediate(v_pnt,op_v_temp,
     &     ops_array,2,orb_info)
      deallocate(ops_array)

      call add_operator(op_v_temp1,op_info)
      idx_v1=idx_oplist2(op_v_temp1,op_info)
      v1_pnt => op_info%op_arr(idx_v1)%op
      call clone_operator(v1_pnt,v_pnt,orb_info)

      call add_operator(op_v_temp2,op_info)
      idx_v2=idx_oplist2(op_v_temp2,op_info)
      v2_pnt => op_info%op_arr(idx_v2)%op
      call clone_operator(v2_pnt,v_pnt,orb_info)
      
      call add_operator(op_v_temp3,op_info)
      idx_v3=idx_oplist2(op_v_temp3,op_info)
      v3_pnt => op_info%op_arr(idx_v3)%op
      call clone_operator(v3_pnt,v_pnt,orb_info)

      ! Expand the temporary operators in terms of available integrals.

      ! Temporary operator 3.
      call add_operator(op_g_temp,op_info)
      idx_gtemp=idx_oplist2(op_g_temp,op_info)
      g_temp_pnt => op_info%op_arr(idx_gtemp)%op
      ndef=3
      allocate(occ_def(ngastp,2,3))
      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/2,0,0,0/)
      occ_def(1:ngastp,2,2) = (/1,1,0,0/)
      occ_def(1:ngastp,1,3) = (/2,0,0,0/)
      occ_def(1:ngastp,2,3) = (/0,2,0,0/)
      call set_uop(g_temp_pnt,op_g_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)
      ! Substitute with actual H integrals.
      call init_formula(form_h)
      form_h_pnt => form_h
      call new_formula_item(form_h_pnt,command_set_target_init,
     &     idx_gtemp)
      form_h_pnt => form_h_pnt%next
      call expand_op_product(form_h_pnt,idx_gtemp,
     &     1d0,1,idxham,-1,-1,
     &     0,0,.false.,op_info)

      call add_operator(op_r_temp,op_info)
      idx_rtemp=idx_oplist2(op_r_temp,op_info)
      r_temp_pnt => op_info%op_arr(idx_rtemp)%op
      ndef=3
      allocate(occ_def(ngastp,2,3))
      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/1,1,0,0/)
      occ_def(1:ngastp,2,2) = (/2,0,0,0/)
      occ_def(1:ngastp,1,3) = (/0,2,0,0/)
      occ_def(1:ngastp,2,3) = (/2,0,0,0/)
      call set_uop(r_temp_pnt,op_r_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)
      ! Substitute with actual R12 integrals.
      call init_formula(form_r)
      form_r_pnt => form_r
      call new_formula_item(form_r_pnt,command_set_target_init,
     &     idx_rtemp)
      form_r_pnt => form_r_pnt%next
      call expand_op_product(form_r_pnt,idx_rtemp,
     &     1d0,1,idxrint,-1,-1,
     &     0,0,.false.,op_info)

      call init_formula(form_gr_temp_3)
      form_gr_temp_3_pnt => form_gr_temp_3
      call new_formula_item(form_gr_temp_3_pnt,
     &     command_set_target_init,idxlcc)
      form_gr_temp_3_pnt => form_gr_temp_3_pnt%next
      ! Call routine which expands the intermediate in terms of available 
      ! integrals. It also forms 'dodgy' contractions, i.e. those where
      ! the arcs form between non-matching faces of the operators. These 
      ! are needed to properly evaluate the intermediates.
      call expand_op_product(form_gr_temp_3_pnt,idxlcc,
     &     1d0,3,(/idx_gtemp,idx_ctemp,idx_rtemp/),
     &     (/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3,1,2,2,3/),3,.true.,op_info)

c dbg
      call print_form_list(luout,form_gr_temp_3,op_info)
      stop
c dbg

      ! Replace G with H and R with R12-INT. Then form the derivative.
      call expand_subexpr(form_gr_temp_3,form_h,op_info)
      call expand_subexpr(form_gr_temp_3,form_r,op_info)      
      
      call form_deriv3(form_v,form_gr_temp_3,
     &     1,idx_ctemp,0,idx_v3,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'V-Temp 3 Formula')
        call print_form_list(luout,form_v,op_info)
      endif  
      stop


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
      call del_operator(idx_v3,op_info)
      call del_operator(idx_v2,op_info)
      call del_operator(idx_v1,op_info)
      call del_operator(idx_v,op_info)
      call del_operator(idx_ftemp,op_info)
      call del_operator(idx_ctemp,op_info)
      call del_operator(idx_sbar,op_info)
      call del_operator(idx_sop,op_info)
      call mem_map(.true.)

      call atim_csw(cpu,sys,wall)
      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'CC-R12 Lagrangian',cpu-cpu0,sys-sys0,wall-wall0)

c      call quit(1,'test','exit')

      return
      end
