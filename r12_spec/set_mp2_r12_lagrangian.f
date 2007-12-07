      subroutine set_mp2_r12_lagrangian(form_cclag,op_info,orb_info,
     &     idxham,idxtbar,idxrba,idxcba,idxtop,idxr12,idxc12,idxlcc)
*-----------------------------------------------------------------------
*     Routine to set up the MP2-R12 Lagrangian.
*-----------------------------------------------------------------------
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
     &     form_lag, form_h, form_t_cr, form_tbar_cbarr
      type(formula_item), pointer ::
     &     form_pnt, form_h_pnt, fl_t_cr_pnt
      integer ::
     &     nterms, idx_h_temp, idx_sop, idx_sbar, ndef, idxrint
      integer, allocatable ::
     &     occ_def(:,:,:)
      type(operator_array), pointer ::
     &     ops_array(:)

      type(operator), pointer::
     &     h_temp_pnt, sop_pnt, sbar_pnt

      integer, external::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(luout,*) '==================================='
        write(luout,*) ' output from set_mp2_r12_lagrangian'
        write(luout,*) '==================================='
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! Add the parts of the Hamiltonian that are required.
      call add_operator(op_h_temp,op_info)
      idx_h_temp = idx_oplist2(op_h_temp,op_info)
      h_temp_pnt => op_info%op_arr(idx_h_temp)%op

      ndef =9
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/1,0,0,0/)
      occ_def(1:ngastp,2,1) = (/1,0,0,0/)
      occ_def(1:ngastp,1,2) = (/0,1,0,0/)
      occ_def(1:ngastp,2,2) = (/0,1,0,0/)
      occ_def(1:ngastp,1,3) = (/0,0,0,1/)
      occ_def(1:ngastp,2,3) = (/0,0,0,1/)
      occ_def(1:ngastp,1,4) = (/0,2,0,0/)
      occ_def(1:ngastp,2,4) = (/2,0,0,0/)
      occ_def(1:ngastp,1,5) = (/0,1,0,1/)
      occ_def(1:ngastp,2,5) = (/2,0,0,0/)
      occ_def(1:ngastp,1,6) = (/0,0,0,2/)
      occ_def(1:ngastp,2,6) = (/2,0,0,0/)
      occ_def(1:ngastp,1,7) = (/2,0,0,0/)
      occ_def(1:ngastp,2,7) = (/0,2,0,0/)
      occ_def(1:ngastp,1,8) = (/2,0,0,0/)
      occ_def(1:ngastp,2,8) = (/0,1,0,1/)
      occ_def(1:ngastp,1,9) = (/2,0,0,0/)
      occ_def(1:ngastp,2,9) = (/0,0,0,2/)
      call set_uop(h_temp_pnt,op_h_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Replace the formal terms with elements of H.
      call init_formula(form_h)
      form_h_pnt => form_h
      call new_formula_item(form_h_pnt,
     &     command_set_target_init,idx_h_temp)
      form_h_pnt => form_h_pnt%next
      call expand_op_product(form_h_pnt,idx_h_temp,
     &     1d0,1,idxham,-1,-1,
     &     0,0,.false.,op_info)
      call print_form_list(luout,form_h,op_info)


      ! Add the formal S=T2+CR2 operator.
      call add_operator(op_sop,op_info)
      idx_sop = idx_oplist2(op_sop,op_info)
      sop_pnt => op_info%op_arr(idx_sop)%op

      if(trir12.ne.0)
     &     call quit(1,'set_mp2_r12_lagrangian','undefined for triples')

      ! Form T2.
      ndef = 1
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/0,2,0,0/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      call set_uop(sop_pnt,op_sop,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Add CR part.
      call join_operator(sop_pnt,op_info%op_arr(idxr12)%op,orb_info)

      ! Replace the formal terms with the predefined operators.
      call init_formula(form_t_cr)
      fl_t_cr_pnt => form_t_cr
      call new_formula_item(fl_t_cr_pnt,command_set_target_init,idx_sop)
      fl_t_cr_pnt => fl_t_cr_pnt%next
      call expand_op_product(fl_t_cr_pnt,idx_sop,
     &     1d0,1,idxtop,-1,-1,
     &     0,0,.false.,op_info)
      do while(associated(fl_t_cr_pnt%next))
        fl_t_cr_pnt => fl_t_cr_pnt%next
      enddo

      call expand_op_product(fl_t_cr_pnt,idx_sop,
     &     1d0,2,(/idxc12,idxr12/),-1,-1,
     &     (/1,2/),1,.false.,op_info)

      call write_title(luout,wst_title,'T2 + CR2')
      call print_form_list(luout,form_t_cr,op_info)


      ! Define Sbar for the projection.
      call add_operator(op_sba,op_info)
      idx_sbar = idx_oplist2(op_sba,op_info)
      sbar_pnt => op_info%op_arr(idx_sbar)%op
      call clone_operator(sbar_pnt,sop_pnt,orb_info)
      sbar_pnt%dagger = .true.

      ! Complete Sbar also.
      call init_formula(form_tbar_cbarr)
      fl_t_cr_pnt => form_tbar_cbarr
      call new_formula_item(fl_t_cr_pnt,
     &     command_set_target_init,idx_sbar)
      fl_t_cr_pnt => fl_t_cr_pnt%next
      call expand_op_product(fl_t_cr_pnt,idx_sbar,
     &     1d0,1,idxtbar,-1,-1,
     &     0,0,.false.,op_info)
      do while(associated(fl_t_cr_pnt%next))
        fl_t_cr_pnt => fl_t_cr_pnt%next
      enddo

      call expand_op_product(fl_t_cr_pnt,idx_sbar,
     &     1d0,2,(/idxrba,idxcba/),-1,-1,
     &     (/1,2/),1,.false.,op_info)

      call write_title(luout,wst_title,'T2BAR + CR2BAR')
      call print_form_list(luout,form_tbar_cbarr,op_info)


      ! The actual formula.
      call init_formula(form_lag)
      form_pnt => form_lag
      call new_formula_item(form_pnt,command_set_target_init,idxlcc)
      form_pnt => form_pnt%next

      ! Expand <0|(1+Sbar) e^{-S} H e^S|0> =
      ! <0| e^{-S} H e^S|0> + 
      call expand_op_product(form_pnt,idxlcc,
     &     1d0,2,(/idx_h_temp,idx_sop/),-1,-1,
     &     (/1,2/),1,.false.,op_info)

      ! Advance.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! <0|Sbar e^{-S} H e^S|0>
      call expand_op_bch(form_pnt,1,idxlcc,
     &     1d0,idx_sbar,idx_h_temp,1d0,idx_sop,1,-1,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'raw formula')
        call print_form_list(luout,form_lag,op_info)
      endif

      call sum_terms(form_lag,op_info)
      call cc_form_post(form_lag,nterms,
     &     idx_sbar,idx_h_temp,idx_sop,op_info)

      ! Replace the formal elements of H, S and Sbar with actual
      ! Hamiltonian, T and CR terms.
      call expand_subexpr(form_lag,form_h,.false.,op_info)
      if(ntest.ge.1000)then
        call write_title(luout,wst_title,'After H replacement.')
        call print_form_list(luout,form_lag,op_info)
      endif

      call expand_subexpr(form_lag,form_t_cr,.false.,op_info)
      if(ntest.ge.1000)then
        call write_title(luout,wst_title,'After S replacement.')
        call print_form_list(luout,form_lag,op_info)
      endif

      call expand_subexpr(form_lag,form_tbar_cbarr,.false.,op_info)
      if(ntest.ge.100)then
        call write_title(luout,wst_title,'Final formula.')
        call print_form_list(luout,form_lag,op_info)
      endif

      ! Assign canonical name and comment.
      form_cclag%label = label_cclg0
      form_cclag%comment = title_cclg0
      ! Write to disc.
      write(name,'(a,".fml")') label_cclg0
      call file_init(form_cclag%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_cclag%fhand,form_lag,title_cclg0)

      ! Delete the temporary operators.
      call del_operator(idx_h_temp,op_info)
      call del_operator(idx_sbar,op_info)
      call del_operator(idx_sop,op_info)

      call atim_csw(cpu,sys,wall)
      if(ntest.ge.100) then
        write(luout,*) 'Number of generated terms: ',nterms
      endif
      call prtim(luout,'MP2-R12 Lagrangian',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end


