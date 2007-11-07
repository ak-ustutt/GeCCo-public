      subroutine form_inter_v(form_lag,idxham,idxc12,idxr12,idxlcc,
     &     idxrint,op_info,orb_info)
*-----------------------------------------------------------------------
*     Subroutine which forms the V-intermediate required for R12 
*     calculations.
*     GWR Sept. 2007
*-----------------------------------------------------------------------
      implicit none

      integer, parameter ::
     &     ntest=100

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

      type(formula_item),intent(inout) ::
     &     form_lag
      integer, intent(in) ::
     &     idxham, idxc12, idxr12, idxlcc, idxrint
      type(operator_info),intent(inout) ::
     &     op_info
      type(orbinf),intent(inout) ::
     &     orb_info

      ! Local variables.
      integer ::
     &     idxopv, idx_ctemp, idx_rtemp, idx_gtemp, ndef, idx_delta
      integer, allocatable ::
     &     occ_def(:,:,:)

      type(formula_item), target ::
     &     form_v, form_gr, form_v_tot, form_h, form_r, form_gr_temp
      type(formula_item), pointer ::
     &     form_gr_pnt, form_v_tot_pnt, form_h_pnt, form_r_pnt,
     &     form_gr_temp_pnt
      
      type(operator), pointer ::
     &     g_temp_pnt, r_temp_pnt, ctemp_pnt
      type(operator_array), pointer ::
     &     ops_array(:)

      integer, external ::
     &     idx_oplist2


      if(.not.explicit)call quit(1,'form_inter_v','R12 not requested')

      if(ntest.ge.100)then
        write(luout,*)'Formation of V-intermediate.'
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

      call form_deriv3(form_v,form_gr,
     &     1,idxc12,0,idxopv,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'V Formula')
        call print_form_list(luout,form_v,op_info)
      endif  

      call factor_out_subexpr(form_lag,form_v,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'Factored R12 Lagrangian')
        call print_form_list(luout,form_lag,op_info)
      endif  

      ! Define a temporary operator, cloned from C12. This is used as
      ! a dummy to allow the parts of V to be formed correctly.
      call add_operator(op_c_temp,op_info)
      idx_ctemp=idx_oplist2(op_c_temp,op_info)
      ctemp_pnt => op_info%op_arr(idx_ctemp)%op
      call clone_operator(ctemp_pnt,op_info%op_arr(idxc12)%op,orb_info)

      ! Generate the derived terms needed to evaluate V.

      ! Initialise the formula in which to hold these terms.
      call init_formula(form_v_tot)
      form_v_tot_pnt => form_v_tot
      call new_formula_item(form_v_tot_pnt,command_set_target_init,
     &     idxopv)
      form_v_tot_pnt => form_v_tot_pnt%next


      ! Part 1.
      call init_formula(form_gr_temp)
      form_gr_temp_pnt => form_gr_temp
      call new_formula_item(form_gr_temp_pnt,
     &     command_set_target_init,idxlcc)
      form_gr_temp_pnt => form_gr_temp_pnt%next

      idx_delta=idx_oplist2(op_del_inter,op_info)

      call inter_full_contr(form_gr_temp_pnt,idxlcc,
     &     1d0,2,(/idx_delta,idx_ctemp/),op_info)

      ! Temporary operator 2.
      call add_operator(op_g_temp,op_info)
      idx_gtemp=idx_oplist2(op_g_temp,op_info)
      g_temp_pnt => op_info%op_arr(idx_gtemp)%op
      ndef=2
      allocate(occ_def(ngastp,2,2))
      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/0,1,0,1/)
      occ_def(1:ngastp,1,2) = (/2,0,0,0/)
      occ_def(1:ngastp,2,2) = (/1,0,0,1/)
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
      ndef=2
      allocate(occ_def(ngastp,2,2))
      occ_def(1:ngastp,1,1) = (/0,1,0,1/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/1,0,0,1/)
      occ_def(1:ngastp,2,2) = (/2,0,0,0/)
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

      form_gr_temp_pnt => form_gr_temp
      do while(associated(form_gr_temp_pnt%next))
        form_gr_temp_pnt => form_gr_temp_pnt%next
      enddo  
      ! Call routine which expands the intermediate in terms of available 
      ! integrals. It also forms 'dodgy' contractions, i.e. those where
      ! the arcs form between non-matching faces of the operators. These 
      ! are needed to properly evaluate the intermediates.
      call expand_op_product(form_gr_temp_pnt,idxlcc,
     &     -1d0,3,(/idx_gtemp,idx_ctemp,idx_rtemp/),
     &     (/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3,1,2,2,3/),3,.true.,op_info)

      call expand_subexpr(form_gr_temp,form_h,.true.,op_info)
      call expand_subexpr(form_gr_temp,form_r,.true.,op_info)

      call del_operator(idx_rtemp,op_info)
      call del_operator(idx_gtemp,op_info)
      call dealloc_formula_list(form_h)
      call dealloc_formula_list(form_r)

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

      form_gr_temp_pnt => form_gr_temp
      do while(associated(form_gr_temp_pnt%next))
        form_gr_temp_pnt => form_gr_temp_pnt%next
      enddo  
      ! Call routine which expands the intermediate in terms of available 
      ! integrals. It also forms 'dodgy' contractions, i.e. those where
      ! the arcs form between non-matching faces of the operators. These 
      ! are needed to properly evaluate the intermediates.
      call expand_op_product(form_gr_temp_pnt,idxlcc,
     &     -0.5d0,3,(/idx_gtemp,idx_ctemp,idx_rtemp/),
     &     (/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3,1,2,2,3/),3,.true.,op_info)

      ! Replace G with H and R with R12-INT. Then form the derivative.
      call expand_subexpr(form_gr_temp,form_h,.true.,op_info)
      call expand_subexpr(form_gr_temp,form_r,.true.,op_info)      
      
      call form_deriv3(form_v_tot,form_gr_temp,
     &     1,idx_ctemp,0,idxopv,op_info)

      if(ntest.ge.100)then
        call print_form_list(luout,form_v_tot,op_info)
      endif  

c      call expand_subexpr(form_v,form_v_tot,op_info)

      call del_operator(idx_rtemp,op_info)
      call del_operator(idx_gtemp,op_info)
      call del_operator(idx_ctemp,op_info)
      call mem_map(.true.)

c      stop

      return
      end
