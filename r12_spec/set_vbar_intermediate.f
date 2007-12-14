*----------------------------------------------------------------------*
      subroutine set_vbar_intermediate(formula_vbint,op_info,orb_info)
*----------------------------------------------------------------------*
*     Generate the formula for the V+-intermediate.
*     GWR November 2007.
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'par_opnames_gen.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'par_formnames_gen.h'

      type(formula), intent(inout), target ::
     &     formula_vbint

      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_vbint
      type(formula_item), pointer ::
     &     form_pnt

      type(operator), pointer ::
     &     g_temp_pnt, r_temp_pnt

      type(formula_item), target ::
     &     form_h, form_r, form_gr_temp, form_gr, form_vb_tot
      type(formula_item), pointer ::
     &     form_h_pnt, form_r_pnt, form_gr_temp_pnt, form_gr_pnt,
     &     form_vb_tot_pnt

      integer ::
     &     nterms, idx_vbint, idxunity, idx_gtemp, idx_rtemp, ndef,
     &     idxham, idxlcc, idx_rint, idxc12

      integer, allocatable ::
     &     occ_def(:,:,:)

      integer,external ::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(luout,*) '==================================='
        write(luout,*) ' output from set_vbar_intermediate '
        write(luout,*) '==================================='
      end if

      call atim_csw(cpu0,sys0,wall0)

      idxham = idx_oplist2(op_ham,op_info)
      idxlcc = idx_oplist2(op_cclg,op_info)
      idx_rint = idx_oplist2(op_rinba,op_info)
      idxc12 = idx_oplist2(op_c12,op_info)

      ! initialize formula
      call init_formula(form_vbint)
      form_pnt => form_vbint

      idx_vbint = idx_oplist2(op_vbar_inter,op_info)

      ! Form the first set of contracted integrals. 
      call init_formula(form_gr_temp)
      form_gr_temp_pnt => form_gr_temp
      call new_formula_item(form_gr_temp_pnt,
     &     command_set_target_init,idxlcc)
      form_gr_temp_pnt => form_gr_temp_pnt%next

      ! Set up a temporary operator to represent g.
      call add_operator(op_g_temp,op_info)
      idx_gtemp=idx_oplist2(op_g_temp,op_info)
      g_temp_pnt => op_info%op_arr(idx_gtemp)%op

      ndef=2
      allocate(occ_def(ngastp,2,2))
      occ_def(1:ngastp,1,1) = (/0,1,0,1/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/1,0,0,1/)
      occ_def(1:ngastp,2,2) = (/2,0,0,0/)
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

      ! Set up a temporary operator to represent R+.
      call add_operator(op_r_temp,op_info)
      idx_rtemp=idx_oplist2(op_r_temp,op_info)
      r_temp_pnt => op_info%op_arr(idx_rtemp)%op
      ndef=2
      allocate(occ_def(ngastp,2,2))
      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/0,1,0,1/)
      occ_def(1:ngastp,1,2) = (/2,0,0,0/)
      occ_def(1:ngastp,2,2) = (/1,0,0,1/)
      call set_uop(r_temp_pnt,op_r_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Substitute with actual R12+ integrals.
      call init_formula(form_r)
      form_r_pnt => form_r
      call new_formula_item(form_r_pnt,command_set_target_init,
     &     idx_rtemp)
      form_r_pnt => form_r_pnt%next
      call expand_op_product(form_r_pnt,idx_rtemp,
     &     1d0,1,idx_rint,-1,-1,
     &     0,0,.false.,op_info)

      ! Call routine which expands the intermediate in terms of available 
      ! integrals. It also forms 'dodgy' contractions, i.e. those where
      ! the arcs form between non-matching faces of the operators. These 
      ! are needed to properly evaluate the intermediates.
      call expand_op_product(form_gr_temp_pnt,idxlcc,
     &     -1d0,3,(/idx_rtemp,idxc12,idx_gtemp/),
c dbg
c     &     1d0,3,(/idx_rtemp,idxc12,idx_gtemp/),
c dbg
     &     (/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3,1,2,2,3/),3,.true.,op_info)

      call expand_subexpr(form_gr_temp,form_h,.true.,op_info)
      call expand_subexpr(form_gr_temp,form_r,.true.,op_info)

      call del_operator(idx_rtemp,op_info)
      call del_operator(idx_gtemp,op_info)
      call dealloc_formula_list(form_h)
      call dealloc_formula_list(form_r)


      ! Temporary operator 2.
      call add_operator(op_g_temp,op_info)
      idx_gtemp=idx_oplist2(op_g_temp,op_info)
      g_temp_pnt => op_info%op_arr(idx_gtemp)%op
      ndef=3
      allocate(occ_def(ngastp,2,3))
      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/1,1,0,0/)
      occ_def(1:ngastp,2,2) = (/2,0,0,0/)
      occ_def(1:ngastp,1,3) = (/0,2,0,0/)
      occ_def(1:ngastp,2,3) = (/2,0,0,0/)
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
      occ_def(1:ngastp,1,2) = (/2,0,0,0/)
      occ_def(1:ngastp,2,2) = (/1,1,0,0/)
      occ_def(1:ngastp,1,3) = (/2,0,0,0/)
      occ_def(1:ngastp,2,3) = (/0,2,0,0/)
      call set_uop(r_temp_pnt,op_r_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Substitute with actual R12+ integrals.
      call init_formula(form_r)
      form_r_pnt => form_r
      call new_formula_item(form_r_pnt,command_set_target_init,
     &     idx_rtemp)
      form_r_pnt => form_r_pnt%next
      call expand_op_product(form_r_pnt,idx_rtemp,
     &     1d0,1,idx_rint,-1,-1,
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
     &     -1d0,3,(/idx_rtemp,idxc12,idx_gtemp/),
c dbg
c     &     1d0,3,(/idx_rtemp,idxc12,idx_gtemp/),
c dbg
     &     (/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3,1,2,2,3/),3,.true.,op_info)

      ! Replace G with H and R with R12-INT. Then form the derivative.
      call expand_subexpr(form_gr_temp,form_h,.true.,op_info)
      call expand_subexpr(form_gr_temp,form_r,.true.,op_info)      

      call form_deriv3(form_pnt,form_gr_temp,
     &     1,idxc12,0,idx_vbint,op_info)

      call del_operator(idx_rtemp,op_info)
      call del_operator(idx_gtemp,op_info)
      call dealloc_formula_list(form_h)
      call dealloc_formula_list(form_r)

      ! Add the formal delta operator.
      do while (associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo  

      idxunity = idx_oplist2(op_unity,op_info)
      ! nvtx = 1
      call new_formula_item(form_pnt,command_add_contribution,idx_vbint)
      call resize_contr(form_pnt%contr,1,0,0,0)
      form_pnt%contr%fac = 1d0
c dbg
c      form_pnt%contr%fac = -1d0
c dbg
      form_pnt%contr%nvtx = 1
      form_pnt%contr%nsupvtx = 1
      form_pnt%contr%idx_res = idx_vbint
      form_pnt%contr%iblk_res = 1
      form_pnt%contr%vertex(1)%idx_op = idxunity
      form_pnt%contr%vertex(1)%iblk_op = 1
      call update_svtx4contr(form_pnt%contr)

      if(ntest.ge.100)then
        write(luout,*) 'V-integrals'
        call print_form_list(luout,form_vbint,op_info)
      endif

      ! write to disc
      formula_vbint%label = label_r12_vbint
      formula_vbint%comment = title_r12_vbint
      write(name,'(a,".fml")') label_r12_vbint
      call file_init(formula_vbint%fhand,name,ftyp_sq_unf,0)
      call write_form_list(formula_vbint%fhand,
     &     form_vbint,title_r12_vbint)

      call dealloc_formula_list(form_vbint)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'V+-interm.',cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
