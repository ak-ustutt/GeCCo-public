*----------------------------------------------------------------------*
      subroutine set_b_intermediate(formula_bint,op_info,orb_info)
*----------------------------------------------------------------------*
*     Generate the formula for the B-intermediate (suitable for 
*     MP2-R12 variant 1A).
*     GWR November 2007
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

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
     &     formula_bint

      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_bint
      type(formula_item), pointer ::
     &     form_pnt

      type(operator), pointer ::
     &     ttr_temp_pnt, rbar_temp_pnt

      type(formula_item), target ::
     &     form_ttr, form_rbar, form_rttr_temp, form_rttr, form_b_tot
      type(formula_item), pointer ::
     &     form_ttr_pnt, form_rbar_pnt, form_rttr_temp_pnt,
     &     form_rttr_pnt, form_b_tot_pnt

      integer ::
     &     nterms, idx_bint, idxunity, ndef, idx_lcc, idx_ttr,
     &     idx_rbar, idx_c12, idx_ctemp, idx_rbar_temp, idx_ttr_temp

      integer, allocatable ::
     &     occ_def(:,:,:)

      integer,external ::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.ge.100) then
        write(luout,*) '==================================='
        write(luout,*) ' output from set_b_intermediate    '
        write(luout,*) '==================================='
      end if

      call atim_csw(cpu0,sys0,wall0)

      idx_lcc = idx_oplist2(op_cclg,op_info)
      idx_ttr = idx_oplist2(op_ttr,op_info)
      idx_rbar = idx_oplist2(op_rinba,op_info)
      idx_c12 = idx_oplist2(op_c12,op_info)

      ! Initialize the formula required for the whole B-intermediate.
      call init_formula(form_bint)
      form_pnt => form_bint
      idx_bint = idx_oplist2(op_b_inter,op_info)


      ! Form the first set of contracted integrals.
      call init_formula(form_rttr_temp)
      form_rttr_temp_pnt => form_rttr_temp
      call new_formula_item(form_rttr_temp_pnt,
     &     command_set_target_init,idx_lcc)
      form_rttr_temp_pnt => form_rttr_temp_pnt%next

      ! Set up a temporary operator to represent R+.
      call add_operator(op_rbar_temp,op_info)
      idx_rbar_temp = idx_oplist2(op_rbar_temp,op_info)
      rbar_temp_pnt => op_info%op_arr(idx_rbar_temp)%op

      ndef = 2
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/0,1,0,1/)
      occ_def(1:ngastp,1,2) = (/2,0,0,0/)
      occ_def(1:ngastp,2,2) = (/1,0,0,1/)
      call set_uop(rbar_temp_pnt,op_rbar_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Substitute with actual R-BAR integrals.
      call init_formula(form_rbar)
      form_rbar_pnt => form_rbar
      call new_formula_item(form_rbar_pnt,command_set_target_init,
     &     idx_rbar_temp)
      form_rbar_pnt = > form_rbar_pnt%next
      call expand_op_product(form_rbar_pnt,idx_rbar_temp,
     &     1d0,1,idx_rbar,-1,-1,
     &     0,0,.false.,op_info)

      ! Set up a temporary operator to represent [T1+T2,r12]
      call add_operator(op_ttr_temp,op_info)
      idx_ttr_temp = idx_oplist2(op_ttr_temp,op_info)
      ttr_temp_pnt => op_info%op_arr(idx_ttr_temp)%op

      ndef = 2
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/0,1,0,1/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/1,0,0,1/)
      occ_def(1:ngastp,2,2) = (/2,0,0,0/)
      call set_uop(ttr_temp_pnt,op_ttr_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Substitute with actual T commutator integrals.
      call init_formula(form_ttr)
      form_ttr_pnt => form_ttr
      call new_formula_item(form_ttr_pnt,command_set_target_init,
     &     idx_ttr_temp)
      form_ttr_pnt = > form_ttr_pnt%next
      call expand_op_product(form_ttr_pnt,idx_ttr_temp,
     &     1d0,1,idx_ttr,-1,-1,
     &     0,0,.false.,op_info)

      ! Combine the necessary integrals into the term required for the
      ! intermediate.
      call expand_op_product(form_rttr_temp_pnt,idx_lcc,
c     &     -1d0,3,(/idx_rbar_temp,idx_c12,idx_ttr_temp/),
c dbg
     &     1d0,3,(/idx_rbar_temp,idx_c12,idx_ttr_temp/),
c dbg
     &     (/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3,1,2,2,3/),3,.true.,op_info)

      call expand_subexpr(form_rttr_temp,form_rbar,.true.,op_info)
      call expand_subexpr(form_rttr_temp,form_ttr,.true.,op_info)

      call del_operator(idx_ttr_temp,op_info)
      call del_operator(idx_rbar_temp,op_info)
      call dealloc_formula_list(form_rbar)
      call dealloc_formula_list(form_ttr)


      ! Form the second set of contracted integrals.

      ! Set up a temporary operator to represent R+.
      call add_operator(op_rbar_temp,op_info)
      idx_rbar_temp = idx_oplist2(op_rbar_temp,op_info)
      rbar_temp_pnt => op_info%op_arr(idx_rbar_temp)%op

      ndef = 3
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/2,0,0,0/)
      occ_def(1:ngastp,2,2) = (/1,1,0,0/)
      occ_def(1:ngastp,1,3) = (/2,0,0,0/)
      occ_def(1:ngastp,2,3) = (/0,2,0,0/)
      call set_uop(rbar_temp_pnt,op_rbar_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Substitute with actual R-BAR integrals.
      call init_formula(form_rbar)
      form_rbar_pnt => form_rbar
      call new_formula_item(form_rbar_pnt,command_set_target_init,
     &     idx_rbar_temp)
      form_rbar_pnt = > form_rbar_pnt%next
      call expand_op_product(form_rbar_pnt,idx_rbar_temp,
     &     1d0,1,idx_rbar,-1,-1,
     &     0,0,.false.,op_info)

      ! Set up a temporary operator to represent [T1+T2,r12]
      call add_operator(op_ttr_temp,op_info)
      idx_ttr_temp = idx_oplist2(op_ttr_temp,op_info)
      ttr_temp_pnt => op_info%op_arr(idx_ttr_temp)%op

      ndef = 3
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/1,1,0,0/)
      occ_def(1:ngastp,2,2) = (/2,0,0,0/)
      occ_def(1:ngastp,1,3) = (/0,2,0,0/)
      occ_def(1:ngastp,2,3) = (/2,0,0,0/)
      call set_uop(ttr_temp_pnt,op_ttr_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Substitute with actual R-BAR integrals.
      call init_formula(form_ttr)
      form_ttr_pnt => form_ttr
      call new_formula_item(form_ttr_pnt,command_set_target_init,
     &     idx_ttr_temp)
      form_ttr_pnt = > form_ttr_pnt%next
      call expand_op_product(form_ttr_pnt,idx_ttr_temp,
     &     1d0,1,idx_ttr,-1,-1,
     &     0,0,.false.,op_info)

      ! Combine the necessary integrals into the term required for the
      ! intermediate.
      form_rttr_temp_pnt = > form_rttr_temp
      do while(associated(form_rttr_temp_pnt%next))
        form_rttr_temp_pnt => form_rttr_temp_pnt%next
      enddo
      call expand_op_product(form_rttr_temp_pnt,idx_lcc,
c     &     -1d0,3,(/idx_rbar_temp,idx_c12,idx_ttr_temp/),
c dbg
     &     1d0,3,(/idx_rbar_temp,idx_c12,idx_ttr_temp/),
c dbg
     &     (/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3,1,2,2,3/),3,.true.,op_info)

      call expand_subexpr(form_rttr_temp,form_rbar,.true.,op_info)
      call expand_subexpr(form_rttr_temp,form_ttr,.true.,op_info)

      ! Form the derivative of the first two terms wrt C12.
      call form_deriv3(form_pnt,form_rttr_temp,
     &     1,idx_c12,0,idx_bint,op_info)

      call del_operator(idx_ttr_temp,op_info)
      call del_operator(idx_rbar_temp,op_info)
      call dealloc_formula_list(form_rbar)
      call dealloc_formula_list(form_ttr)


      ! Add the formal delta operator.
      do while (associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo  

      idxunity = idx_oplist2(op_unity,op_info)
      ! nvtx = 1
      call new_formula_item(form_pnt,command_add_contribution,idx_bint)
      call resize_contr(form_pnt%contr,1,0,0)
c      form_pnt%contr%fac = 1d0
c dbg
      form_pnt%contr%fac = -1d0
c dbg
      form_pnt%contr%nvtx = 1
      form_pnt%contr%nsupvtx = 1
      form_pnt%contr%idx_res = idx_bint
      form_pnt%contr%iblk_res = 1
      form_pnt%contr%vertex(1)%idx_op = idxunity
      form_pnt%contr%vertex(1)%iblk_op = 1
      call update_svtx4contr(form_pnt%contr)

      if(ntest.ge.100)then
        write(luout,*) 'B-intermediate: '
        call print_form_list(luout,form_bint,op_info)
      endif

      ! write to disc
      formula_bint%label = label_r12_bint
      formula_bint%comment = title_r12_bint
      write(name,'(a,".fml")') label_r12_bint
      call file_init(formula_bint%fhand,name,ftyp_sq_unf,0)
      call write_form_list(formula_bint%fhand,form_bint,title_r12_bint)

      call dealloc_formula_list(form_bint)


      call atim_csw(cpu,sys,wall)
      call prtim(luout,'B-interm.',cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
