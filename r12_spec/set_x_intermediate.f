*----------------------------------------------------------------------*
      subroutine set_x_intermediate(formula_xint,op_info,orb_info)
*----------------------------------------------------------------------*
*     Generate the formula for the X-intermediate (suitable for 
*     MP2-R12 variants 1A' and 1B).
*     GWR December 2007.
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
     &     formula_xint

      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      ! Local variables:
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_xint, form_rbarr_temp, form_rbar, form_r
      type(formula_item), pointer ::
     &     form_pnt, form_rbarr_temp_pnt, form_rbar_pnt, form_r_pnt

      type(operator), pointer ::
     &     r_temp_pnt, rbar_temp_pnt


      integer ::
     &     idxlcc, idxrint, idxrbar, idxc12, idx_xint, idx_rbar_temp,
     &     ndef, idx_r_temp, idxr12sq
      integer, allocatable ::
     &     occ_def(:,:,:)
      integer, external ::
     &     idx_oplist2

      ! For timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if(ntest.ge.100)then
        write(luout,*)'================================'
        write(luout,*)' Output from set_x_intermediate '
        write(luout,*)'================================'
      endif

      call atim_csw(cpu0,sys0,wall0)

      idxlcc = idx_oplist2(op_cclg,op_info)
      idxrint = idx_oplist2(op_rint,op_info)
      idxrbar = idx_oplist2(op_rinba,op_info)
      idxc12 = idx_oplist2(op_c12,op_info)

      ! Initialise the formula required for the X-intermediate.
      call init_formula(form_xint)
      form_pnt => form_xint
      idx_xint = idx_oplist2(op_x_inter,op_info)

      ! Form the first set of contracted integrals.

      ! Initialise a temporary formula.
      call init_formula(form_rbarr_temp)
      form_rbarr_temp_pnt => form_rbarr_temp
      call new_formula_item(form_rbarr_temp_pnt,
     &     command_set_target_init,idxlcc)
      form_rbarr_temp_pnt => form_rbarr_temp_pnt%next

      ! Form a temporary operator to represent R+.
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

      ! Substitute with the actual R-BAR integrals.
      call init_formula(form_rbar)
      form_rbar_pnt => form_rbar
      call new_formula_item(form_rbar_pnt,command_set_target_init,
     &     idx_rbar_temp)
      form_rbar_pnt = > form_rbar_pnt%next
      call expand_op_product(form_rbar_pnt,idx_rbar_temp,
     &     1d0,1,idxrbar,-1,-1,
     &     0,0,.false.,op_info)

      ! Form a temporary operator to represent R.
      call add_operator(op_r_temp,op_info)
      idx_r_temp = idx_oplist2(op_r_temp,op_info)
      r_temp_pnt => op_info%op_arr(idx_r_temp)%op

      ndef = 2
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/0,1,0,1/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/1,0,0,1/)
      occ_def(1:ngastp,2,2) = (/2,0,0,0/)
      call set_uop(r_temp_pnt,op_r_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Substitute with the actual R-integrals.
      call init_formula(form_r)
      form_r_pnt => form_r
      call new_formula_item(form_r_pnt,command_set_target_init,
     &     idx_r_temp)
      form_r_pnt = > form_r_pnt%next
      call expand_op_product(form_r_pnt,idx_r_temp,
     &     1d0,1,idxrint,-1,-1,
     &     0,0,.false.,op_info)

      ! Combine the necessary integrals into the term required for the
      ! intermediate.
      call expand_op_product(form_rbarr_temp_pnt,idxlcc,
     &     -1d0,3,(/idx_rbar_temp,idxc12,idx_r_temp/),
     &     (/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3,1,2,2,3/),3,.true.,op_info)

      call expand_subexpr(form_rbarr_temp,form_rbar,.true.,op_info)
      call expand_subexpr(form_rbarr_temp,form_r,.true.,op_info)

      ! Delete temporary items.
      call del_operator(idx_r_temp,op_info)
      call del_operator(idx_rbar_temp,op_info)
      call dealloc_formula_list(form_r)
      call dealloc_formula_list(form_rbar)


      ! Form the second set of contracted integrals.

      ! Form a temporary operator to represent R+.
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

      ! Substitute with the actual R-BAR integrals.
      call init_formula(form_rbar)
      form_rbar_pnt => form_rbar
      call new_formula_item(form_rbar_pnt,command_set_target_init,
     &     idx_rbar_temp)
      form_rbar_pnt = > form_rbar_pnt%next
      call expand_op_product(form_rbar_pnt,idx_rbar_temp,
     &     1d0,1,idxrbar,-1,-1,
     &     0,0,.false.,op_info)

      ! Form a temporary operator to represent R.
      call add_operator(op_r_temp,op_info)
      idx_r_temp = idx_oplist2(op_r_temp,op_info)
      r_temp_pnt => op_info%op_arr(idx_r_temp)%op

      ndef = 3
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/1,1,0,0/)
      occ_def(1:ngastp,2,2) = (/2,0,0,0/)
      occ_def(1:ngastp,1,3) = (/0,2,0,0/)
      occ_def(1:ngastp,2,3) = (/2,0,0,0/)
      call set_uop(r_temp_pnt,op_r_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Substitute with the actual R-integrals.
      call init_formula(form_r)
      form_r_pnt => form_r
      call new_formula_item(form_r_pnt,command_set_target_init,
     &     idx_r_temp)
      form_r_pnt = > form_r_pnt%next
      call expand_op_product(form_r_pnt,idx_r_temp,
     &     1d0,1,idxrint,-1,-1,
     &     0,0,.false.,op_info)

      ! Combine the necessary integrals into the term required for the
      ! intermediate.
      form_rbarr_temp_pnt = > form_rbarr_temp
      do while(associated(form_rbarr_temp_pnt%next))
        form_rbarr_temp_pnt => form_rbarr_temp_pnt%next
      enddo
      call expand_op_product(form_rbarr_temp_pnt,idxlcc,
     &     -1d0,3,(/idx_rbar_temp,idxc12,idx_r_temp/),
     &     (/-1,-1,-1/),(/-1,-1,-1/),
     &     (/1,3,1,2,2,3/),3,.true.,op_info)

      call expand_subexpr(form_rbarr_temp,form_rbar,.true.,op_info)
      call expand_subexpr(form_rbarr_temp,form_r,.true.,op_info)

      ! Delete temporary items.
      call del_operator(idx_r_temp,op_info)
      call del_operator(idx_rbar_temp,op_info)
      call dealloc_formula_list(form_r)
      call dealloc_formula_list(form_rbar)


      ! Form the derivative of the first two terms wrt C12.
      call form_deriv3(form_pnt,form_rbarr_temp,
     &     1,idxc12,0,idx_xint,op_info)

      ! Add the R12^2 operator.
      do while (associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo  

      idxr12sq = idx_oplist2(op_f2,op_info)
      ! nvtx = 1
      call new_formula_item(form_pnt,command_add_contribution,idx_xint)
      call resize_contr(form_pnt%contr,1,0,0,0)
      form_pnt%contr%fac = 1d0
      form_pnt%contr%nvtx = 1
      form_pnt%contr%nsupvtx = 1
      form_pnt%contr%idx_res = idx_xint
      form_pnt%contr%iblk_res = 1
      form_pnt%contr%vertex(1)%idx_op = idxr12sq
      form_pnt%contr%vertex(1)%iblk_op = 1
      call update_svtx4contr(form_pnt%contr)


      ! Print out final result if required.
      if(ntest.ge.100)then
        write(luout,*) 'X-intermediate: '
        call print_form_list(luout,form_xint,op_info)
      endif

      ! Write formula to disc.
      formula_xint%label = label_r12_xint
      formula_xint%comment = title_r12_xint
      write(name,'(a,".fml")') label_r12_xint
      call file_init(formula_xint%fhand,name,ftyp_sq_unf,0)
      call write_form_list(formula_xint%fhand,form_xint,title_r12_xint)

      call dealloc_formula_list(form_xint)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'X-interm.',cpu-cpu0,sys-sys0,wall-wall0)

c dbg
c      stop
c dbg

      return
      end


