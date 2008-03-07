      subroutine set_mp2_r12_lagrangian(form_mpr12,
     &     title,label,nlabels,
     &     op_info,orb_info)
*-----------------------------------------------------------------------
*     Routine to set up the MP2-R12 Lagrangian.
*-----------------------------------------------------------------------
      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, intent(in) ::
     &     nlabels
      character(*), intent(in) ::
     &     label(nlabels), title

      type(formula), intent(inout), target ::
     &     form_mpr12

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf), intent(inout) ::
     &     orb_info

      ! local constants
      character(3), parameter ::
     &     op_h_temp = '_H_',
     &     op_sop    = '_S_',
     &     op_sba    = '_SB'

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     form_lag, form_h, form_t_cr, form_tbar_cbarr
      type(formula_item), pointer ::
     &     form_pnt, form_h_pnt, fl_t_cr_pnt
      integer ::
     &     idxham,idxtbar,idxtop,idxlag,idxrba,idxcba,idxr12,idxc12
      integer ::
     &     nterms, idx_h_temp, idx_sop, idx_sbar, ndef, ilabel, idx
      integer, allocatable ::
     &     occ_def(:,:,:)

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

      ! get indices
      if (nlabels.ne.8) then
        write(luout,*) 'nlabels = ',nlabels
        call quit(1,'set_mp2_r12_lagrangian',
     &     'I expect exactly 8 labels')
      end if
      do ilabel = 1, nlabels
        idx = idx_oplist2(label(ilabel),op_info)
        if (idx.le.0)
     &       call quit(1,'set_mp2_r12_lagrangian',
     &       'label not on list: '//trim(label(ilabel)))
        if (ilabel.eq.1) idxlag = idx
        if (ilabel.eq.2) idxham = idx
        if (ilabel.eq.3) idxr12 = idx
        if (ilabel.eq.4) idxrba = idx
        if (ilabel.eq.5) idxtbar = idx
        if (ilabel.eq.6) idxcba = idx
        if (ilabel.eq.7) idxtop = idx
        if (ilabel.eq.8) idxc12 = idx
      end do

      ! Add the parts of the Hamiltonian that are required.
      call add_operator(op_h_temp,op_info)
      idx_h_temp = idx_oplist2(op_h_temp,op_info)
      h_temp_pnt => op_info%op_arr(idx_h_temp)%op

      ndef = 11
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
      occ_def(1:ngastp,1,10) = (/0,1,0,0/)
      occ_def(1:ngastp,2,10) = (/0,0,0,1/)
      occ_def(1:ngastp,1,11) = (/0,0,0,1/)
      occ_def(1:ngastp,2,11) = (/0,1,0,0/)
      call set_uop(h_temp_pnt,op_h_temp,.false.,
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
      if (ntest.ge.1000) then
        write(luout,*) 'H list:'
        call print_form_list(luout,form_h,op_info)
      end if

      ! Add the formal S=T2+CR2 operator.
      call add_operator(op_sop,op_info)
      idx_sop = idx_oplist2(op_sop,op_info)
      sop_pnt => op_info%op_arr(idx_sop)%op

      ! Form T2.
      ndef = 1
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/0,2,0,0/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      call set_uop(sop_pnt,op_sop,.false.,
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

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'T2 + CR2')
        call print_form_list(luout,form_t_cr,op_info)
      end if

      ! Define Sbar for the projection.
      call add_operator(op_sba,op_info)
      idx_sbar = idx_oplist2(op_sba,op_info)
      sbar_pnt => op_info%op_arr(idx_sbar)%op
      call clone_operator(sbar_pnt,sop_pnt,.true.,orb_info)
c      sbar_pnt%dagger = .true.

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

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'T2BAR + CR2BAR')
        call print_form_list(luout,form_tbar_cbarr,op_info)
      end if

      ! The actual formula.
      call init_formula(form_lag)
      form_pnt => form_lag
      call new_formula_item(form_pnt,command_set_target_init,idxlag)
      form_pnt => form_pnt%next

      ! Expand <0|(1+Sbar) e^{-S} H e^S|0> =
      ! <0| e^{-S} H e^S|0> + 
      call expand_op_product(form_pnt,idxlag,
     &     1d0,2,(/idx_h_temp,idx_sop/),-1,-1,
     &     (/1,2/),1,.false.,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'raw formula 1')
        call print_form_list(luout,form_lag,op_info)
      endif

      ! Advance.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! <0|Sbar e^{-S} H e^S|0>
      call expand_op_bch(form_pnt,1,idxlag,
     &     1d0,idx_sbar,idx_h_temp,1d0,idx_sop,1,-1,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'raw formula')
        call print_form_list(luout,form_lag,op_info)
      endif

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
      if(ntest.ge.1000)then
        call write_title(luout,wst_title,'After S+ replacement.')
        call print_form_list(luout,form_lag,op_info)
      endif

      ! this runs counter the philosophy of the approach
c      ! Add extra, artificial terms for approximations beyond A.
c      if(trim(r12_apprx).ne.'A')then
c
c        call del_operator(idx_h_temp,op_info)
c        call add_operator(op_h_temp,op_info)
c        idx_h_temp = idx_oplist2(op_h_temp,op_info)
c        h_temp_pnt => op_info%op_arr(idx_h_temp)%op
c
c        ndef = 1
c        allocate(occ_def(ngastp,2,ndef))
c        occ_def(1:ngastp,1,1) = (/1,0,0,0/)
c        occ_def(1:ngastp,2,1) = (/1,0,0,0/)
c        call set_uop(h_temp_pnt,op_h_temp,.false.,
c     &       occ_def,ndef,orb_info)
c        deallocate(occ_def)
c        
c       ! Replace the formal terms with elements of H.
c        call dealloc_formula_list(form_h)
c        call init_formula(form_h)
c        form_h_pnt => form_h
c        call new_formula_item(form_h_pnt,
c     &       command_set_target_init,idx_h_temp)
c        form_h_pnt => form_h_pnt%next
c        call expand_op_product(form_h_pnt,idx_h_temp,
c     &       1d0,1,idxham,-1,-1,
c     &       0,0,.false.,op_info)
c
c        do while(associated(form_pnt%next))
c          form_pnt => form_pnt%next
c        enddo
c
c        call expand_op_product(form_pnt,idxlag,
c     &       -0.5d0,5,(/idxrba,idx_h_temp,idxcba,idxc12,idxr12/),-1,-1,
c     &       (/1,2,2,3,1,3,1,5,3,4,4,5/),6,.false.,op_info)
c
c        do while(associated(form_pnt%next))
c          form_pnt => form_pnt%next
c        enddo
c
c        call expand_op_product(form_pnt,idxlag,
c     &       -0.5d0,5,(/idxrba,idxcba,idxc12,idx_h_temp,idxr12/),-1,-1,
c     &       (/1,2,2,3,1,5,3,4,3,5,4,5/),6,.false.,op_info)
c
c        call expand_subexpr(form_lag,form_h,.false.,op_info)
c      endif
c
      if(ntest.ge.100)then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,form_lag,op_info)
      endif

      ! Assign canonical name and comment.
      form_mpr12%comment = title
      ! Write to disc.
      write(name,'(a,".fml")') trim(form_mpr12%label)
      call file_init(form_mpr12%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_mpr12%fhand,form_lag,title)

      ! Delete linked lists (!)
      call dealloc_formula_list(form_h)
      call dealloc_formula_list(form_t_cr)
      call dealloc_formula_list(form_tbar_cbarr)
      call dealloc_formula_list(form_lag)

      ! Delete the temporary operators.
      call del_operator(op_h_temp,op_info)
      call del_operator(op_sba,op_info)
      call del_operator(op_sop,op_info)

      call atim_csw(cpu,sys,wall)
      if(ntest.ge.100) then
        write(luout,*) 'Number of generated terms: ',nterms
      endif
      call prtim(luout,'MP2-R12 Lagrangian',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end


