*----------------------------------------------------------------------*
      subroutine set_ccpt_lagrangian(form_lag,
     &     title,label_op,nlabels,ansatz,mode,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*     set up sequence of operators, integrals and contractions that
*     defines a CC-Lagrangian within the chosen operator space 
*
*     written by andreas, june 2008
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
      include 'ifc_input.h'

      type(formula), intent(inout), target ::
     &     form_lag

      character*(*), intent(in) ::
     &     title,
     &     label_op(*)
      integer, intent(in) ::
     &     nlabels

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf), intent(inout) ::
     &     orb_info

      ! local constants
      character(3), parameter ::
     &     op_f_temp = '_F_',
     &     op_h_temp = '_H_',
     &     op_sop    = '_S_',
     &     op_spt    = '_SP'

      ! local variables

      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist_lag, flist_t_r, flist_t_r_pt
      type(formula_item), pointer ::
     &     flist_pnt
      integer, intent(in) ::
     &     ansatz
      character*(*) ::
     &     mode

      logical ::
     &     r12fix,truncate
      integer ::
     &     nterms, ilabel, idx, ndef, 
     &     idxham,idxtbar,idxtop,idxtpt,idxtptbar,idxlcc,
     &     idxr12,idxc12,idxcpp12,idxc12_pt,idxcpp12_pt,
     &     idx_f_temp,idx_h_temp,
     &     idxsop,idxspt,r12op,trunc_type
      integer, allocatable ::
     &     occ_def(:,:,:)

      type(operator), pointer::
     &     f_temp_pnt, h_temp_pnt, sop_pnt, spt_pnt

      integer, external ::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &     'Setting up CCPT-Lagrangian')
      end if

      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','r12op',ival=r12op)
      call get_argument_value('method.R12','trunc',ival=trunc_type)
      truncate = trunc_type.ge.0
      if (is_keyword_set('method.truncate').gt.0) then
        truncate = is_keyword_set('method.truncate').gt.0
        if(truncate)
     &     call get_argument_value('method.truncate','trunc_type',
     &                              ival=trunc_type)
      end if

      call atim_csw(cpu0,sys0,wall0)

      do ilabel = 1, nlabels
        idx = idx_oplist2(label_op(ilabel),op_info)
        if (idx.le.0)
     &       call quit(1,'set_ccpt_lagrangian',
     &       'label not on list: '//trim(label_op(ilabel)))
        if (ilabel.eq.1) idxlcc = idx
        if (ilabel.eq.2) idxham = idx
        if (ilabel.eq.3) idxtbar = idx
        if (ilabel.eq.4) idxtptbar = idx
        if (ilabel.eq.5) idxtop = idx
        if (ilabel.eq.6) idxtpt = idx
        if (ilabel.eq.7) idxr12 = idx
        if (ilabel.eq.8) idxc12 = idx
        if (ilabel.eq.9) idxc12_pt = idx
        if (ilabel.eq.10) idxcpp12 = idx
        if (ilabel.eq.11) idxcpp12_pt = idx
      end do

      call add_operator(op_f_temp,op_info)
      idx_f_temp = idx_oplist2(op_f_temp,op_info)
      f_temp_pnt => op_info%op_arr(idx_f_temp)%op

      if (ansatz.eq.0) then
c        ndef = 4
c        allocate(occ_def(ngastp,2,ndef))
c        occ_def(1:ngastp,1,1) = (/1,0,0,0/)
c        occ_def(1:ngastp,2,1) = (/1,0,0,0/)
c        occ_def(1:ngastp,1,2) = (/1,0,0,0/)
c        occ_def(1:ngastp,2,2) = (/0,1,0,0/)
c        occ_def(1:ngastp,1,3) = (/1,0,0,0/)
c        occ_def(1:ngastp,2,3) = (/0,1,0,0/)
c        occ_def(1:ngastp,1,4) = (/0,1,0,0/)
c        occ_def(1:ngastp,2,4) = (/0,1,0,0/)
c        call set_uop(f_temp_pnt,op_f_temp,.false.,
c     &       occ_def,ndef,orb_info)
c        deallocate(occ_def)
        if (mode(1:4).eq.'EXT ') then
          call set_hop(f_temp_pnt,op_f_temp,.false.,
     &         1,1,2,.true.,orb_info)
        else if (mode(4:4).eq.'0') then
          ndef = 3
          allocate(occ_def(ngastp,2,ndef))
          occ_def(1:ngastp,1,1) = (/1,0,0,0/)
          occ_def(1:ngastp,2,1) = (/1,0,0,0/)
          occ_def(1:ngastp,1,2) = (/0,1,0,0/)
          occ_def(1:ngastp,2,2) = (/0,1,0,0/)
          occ_def(1:ngastp,1,3) = (/0,0,0,1/)
          occ_def(1:ngastp,2,3) = (/0,0,0,1/)
          call set_uop(f_temp_pnt,op_f_temp,.false.,
     &         occ_def,ndef,orb_info)
          deallocate(occ_def)
        else
          call set_hop(f_temp_pnt,op_f_temp,.false.,
     &         1,1,2,.false.,orb_info)
        end if

        idxsop = idxtop
        idxspt = idxtpt
      else

        call set_hop(f_temp_pnt,op_f_temp,.false.,
     &       1,1,2,.true.,orb_info)

        ! Definition of the S = T+R operator.
        call add_operator(op_sop,op_info)
        idxsop = idx_oplist2(op_sop,op_info)
        sop_pnt => op_info%op_arr(idxsop)%op
        if (mode(5:8).eq.'NOR2') then ! no coupling to R2?
          call clone_operator(sop_pnt,op_info%op_arr(idxtop)%op,.false.,
     &         orb_info)
        else
          ! set R12 part:
          call set_r12gem(sop_pnt,op_sop,0,
     &         2 , 2 ,ansatz,orb_info)
          ! join with T
          call join_operator(sop_pnt,op_info%op_arr(idxtop)%op,orb_info)
        end if

        call init_formula(flist_t_r)

        call set_t_r(flist_t_r,.false.,idxsop,idxtop,
     &               idxr12,idxc12,idxcpp12,
     &               r12op,r12fix,op_info)

c        if (ntest.ge.1000) then
          call write_title(luout,wst_title,'T + CR')
          call print_form_list(luout,flist_t_r,op_info)
c        end if

        ! Definition of the S(pt) = T(pt)+[R,T] operator.
        call add_operator(op_spt,op_info)
        idxspt = idx_oplist2(op_spt,op_info)
        spt_pnt => op_info%op_arr(idxspt)%op

        if (r12op.gt.0) then
          ! set R12 part:
          call set_r12gem(spt_pnt,op_spt,0,
     &         3 , 3 ,ansatz,orb_info)
          ! join with T
          call join_operator(spt_pnt,op_info%op_arr(idxtpt)%op,orb_info)
        else
          call clone_operator(spt_pnt,op_info%op_arr(idxtpt)%op,.false.,
     &         orb_info)
        end if

        call init_formula(flist_t_r_pt)

        call set_t_r(flist_t_r_pt,.false.,idxspt,idxtpt,
     &               idxr12,idxc12_pt,idxcpp12_pt,
     &               r12op,r12fix,op_info)


c        if (ntest.ge.1000) then
          call write_title(luout,wst_title,'T + CR (PT)')
          call print_form_list(luout,flist_t_r_pt,op_info)
c        end if


      end if

      call add_operator(op_h_temp,op_info)
      idx_h_temp = idx_oplist2(op_h_temp,op_info)
      h_temp_pnt => op_info%op_arr(idx_h_temp)%op

      if (mode(4:4).eq.'0') then
        call set_hop(h_temp_pnt,op_h_temp,.false.,
     &       2,2,2,ansatz.ne.0,orb_info)
      else
        call set_hop(h_temp_pnt,op_h_temp,.false.,
     &       1,2,2,ansatz.ne.0,orb_info)
      end if

      ! initialize formula
      call init_formula(flist_lag)
      flist_pnt => flist_lag
      ! put [INIT] at the beginning
      call new_formula_item(flist_pnt,command_set_target_init,idxlcc)
      flist_pnt => flist_pnt%next
      ! set <0|T^+[H,T(pt)]|0>
      call expand_op_product2(flist_pnt,idxlcc,
     &     1d0,5,4,
     &     (/idxlcc,-idxsop,idx_h_temp,idxspt,idxlcc/),
     &     (/1     ,2      ,3     ,4     ,1/),
     &     -1,-1,
     &     (/3,4/),1,
     &     0,0,
     &     0,0,
     &     op_info)

      ! set <0|T(pt)^+[H,T]|0>
      ! advance pointer
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      call expand_op_product2(flist_pnt,idxlcc,
     &     1d0,5,4,
     &     (/idxlcc,-idxspt,idx_h_temp,idxsop,idxlcc/),
     &     (/1     ,2      ,3     ,4     ,1/),
     &     -1,-1,
     &     (/3,4/),1,
     &     0,0,
     &     0,0,
     &     op_info)

      ! set <0|T(pt)^+[F,T(pt)]|0>
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      call expand_op_product2(flist_pnt,idxlcc,
     &     1d0,5,4,
     &     (/idxlcc,-idxspt,idx_f_temp,idxspt,idxlcc/),
     &     (/1     ,2      ,3     ,4     ,1/),
     &     -1,-1,
     &     (/3,4/),1,
     &     0,0,
     &     0,0,
     &     op_info)

      if (ntest.ge.100) then
        call write_title(luout,wst_title,'raw formula')
        call print_form_list(luout,flist_lag,op_info)
      end if

      call form_op_replace(op_info%op_arr(idx_f_temp)%op%name,
     &                     op_info%op_arr(idxham)%op%name,.true.,
     &     flist_lag,op_info)
      call form_op_replace(op_info%op_arr(idx_h_temp)%op%name,
     &                     op_info%op_arr(idxham)%op%name,.true.,
     &     flist_lag,op_info)

      if (ansatz.gt.0) then
        call expand_subexpr(flist_lag,flist_t_r,.false.,op_info)
        call expand_subexpr(flist_lag,flist_t_r_pt,.false.,op_info)

        call sum_terms(flist_lag,op_info)

        call transpose_formula(flist_t_r,op_info)
        call transpose_formula(flist_t_r_pt,op_info)

        call expand_subexpr(flist_lag,flist_t_r,.false.,op_info)
        call expand_subexpr(flist_lag,flist_t_r_pt,.false.,op_info)

        call sum_terms(flist_lag,op_info)

        ! Produce truncated expansions.
        if (truncate)
     &     call r12_truncation(flist_lag,trunc_type,
     &     idxr12,idxham,idxtbar,idxtop,op_info)

c        call truncate_form(flist_lag,op_info)

        ! replace T12 -> T
        if (r12fix.and.r12op.gt.0) then
          if (r12op.ne.2) then
            call form_op_replace(op_info%op_arr(idxc12)%op%name,
     &                         op_info%op_arr(idxtop)%op%name,.true.,
     &           flist_lag,op_info)
            call form_op_replace(op_info%op_arr(idxc12_pt)%op%name,
     &                         op_info%op_arr(idxtop)%op%name,.true.,
     &           flist_lag,op_info)
          end if
          if (r12op.gt.1) then
            call form_op_replace(op_info%op_arr(idxcpp12)%op%name,
     &                       op_info%op_arr(idxtop)%op%name,
     &           flist_lag,op_info)
            call form_op_replace(op_info%op_arr(idxcpp12_pt)%op%name,
     &                       op_info%op_arr(idxtpt)%op%name,
     &           flist_lag,op_info)
          end if
        end if

      end if

      if (ntest.ge.100) then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,flist_lag,op_info)
      end if

      ! assign comment
      form_lag%comment = trim(title)
      ! write to disc
      write(name,'(a,".fml")') trim(form_lag%label)
      call file_init(form_lag%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_lag%fhand,flist_lag,
     &     form_lag%comment)

      call dealloc_formula_list(flist_lag)
      call dealloc_formula_list(flist_t_r)
      call dealloc_formula_list(flist_t_r_pt)

      call del_operator(op_f_temp,op_info)

      call atim_csw(cpu,sys,wall)
      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'CCPT Lagrangian',cpu-cpu0,sys-sys0,wall-wall0)

      end
