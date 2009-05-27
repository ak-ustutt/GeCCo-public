*----------------------------------------------------------------------*
      subroutine set_r12_lagrangian(form_cclag,
     &     title,label,nlabels,ansatz,
     &     op_info,orb_info,form_info)
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
      include 'mdef_formula_info.h'
      include 'def_formula_item.h'
      include 'ifc_input.h'

      type(formula), intent(inout), target ::
     &     form_cclag
      type(formula_info), intent(inout) ::
     &     form_info
      integer, intent(in) ::
     &     ansatz, nlabels
      character(*), intent(in) ::
     &     label(nlabels), title

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf), intent(inout) ::
     &     orb_info


      ! local constants
      character(3), parameter ::
     &     op_sop    = '_S_',
     &     op_sba    = '_SB',
     &     op_scr    = '_T_',
     &     op_scrbar = '_TB'

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist_lag, flist_t_cr, flist_tbar_cbarr
      type(formula_item), pointer ::
     &     flist_pnt, fl_t_cr_pnt

      integer ::
     &     nterms, idx_sop, idx_sbar, ndef, idxrint, ilabel, idx,
     &     idx_scr,idx_scrbar, r12op,
     &     idxham,idxtbar,idxtop,idxlcc,idxrba,idxcbar,idxr12,idxc12,
     &     idxcpp12, idxcppbar,
     &     iblk_xxhp, iblk_pxhp, iblk_xxpp, iblk_pxpp,
     &     min_rank, max_rank, iprint, trunc_type, trunc_t1x, h0_t1x,
     &     occ(ngastp,2)
      logical ::
     &     r12fix, truncate, l_h0d
      integer ::
     &     extend

      type(operator), pointer::
     &     sop_pnt, sbar_pnt, scr_pnt, scrbar_pnt

      integer, external::
     &     idx_oplist2, max_rank_op, iblk_occ

      character(len=8) ::
     &     trmode

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(luout,*) '==============================='
        write(luout,*) ' output from set_r12_lagrangian'
        write(luout,*) '==============================='
        write(luout,*) ' ansatz = ',ansatz
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! Are we fixing the F12 amplitudes?
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','extend',ival=extend)
      call get_argument_value('method.R12','r12op',ival=r12op)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.R12','trunc',ival=trunc_type)
      truncate = trunc_type.ge.0
      if (is_keyword_set('method.truncate').gt.0) then
        truncate = is_keyword_set('method.truncate').gt.0
        if(truncate)
     &     call get_argument_value('method.truncate','trunc_type',
     &                              ival=trunc_type)
      end if
      call get_argument_value('method.R12','T1ext',ival=trunc_t1x)
      call get_argument_value('method.R12','H0_T1ext',ival=h0_t1x)
      call get_argument_value('method.R12','H0d',lval=l_h0d)
      
      if (extend.gt.0) call quit(1,'set_r12_lagrangian',
     &     'do not use "extend" for CC (use "r12op" instead)!')

c      r12fix = r12fix .or. extend.gt.0

      ! get indices
c      if (nlabels.ne.8.and..not.r12fix) then
c        write(luout,*) 'nlabels = ',nlabels
c        call quit(1,'set_r12_lagrangian',
c     &     'I expect exactly 8 labels')
c      end if
c      if (nlabels.lt.6.and.r12fix) then
c        write(luout,*) 'nlabels = ',nlabels
c        call quit(1,'set_mp2_r12_lagrangian fixed amp.',
c     &     'I expect > 6 labels')c
c      end if

      idxcbar = -1
      idxc12  = -1
      idxcppbar = -1
      idxcpp12  = -1
      do ilabel = 1, nlabels
        idx = idx_oplist2(label(ilabel),op_info)
        if (idx.le.0)
     &       call quit(1,'set_mp2_r12_lagrangian',
     &       'label not on list: '//trim(label(ilabel)))
        if (ilabel.eq.1) idxlcc = idx
        if (ilabel.eq.2) idxham = idx
        if (ilabel.eq.3) idxr12 = idx
        if (ilabel.eq.4) idxrba = idx
        if (ilabel.eq.5) idxtbar = idx
        if (ilabel.eq.6) idxtop = idx
        if (r12op.ne.2) then
          if (ilabel.eq.7) idxcbar = idx
          if (ilabel.eq.8) idxc12 = idx
        else
          if (ilabel.eq.7) idxcppbar = idx
          if (ilabel.eq.8) idxcpp12 = idx
        end if
        if (r12op.ne.2) then
          if (ilabel.eq.9) idxcppbar = idx
          if (ilabel.eq.10) idxcpp12 = idx
        end if
      end do

      ! Definition of the S=T+CR operator.
      call add_operator(op_sop,op_info)
      idx_sop = idx_oplist2(op_sop,op_info)
      sop_pnt => op_info%op_arr(idx_sop)%op

      min_rank = 2

      ! set CR part:
      call set_r12gem(sop_pnt,op_sop,0,
     &     min_rank,max_rank,ansatz,orb_info)

      ! join with T
      call join_operator(sop_pnt,op_info%op_arr(idxtop)%op,orb_info)

      ! define Sbar for the projection
      call add_operator(op_sba,op_info)
      idx_sbar = idx_oplist2(op_sba,op_info)
      sbar_pnt => op_info%op_arr(idx_sbar)%op

      call clone_operator(sbar_pnt,sop_pnt,.true.,orb_info)
c      sbar_pnt%dagger = .true.

      ! combine the C and R operators (S=T+CR).
      call init_formula(flist_t_cr)
      call set_t_r(flist_t_cr,.false.,.false.,
     &     idx_sop,idxtop,
     &     idxr12,-1,idxc12,idxcpp12,
     &     r12op,r12fix,op_info)

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'T + CR')
        call print_form_list(luout,flist_t_cr,op_info)
      end if

      ! Must also form SBAR.
      call init_formula(flist_tbar_cbarr)
      call set_t_r(flist_tbar_cbarr,.true.,.false.,
     &     idx_sbar,idxtbar,
     &     idxr12,-1,idxcbar,idxcppbar,
     &     r12op,r12fix,op_info)

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'TBAR + R CBAR')
        call print_form_list(luout,flist_tbar_cbarr,op_info)
      end if

      ! and now: the actual formula
      ! initialize formula
      call init_formula(flist_lag)
      flist_pnt => flist_lag
      ! put [INIT] at the beginning
      call new_formula_item(flist_pnt,command_set_target_init,idxlcc)
      flist_pnt => flist_pnt%next

      ! expand <0|(1+Sbar) e^{-S} H e^S|0> =
      ! <0| e^{-S} H e^S|0> +
      call expand_op_bch(flist_pnt,2,idxlcc,
     &     1d0,-1,idxham,1d0,idx_sop,1,-1,op_info)

      ! advance pointer
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! <0|Sbar e^{-S} H e^S|0>
      call expand_op_bch(flist_pnt,4,idxlcc,
     &     1d0,idx_sbar,idxham,1d0,idx_sop,1,-1,op_info)

      ! insert here procedure to produce approx. expansions      
      ! ...
      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'raw formula')
        call print_form_list(luout,flist_lag,op_info)
      end if

      ! replace S by T+CR
      call expand_subexpr(flist_lag,flist_t_cr,.false.,op_info)

      ! sum up duplicate terms (due to S->T+CR replacement)
      call sum_terms(flist_lag,op_info)

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'after replacing S')
        call print_form_list(luout,flist_lag,op_info)
      end if

      ! replace Sbar by Tbar + R^t CBAR
      call expand_subexpr(flist_lag,flist_tbar_cbarr,.false.,op_info)

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'after replacing S')
        call print_form_list(luout,flist_lag,op_info)
      end if

      ! Produce truncated expansions if required.
      if (truncate)
     &     call r12_truncation(flist_lag,trunc_type,h0_t1x,
     &     idxr12,idxham,idxtbar,idxtop,idxcbar,idxc12,op_info)
      if (trunc_t1x.gt.0.and.h0_t1x.ne.-1) then
        trmode = '        '
        write(trmode,'("ord",i1," ",i1)') trunc_t1x, h0_t1x
        if (l_h0d) trmode(7:7) = 'd'
        call t1x_r12_truncation(flist_lag,trmode,
     &       idxr12,idxtbar,idxham,idxtop,op_info)
      end if

      ! sum up duplicate terms (due to S->T+CR replacement)
      call sum_terms(flist_lag,op_info)

      ! replace T12 -> T
      if (r12fix.and.r12op.gt.0) then
        if (r12op.ne.2) then
          call form_op_replace(op_info%op_arr(idxc12)%op%name,
     &                       op_info%op_arr(idxtop)%op%name,.true.,
     &     flist_lag,op_info)
          call form_op_replace(op_info%op_arr(idxcbar)%op%name,
     &                       op_info%op_arr(idxtbar)%op%name,.true.,
     &     flist_lag,op_info)
        end if
        if (r12op.gt.1) then
          call form_op_replace(op_info%op_arr(idxcpp12)%op%name,
     &                       op_info%op_arr(idxtop)%op%name,.true.,
     &     flist_lag,op_info)
          call form_op_replace(op_info%op_arr(idxcppbar)%op%name,
     &                       op_info%op_arr(idxtbar)%op%name,.true.,
     &     flist_lag,op_info)
        end if
cc        call form_op_replace(op_scr,op_info%op_arr(idxtop)%op%name,
cc     &     flist_lag,op_info)
      end if

      ! post_processing and term counting:
      if(.not.r12fix.and.r12op.le.1)then
        iprint = iprlvl
        call r12_form_post(flist_lag,nterms,
     &       idxtbar,idxcbar,idxham,idxtop,idxc12, iprint,
     &       op_info)
      else if (r12fix.and.r12op.le.1) then
        iprint = iprlvl
        call cc_form_post(flist_lag,nterms,
     &       idxtbar,idxham,idxtop, iprint,
     &       op_info)
      endif

      if (ntest.ge.100) then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,flist_lag,op_info)
      end if

      ! assign comment
      form_cclag%comment = trim(title)
      ! write to disc
      write(name,'(a,".fml")') trim(form_cclag%label)
      call file_init(form_cclag%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_cclag%fhand,flist_lag,
     &     form_cclag%comment)

      call dealloc_formula_list(flist_t_cr)
      call dealloc_formula_list(flist_tbar_cbarr)
      call dealloc_formula_list(flist_lag)

      ! remove the formal operators
      call del_operator(op_sba,op_info)
      call del_operator(op_sop,op_info)

      call atim_csw(cpu,sys,wall)
      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'CC-R12 Lagrangian',cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
