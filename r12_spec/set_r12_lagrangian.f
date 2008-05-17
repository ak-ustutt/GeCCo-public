*----------------------------------------------------------------------*
      subroutine set_r12_lagrangian(form_cclag,
     &     title,label,nlabels,ansatz,
     &     op_info,orb_info)
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
     &     ntest = 1000

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
     &     form_cclag
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
     &     op_scr    = '_T_'

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist_lag, flist_t_cr, flist_tbar_cbarr
      type(formula_item), pointer ::
     &     flist_pnt, fl_t_cr_pnt

      integer ::
     &     nterms, idx_sop, idx_sbar, ndef, idxrint, ilabel, idx,
     &     idx_scr,
     &     idxham,idxtbar,idxtop,idxlcc,idxrba,idxcbar,idxr12,idxc12,
     &     min_rank, max_rank, iprint
      logical ::
     &     r12fix
      integer ::
     &     extend

      type(operator), pointer::
     &     sop_pnt, sbar_pnt, scr_pnt

      integer, external::
     &     idx_oplist2, max_rank_op

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
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      
      r12fix = r12fix .or. extend.gt.0

      ! get indices
      if (nlabels.ne.8.and..not.r12fix) then
        write(luout,*) 'nlabels = ',nlabels
        call quit(1,'set_r12_lagrangian',
     &     'I expect exactly 8 labels')
      end if
      if (nlabels.ne.6.and.r12fix) then
        write(luout,*) 'nlabels = ',nlabels
        call quit(1,'set_mp2_r12_lagrangian fixed amp.',
     &     'I expect exactly 6 labels')
      end if

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
        if (ilabel.eq.7) idxcbar = idx
        if (ilabel.eq.8) idxc12 = idx
      end do

      ! Definition of the S=T+CR operator.
      call add_operator(op_sop,op_info)
      idx_sop = idx_oplist2(op_sop,op_info)
      sop_pnt => op_info%op_arr(idx_sop)%op

      min_rank = 2
c      max_rank = 2
c      if(.not.r12fix)then
c        max_rank = max_rank_op('A',op_info%op_arr(idxc12)%op,.false.)
c      endif

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
      fl_t_cr_pnt => flist_t_cr
      call new_formula_item(fl_t_cr_pnt,command_set_target_init,idx_sop)
      fl_t_cr_pnt => fl_t_cr_pnt%next
      call expand_op_product(fl_t_cr_pnt,idx_sop,
     &     1d0,1,idxtop,-1,-1,
     &     0,0,.false.,op_info)
      do while(associated(fl_t_cr_pnt%next))
        fl_t_cr_pnt => fl_t_cr_pnt%next
      end do

      ! in order to use different names for T connected with geminal:
      if (extend.gt.0) then
        call add_operator(op_scr,op_info)
        idx_scr = idx_oplist2(op_scr,op_info)
        scr_pnt => op_info%op_arr(idx_scr)%op
        call clone_operator(scr_pnt,op_info%op_arr(idxtop)%op,
     &       .false.,orb_info)
      end if

      ! Form of the R12 part depends on whether the amplitudes are fixed.
      if(.not.r12fix)then
        call expand_op_product(fl_t_cr_pnt,idx_sop,
     &       1d0,2,(/idxc12,idxr12/),-1,-1,
     &       (/1,2/),1,.false.,op_info)
      else
        call expand_op_product(fl_t_cr_pnt,idx_sop,
     &       1d0,1,idxr12,-1,-1,
     &       0,0,.false.,op_info)
      endif

      if(extend.gt.0)then
        do while(associated(fl_t_cr_pnt%next))
          fl_t_cr_pnt => fl_t_cr_pnt%next
        enddo

        call expand_op_product2(fl_t_cr_pnt,idx_sop,
     &       1d0,4,3,
     &       (/idx_sop,idxr12,idx_scr,idx_sop/),
     &       (/1      ,2     ,3       ,1     /),
     &       (/1,1,1,1/),-1,
     &       (/2,3/),1,
     &       0,0,
     &       0,0,
     &       op_info)
      endif

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'T + CR')
        call print_form_list(luout,flist_t_cr,op_info)
      end if

      ! Must also form SBAR.
      call init_formula(flist_tbar_cbarr)
      fl_t_cr_pnt => flist_tbar_cbarr
      call new_formula_item(fl_t_cr_pnt,
     &                      command_set_target_init,idx_sbar)
      fl_t_cr_pnt => fl_t_cr_pnt%next
      call expand_op_product(fl_t_cr_pnt,idx_sbar,
     &     1d0,1,idxtbar,-1,-1,
     &     0,0,.false.,op_info)
      do while(associated(fl_t_cr_pnt%next))
        fl_t_cr_pnt => fl_t_cr_pnt%next
      end do

      if(.not.r12fix)then
        call expand_op_product2(fl_t_cr_pnt,idx_sbar,
     &       1d0,4,3,
     &       (/idx_sbar,-idxr12,idxcbar,idx_sbar/),(/1,2,3,1/),
     &       -1,-1,
     &       (/2,3/),1,
     &       0,0,
     &       0,0,
     &       op_info)
      else
        call expand_op_product2(fl_t_cr_pnt,idx_sbar,
     &       1d0,3,2,
     &       (/idx_sbar,-idxr12,idx_sbar/),(/1,2,1/),
     &       -1,-1,
     &       0,0,
     &       0,0,
     &       0,0,
     &       op_info)
      endif
c      call expand_op_product(fl_t_cr_pnt,idx_sbar,
c     &     1d0,2,(/idxrba,idxcbar/),-1,-1,
c     &     (/1,2/),1,.false.,op_info)

      if(extend.gt.0)then
        do while(associated(fl_t_cr_pnt%next))
          fl_t_cr_pnt => fl_t_cr_pnt%next
        enddo

        call expand_op_product2(fl_t_cr_pnt,idx_sbar,
     &       1d0,4,3,
     &       (/idx_sbar,idxtbar,-idxr12,idx_sbar/),(/1,2,3,1/),
     &       (/1,1,1,1/),-1,
     &       (/2,3/),1,
     &       0,0,
     &       0,0,
     &       op_info)
      endif

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

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'after replacing S')
        call print_form_list(luout,flist_lag,op_info)
      end if

      call sum_terms(flist_lag,op_info)

      ! replace Sbar by Tbar + R^t CBAR
      call expand_subexpr(flist_lag,flist_tbar_cbarr,.false.,op_info)

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'after replacing S')
        call print_form_list(luout,flist_lag,op_info)
      end if

      ! Produce truncated expansions.
      call truncate_form(flist_lag,op_info)

      ! rename _T_ -> T
      if (extend.gt.0) then
        call form_op_replace(op_scr,op_info%op_arr(idxtop)%op%name,
     &     flist_lag,op_info)
      end if

      ! sum up duplicate terms (due to S->T+CR replacement)
      call sum_terms(flist_lag,op_info)

      ! post_processing and term counting:
      if(.not.r12fix)then
        iprint = iprlvl
        call r12_form_post(flist_lag,nterms,
     &       idxtbar,idxcbar,idxham,idxtop,idxc12, iprint,
     &       op_info)
      endif

c      if (ntest.ge.100) then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,flist_lag,op_info)
c      end if

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
      if (extend.gt.0) call del_operator(op_scr,op_info)

      call atim_csw(cpu,sys,wall)
      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'CC-R12 Lagrangian',cpu-cpu0,sys-sys0,wall-wall0)

c dbg
c      stop 'testing'
c dbg
      return
      end
