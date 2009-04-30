*----------------------------------------------------------------------*
      subroutine set_r12_metric(form_metric,
     &     title,label,nlabels,ansatz,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*     set up 
*          S_0 = <0|(1+Tbar)e^{-T}Te^{T}|0> = <0|Tbar T|0>
*     
*     which defines the metric of the CC approach needed for
*     response equations. Not necessary for normal CC which uses
*     a unit metric, but for R12 a modified metric occurs due to
*     non-unit R^+R integrals (X-intermediate).
*
*     andreas, july 2008
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
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'ifc_input.h'

      type(formula), intent(inout), target ::
     &     form_metric
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
     &     flist_metric, flist_t_cr, flist_tbar_cbarr
      type(formula_item), pointer ::
     &     flist_pnt, fl_t_cr_pnt

      integer ::
     &     nterms, idx_sop, idx_sbar, ndef, idxrint, ilabel, idx,
     &     idx_scr,idx_scrbar, r12op,
     &     idxham,idxtbar,idxtop,idxmet,idxrba,idxcbar,idxr12,idxc12,
     &     idxcpp12, idxcppbar,
     &     iblk_xxhp, iblk_pxhp, iblk_xxpp, iblk_pxpp,
     &     min_rank, max_rank, iprint,
     &     occ(ngastp,2)
      logical ::
     &     r12fix, opt
      integer ::
     &     extend

      type(operator), pointer::
     &     sop_pnt, sbar_pnt, scr_pnt, scrbar_pnt

      integer, external::
     &     idx_oplist2, max_rank_op, iblk_occ

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.eq.100) then
        write(luout,*) '==============================='
        write(luout,*) ' output from set_r12_metricrangian'
        write(luout,*) '==============================='
        write(luout,*) ' ansatz = ',ansatz
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! Are we fixing the F12 amplitudes?
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','extend',ival=extend)
      call get_argument_value('method.R12','r12op',ival=r12op)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.R12','opt',lval=opt)
      
      if (extend.gt.0) call quit(1,'set_r12_metric',
     &     'do not use "extend" for CC (use "r12op" instead)!')

      do ilabel = 1, nlabels
        idx = idx_oplist2(label(ilabel),op_info)
        if (idx.le.0)
     &       call quit(1,'set_r12_metric',
     &       'label not on list: '//trim(label(ilabel)))
        if (ilabel.eq.1) idxmet = idx
        if (ilabel.eq.2) idxr12 = idx
        if (ilabel.eq.3) idxtbar = idx
        if (ilabel.eq.4) idxtop = idx
        if (r12op.ne.2) then
          if (ilabel.eq.5) idxcbar = idx
          if (ilabel.eq.6) idxc12 = idx
        else
          if (ilabel.eq.5) idxcppbar = idx
          if (ilabel.eq.6) idxcpp12 = idx
        end if
        if (r12op.ne.2) then
          if (ilabel.eq.7) idxcppbar = idx
          if (ilabel.eq.8) idxcpp12 = idx
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
      call init_formula(flist_metric)
      flist_pnt => flist_metric
      ! put [INIT] at the beginning
      call new_formula_item(flist_pnt,command_set_target_init,idxmet)
      flist_pnt => flist_pnt%next

      ! expand <0|Sbar S|0>
      call expand_op_product2(flist_pnt,idxmet,
     &     1d0,4,3,
     &     (/idxmet,idx_sbar,idx_sop,idxmet/),
     &     (/1     ,2       ,3      ,1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     0,0,
     &     op_info)

      ! replace S by T+CR
      call expand_subexpr(flist_metric,flist_t_cr,.false.,op_info)

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'after replacing S')
        call print_form_list(luout,flist_metric,op_info)
      end if

      ! replace Sbar by Tbar + R^t CBAR
      call expand_subexpr(flist_metric,flist_tbar_cbarr,.false.,op_info)

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'after replacing SBAR')
        call print_form_list(luout,flist_metric,op_info)
      end if

      ! delete redundant operator blocks (if more than one block version)
      if (opt) then
        call r12_opt_truncation(flist_metric,idxtop,idxc12,op_info)
        call r12_opt_truncation(flist_metric,idxtbar,idxcbar,op_info)
      end if

      ! sum up duplicate terms (due to S->T+CR replacement)
      call sum_terms(flist_metric,op_info)

      ! replace T12 -> T
      if (r12fix.and.r12op.gt.0) then
        if (r12op.ne.2) then
          call form_op_replace(op_info%op_arr(idxc12)%op%name,
     &                       op_info%op_arr(idxtop)%op%name,.true.,
     &     flist_metric,op_info)
          call form_op_replace(op_info%op_arr(idxcbar)%op%name,
     &                       op_info%op_arr(idxtbar)%op%name,.true.,
     &     flist_metric,op_info)
        end if
        if (r12op.gt.1) then
          call form_op_replace(op_info%op_arr(idxcpp12)%op%name,
     &                       op_info%op_arr(idxtop)%op%name,.true.,
     &     flist_metric,op_info)
          call form_op_replace(op_info%op_arr(idxcppbar)%op%name,
     &                       op_info%op_arr(idxtbar)%op%name,.true.,
     &     flist_metric,op_info)
        end if
      end if

      if (ntest.ge.100) then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,flist_metric,op_info)
      end if

      ! assign comment
      form_metric%comment = trim(title)
      ! write to disc
      write(name,'(a,".fml")') trim(form_metric%label)
      call file_init(form_metric%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_metric%fhand,flist_metric,
     &     form_metric%comment)

      call dealloc_formula_list(flist_t_cr)
      call dealloc_formula_list(flist_tbar_cbarr)
      call dealloc_formula_list(flist_metric)

      ! remove the formal operators
      call del_operator(op_sba,op_info)
      call del_operator(op_sop,op_info)

      call atim_csw(cpu,sys,wall)
c      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'CC-R12 metric',cpu-cpu0,sys-sys0,wall-wall0)

c dbg
c      if (r12op.gt.0) stop 'testing'
c dbg
      return
      end
