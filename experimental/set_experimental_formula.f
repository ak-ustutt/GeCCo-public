*----------------------------------------------------------------------*
      subroutine set_experimental_formula(form_exp,
     &     title,label,nlabels,ansatz,approx,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     template for setting up formulae
*     ansatz, approx can be freely used
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
     &     form_exp
      integer, intent(in) ::
     &     ansatz, nlabels
      character(*), intent(in) ::
     &     label(nlabels), title, approx

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf), intent(inout) ::
     &     orb_info


      ! local constants
c      character(3), parameter ::
c     &     op_sop    = '_S_',
c     &     op_sba    = '_SB',
c     &     op_scr    = '_T_',
c     &     op_scrbar = '_TB'

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist_exp, flist_t_cr, flist_tbar_cbarr
      type(formula_item), pointer ::
     &     flist_pnt, fl_t_cr_pnt

      integer ::
     &     nterms, idx_sop, idx_sbar, ndef, idxrint, ilabel, idx,
     &     idx_scr,idx_scrbar, r12op,
     &     idxham,idx_l,idx_t,idx_resplag,idxrba,idxcbar,idx_op,
     &     idxc12,idxcpp12,idxcppbar,
     &     iblk_xxhp, iblk_pxhp, iblk_xxpp, iblk_pxpp,
     &     min_rank, max_rank, iprint, fact, iterm,
     &     occ(ngastp,2)
 
      type(operator), pointer::
     &     sop_pnt, sbar_pnt, scr_pnt, scrbar_pnt

      integer, external::
     &     idx_oplist2, max_rank_op, iblk_occ

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      integer, parameter ::
     &     maxterm = 4

      integer ::
     &     op_idx(5+maxterm), op_num(5+maxterm),
     &     connect(2*maxterm+1), switch

      integer, external ::
     &      factorial

      if (ntest.eq.100) then
        call write_title(luout,wst_dbg_subr,
     &       'output from set_experimental')
        write(luout,*) ' ansatz = ',ansatz
        write(luout,*) ' approx = ',trim(approx)
      end if

      call atim_csw(cpu0,sys0,wall0)
      ! transform labels into indices in this way
      do ilabel = 1, nlabels
        idx = idx_oplist2(label(ilabel),op_info)
        if (idx.le.0)
     &       call quit(1,'set_exp',
     &       'label not on list: '//trim(label(ilabel)))
        if (ilabel.eq.1) idx_resplag = idx
        if (ilabel.eq.2) idx_op = idx
        if (ilabel.eq.3) idx_t = idx
      end do

      ! initialize formula
      call init_formula(flist_exp)
      flist_pnt => flist_exp
      ! put [INIT] at the beginning
      call new_formula_item(flist_pnt,
     &    command_set_target_init,idx_resplag)
c      flist_pnt => flist_pnt%next

      ! expand <0|Tbar exp(-T) Op exp(T) |0>
      op_idx(1) = idx_resplag
      op_idx(2) = -idx_t
      op_idx(3) = idx_op
      op_idx(4) = idx_resplag
      op_num(1:4) = (/1,2,3,1/)
      connect(1) = 0
      switch = 1
      do iterm = 1,maxterm+1
        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        end do
        fact = factorial(iterm-1)
   
        call expand_op_product2(flist_pnt,idx_resplag,
     &       1d0/fact,iterm+3,iterm+2,
     &       op_idx(1:iterm+3),
     &       op_num(1:iterm+3),
     &       -1,-1,
     &       connect(1:2*iterm-switch),iterm-1,
     &       0,0,
     &       0,0,
     &       op_info)
        op_idx(iterm+3) = idx_t
        op_idx(iterm+4) = idx_resplag
        op_num(iterm+3) = iterm+3
        op_num(iterm+4) = 1
        connect(2*iterm-1) = 3
        connect(2*iterm) = iterm+3
        switch = 2
      end do

c      ! formulae are set up e.g. by ...
c      ! expand <0|T2bar [F,T2]|0>
c      call expand_op_product2(flist_pnt,idx_mp2lag,
c     &     1d0,5,4,
c     &     (/idx_mp2lag,-idx_t2,idx_fock,idx_t2,idx_mp2lag/),
c     &     (/1         ,2      ,3       ,4     ,1/),
c     &     -1,-1,
c     &     (/3,4/),1,
c     &     0,0,
c     &     0,0,
c     &     op_info)

c      do while(associated(flist_pnt%next))
c        flist_pnt => flist_pnt%next
c      end do
c      ! expand <0|PHI T2|0>
c      call expand_op_product2(flist_pnt,idx_mp2lag,
c     &     1d0,4,3,
c     &     (/idx_mp2lag,idx_phi,idx_t2,idx_mp2lag/),
c     &     (/1         ,2      ,3     ,1         /),
c     &     -1,-1,
c     &     (/2,3/),1,
c     &     0,0,
c     &     0,0,
c     &     op_info)

c      do while(associated(flist_pnt%next))
c        flist_pnt => flist_pnt%next
c      end do
c      ! expand <0|T2bar PHI|0>
c      call expand_op_product2(flist_pnt,idx_mp2lag,
c     &     1d0,4,3,
c     &     (/idx_mp2lag,-idx_t2,idx_phi,idx_mp2lag/),
c     &     (/1         ,2      ,3     ,1         /),
c     &     -1,-1,
c     &     0,0,
c     &     0,0,
c     &     0,0,
c     &     op_info)

c      if (ntest.ge.100) then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,flist_exp,op_info)
c      end if

      ! assign comment
      form_exp%comment = trim(title)
      ! write to disc
      write(name,'(a,".fml")') trim(form_exp%label)
      call file_init(form_exp%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_exp%fhand,flist_exp,
     &     form_exp%comment)

      call dealloc_formula_list(flist_exp)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'set_exp',cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
