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
     &     idxham,idxtbar,idxtop,idxmet,idxrba,idxcbar,idxr12,idxc12,
     &     idxcpp12, idxcppbar,
     &     iblk_xxhp, iblk_pxhp, iblk_xxpp, iblk_pxpp,
     &     min_rank, max_rank, iprint,
     &     occ(ngastp,2)
 
      type(operator), pointer::
     &     sop_pnt, sbar_pnt, scr_pnt, scrbar_pnt

      integer, external::
     &     idx_oplist2, max_rank_op, iblk_occ

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

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
        if (ilabel.eq.1) idxmet = idx
        if (ilabel.eq.2) idxr12 = idx
        if (ilabel.eq.3) idxtbar = idx
        if (ilabel.eq.4) idxtop = idx
      end do

      ! initialize formula
      call init_formula(flist_exp)
      flist_pnt => flist_exp
      ! put [INIT] at the beginning
      call new_formula_item(flist_pnt,command_set_target_init,idxmet)
      flist_pnt => flist_pnt%next

      ! formulae are set up e.g. by ...
      ! expand <0|Sbar S|0>
c      call expand_op_product2(flist_pnt,idxmet,
c     &     1d0,4,3,
c     &     (/idxmet,idx_sbar,idx_sop,idxmet/),
c     &     (/1     ,2       ,3      ,1/),
c     &     -1,-1,
c     &     0,0,
c     &     0,0,
c     &     0,0,
c     &     op_info)

      if (ntest.ge.100) then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,flist_exp,op_info)
      end if

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
