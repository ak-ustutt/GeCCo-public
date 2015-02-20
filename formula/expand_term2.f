*----------------------------------------------------------------------*
      subroutine expand_term2(fl_expand,nterms,
     &                        njoined,f_term,fpl_intm,force,op_info)
*----------------------------------------------------------------------*
*     expand O1.O2...Int...On (on f_term as formula list) to
*            O1.O2...I1aI1b...On + O1.O2...I2...On + ...
*     where the definition Int = I1aI1b... + I2... + ...
*     is given as formula pointer list
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula_item_array.h'
      include 'def_formula_item_list.h'

      integer, parameter ::
     &     ntest = 00
      
      type(formula_item), target, intent(out) ::
     &     fl_expand
      integer, intent(out) ::
     &     nterms
      type(formula_item), target, intent(in) ::
     &     f_term
      type(formula_item_list), target, intent(in) ::
     &     fpl_intm
      type(operator_info) ::
     &     op_info
      logical, intent(in)::
     &     force
      integer, intent(in) ::
     &     njoined

      type(contraction) ::
     &     proto
      logical ::
     &     adj_intm
      integer ::
     &     nvtx, narc, narc0, ivtx, jvtx, kvtx, iarc,
     &     iop_intm, iblk_intm, iblk, iadd
      integer ::
     &     occ_temp(ngastp,2)
      type(contraction), pointer ::
     &     term, intm
      type(contraction) ::
     &     intm_spl, term_rem
      type(formula_item_list), pointer ::
     &     fpl_intm_pnt
      type(formula_item), pointer ::
     &     fl_expand_pnt
      type(operator), pointer ::
     &     op_intm
      integer, pointer ::
     &     ivtx_term_reo(:), ivtx_intm_reo(:),
     &     occ_vtx(:,:,:), svmap(:), vtxmap(:),
     &     ipos_vtx(:)
      logical, pointer ::
     &     fix_vtx(:)

      integer, external ::
     &     vtx_in_contr, ifac, idxlist, ielsqsum
      call quit(1,'expand_term2','call to obsolete routine')

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'Here speaks expand_term2')
      end if

      iop_intm  = fpl_intm%item%contr%idx_res
      iblk_intm = fpl_intm%item%contr%iblk_res
      adj_intm  = fpl_intm%item%contr%dagger
      op_intm => op_info%op_arr(iop_intm)%op
c      njoined = op_intm%njoined

      allocate(ipos_vtx(njoined))

      term => f_term%contr
      call get_vtx_in_contr(ipos_vtx,iop_intm,adj_intm,njoined,1,term)

      iblk = term%vertex(ipos_vtx(1))%iblk_op
      if (njoined.gt.1)
     &     iblk = (iblk-1)/njoined + 1

      if (iblk.ne.iblk_intm)
     &     call quit(1,'expand_term2','inconsistency')

      call init_contr(intm_spl)
      call set_primitive_contr2(intm_spl,
     &     1d0,iop_intm,iblk_intm,
     &     iop_intm,iblk_intm,.false.,
     &     op_info)
      
      call init_contr(term_rem)
      call split_contr2(.false.,term_rem,intm_spl,term,op_info)

      call dealloc_contr(intm_spl)
      stop 'testing'

      nterms = 0
      fpl_intm_pnt => fpl_intm
      fl_expand_pnt => fl_expand
      ! loop over intermediate blocks
      do

        intm => fpl_intm_pnt%item%contr

c        call join_contr2(term_rem,intm,op_info)

        if (.not.associated(fpl_intm_pnt%next)) exit
        fpl_intm_pnt => fpl_intm_pnt%next
      end do

      call dealloc_contr(term_rem)

      call dealloc_contr(proto)

      deallocate(ipos_vtx)

      return
      end
