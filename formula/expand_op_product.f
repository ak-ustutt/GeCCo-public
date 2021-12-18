*----------------------------------------------------------------------*
      subroutine expand_op_product(form_list,idx_res,
     &                             fac,nops,idx_op,
     &                             iblk_min_in, iblk_max_in,
     &                             connect,nconnect,force,
     &                             op_info)
*----------------------------------------------------------------------*
*     given a list of operator indices and a result operator generate 
*     all contractions arising from
*
*        fac * <C_res|Op(1)Op(2) ....Op(n)|A_res>
*
*     idx_res: operator describing resulting "shape" of contraction
*     idx_op(nops): operators to be contracted
*     connect(2,nconnect): list of operator pairs that must be connected
*                          e.g. (1,3) if first and third operator should
*                          be connected.
*     force is set true if we want the connection of the first and last 
*     vertices in a non-standard way. If this is so then the first 
*     element of connect should be (1,nconnect).
*
*     andreas, june 2007
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'ifc_operators.h'
      include 'ifc_baserout.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(inout), target ::
     &     form_list
      real(8), intent(in) ::
     &     fac
      integer, intent(in) ::
     &     nops, nconnect,
     &     idx_res, idx_op(nops), iblk_min_in(nops), iblk_max_in(nops),
     &     connect(2,nconnect)
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in)::
     &     force

      type(formula_item), pointer ::
     &     form_pnt
      type(contraction) ::
     &     proto
      type(operator), pointer ::
     &     op_res, op
      logical ::
     &     fix_vtx(nops)
      logical ::
     &     ok
      integer ::
     &     nvtx, narc, iarc, iblk_res, iop, jop
      integer ::
     &     iblk_min(nops), iblk_max(nops), iblk_op(nops),
     &     occ_test(ngastp,2), occ_temp(ngastp,2),
     &     idx_op_vtx(nops+2), idx_sv_vtx(nops+2),
     &     iblk_min2(nops+2), iblk_max2(nops+2),
     &     connect2(2,nconnect)
      integer, pointer ::
     &     occ_vtx(:,:,:), neqv(:), idx_eqv(:,:), iop_typ(:)

      integer, external ::
     &     vtx_type

      ! passing arguments to expand_op_product2
      if (force) call quit(1,'obsolete "force" feature?')
      idx_op_vtx(2:nops+1) = idx_op(1:nops)
      idx_op_vtx(1) = idx_res
      idx_op_vtx(nops+2) = idx_res
      do iop = 1,nops+1
        idx_sv_vtx(iop) = iop
      end do
      idx_sv_vtx(nops+2) = 1
      iblk_min2(2:nops+1) = iblk_min_in(1:nops)
      iblk_max2(2:nops+1) = iblk_max_in(1:nops)
      if (iblk_min2(2).lt.0) then
        iblk_min2(1) = -1
      else
        iblk_min2(1) = 1
        iblk_min2(nops+2) = 1
      end if
      if (iblk_max2(2).lt.0) then
        iblk_max2(1) = -1
      else
        iblk_max2(1) = op_info%op_arr(idx_res)%op%n_occ_cls
        iblk_max2(nops+2) = op_info%op_arr(idx_res)%op%n_occ_cls
      end if
      connect2(1:2,1:nconnect) = connect(1:2,1:nconnect) + 1
      call expand_op_product2(form_list,idx_res,
     &       fac,nops+2,nops+1,
     &       idx_op_vtx,idx_sv_vtx,
     &       iblk_min2,iblk_max2,
     &       connect2,nconnect,
     &       0,0,
     &       0,0,
     &       .false.,op_info)

******OBSOLETE*********************************************************
*     obsolete part of subroutine
*     now deleted
***********************************************************************

      return
      end

