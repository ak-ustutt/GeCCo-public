*----------------------------------------------------------------------*
      subroutine expand_op_bch(form_res,max_n,idx_res,
     &     fac,idx_proj,idx_a,fac_b,idx_b,iblk_b_min,iblk_b_max,op_info)
*----------------------------------------------------------------------*
*     generate
*       f <res_C| Proj e^{-bB} A e^{bB}   |res_A> =
*      f <res_C| Proj (A + b[A,B] + 1/2b^2[[A,B],B] + ...) |res_A>
*
*     up to max_n-fold commutator
*
*     usually, either idxop_res will be a scalar and Proj generates
*     a projection manifold (e.g. yield directly the CC-Lagrangian)
*     or idxop_res is chosen to be the projection (yield e.g. directly
*     the CC-residual)
*
*     andreas, June 2007
*
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'ifc_operators.h'

      type(formula_item), target ::
     &     form_res
      real(8), intent(in) ::
     &     fac, fac_b
      integer, intent(in) ::
     &     idx_a, idx_b, iblk_b_min, iblk_b_max,
     &     idx_proj, idx_res, max_n
      type(operator_info), intent(in) ::
     &     op_info
      
      logical ::
     &     pure_ex, pure_dx
      real(8) ::
     &     fac_t
      integer ::
     &     iop_typ_b, n_proj, n
      type(contraction) ::
     &     proto
      type(formula_item), pointer ::
     &     fres_pnt
      integer, pointer ::
     &     connect(:,:), idx_op(:), iblk_min(:), iblk_max(:)
      type(formula_item), pointer ::
     &     form_pnt

      integer, external ::
     &     vtx_type

      form_pnt => form_res

      ! check type of operator b
      iop_typ_b = vtx_type(op_info%op_arr(idx_b)%op)

      if (iop_typ_b.ne.vtxtyp_ph_ex.and.iop_typ_b.ne.vtxtyp_ph_dx)
     &     call quit(1,'expand_op_bch',
     &     'currently, only pure EX/DX operators allowed')

      pure_ex = iop_typ_b.eq.vtxtyp_ph_ex
      pure_dx = iop_typ_b.eq.vtxtyp_ph_dx

      allocate(connect(2,max_n),idx_op(max_n+2),
     &     iblk_min(max_n+2),iblk_max(max_n+2))
      if (idx_proj.gt.0) then
        idx_op(1) = idx_proj
        iblk_min(1) = 1
        iblk_max(1) = 0 ! => all blocks
        idx_op(2) = idx_a
        iblk_min(2) = 1
        iblk_max(2) = 0 ! => all blocks
        n_proj = 1
      else
        idx_op(1) = idx_a
        iblk_min(1) = 1
        iblk_max(1) = 0 ! => all blocks
        n_proj = 0
      end if
      fac_t = fac
      do n = 0, max_n
        if (n.gt.0) then
          fac_t = fac_t*fac
          if (pure_ex) then
            connect(1,n) = n_proj+1
            connect(2,n) = n_proj+n+1
            idx_op(n_proj+n+1) = idx_b
            iblk_min(n_proj+n+1) = iblk_b_min
            iblk_max(n_proj+n+1) = iblk_b_max
          else if (pure_dx) then
            connect(1,n) = n_proj+1
            connect(2,n) = n_proj+n+1
            idx_op(n_proj+n) = idx_b
            iblk_min(n_proj+n) = iblk_b_min
            iblk_max(n_proj+n) = iblk_b_max
            idx_op(n_proj+n+1) = idx_a
            iblk_min(n_proj+n+1) = 1
            iblk_max(n_proj+n+1) = 0
          end if
        end if
c dbg
c        print *,'n = ',n
c dbg

        call expand_op_product(form_pnt,idx_res,
     &       fac_t,n_proj+n+1,idx_op,
     &       iblk_min,iblk_max,
     &       connect,n,
     &       op_info)

        ! check whether any term was created
        if (form_pnt%command.eq.command_end_of_formula) exit
        do while(form_pnt%command.ne.command_end_of_formula)
          form_pnt => form_pnt%next
        end do

      end do

      return
      end
