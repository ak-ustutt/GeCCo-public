*----------------------------------------------------------------------*
      subroutine topo_contr2(ieqvfac,resort,vtx_reo,
     &     contr,fix_vtx,op_info)
*----------------------------------------------------------------------*
*     analyze the topology of the contraction contr
*     ieqvfac is the number of identical permutations of the vertex
*     sequence
*     fix_vtx means: do not consider permutations of this node for
*       evaluation of ieqvfac; note that this does not affect the
*       resort array, i.e. fixed nodes will be sorted to their canonical
*       position as well
*     resort is set true if the sequence is not in canonical order
*     vtx_reo gives the reodering information:
*       vtx_reo(idx_reo) = idx_ori
*
*     version using more obvious topo-matrix approach
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     ieqvfac, vtx_reo(*)
      logical, intent(out) ::
     &     resort
      type(contraction), intent(in) ::
     &     contr
      logical, intent(in) ::
     &     fix_vtx(*)
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     resort_non_fixed
      integer ::
     &     nvtx, nj, ir, ivtx, jvtx, idx
      integer, pointer ::
     &     svertex(:), neqv(:), idx_eqv(:,:), rank(:)
      integer(8), pointer ::
     &     vtx(:), topo(:,:), xlines(:,:)

      integer, external ::
     &     get_eqvfac

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'topo_contr2')
      end if

      nj = op_info%op_arr(contr%idx_res)%op%njoined
      nvtx = contr%nvtx

      allocate(vtx(nvtx),topo(nvtx,nvtx),xlines(nvtx,nj),svertex(nvtx),
     &     neqv(nvtx),idx_eqv(nvtx,nvtx),rank(nvtx))

      call pack_contr(svertex,vtx,topo,xlines,contr,nj)

      call set_eqv_map(neqv,idx_eqv,
     &                vtx,svertex,topo,xlines,nvtx,nj)

      ieqvfac = get_eqvfac(neqv,fix_vtx,nvtx)

      call rank_vtxs(rank,vtx,topo,xlines,nvtx,nj)

      ! modify rank: same values for equivalent vertices
      do ivtx = 1, nvtx
        if (neqv(ivtx).le.1) cycle
        ir = rank(ivtx)
        do idx = 1, neqv(ivtx)
          jvtx = idx_eqv(idx,ivtx)
          rank(jvtx) = ir
        end do
      end do
      
      call standard_order(vtx_reo,rank,
     &                    svertex,topo,xlines,nvtx,nj)

      ! check whether resorting is necessary and 
      ! whether any non-fixed vertex was to change place
      resort = .false.
      resort_non_fixed = .false.
      do ivtx = 1, nvtx
        resort = resort.or.vtx_reo(ivtx).ne.ivtx
        resort_non_fixed = resort_non_fixed.or.
     &       (vtx_reo(ivtx).ne.ivtx.and..not.fix_vtx(ivtx))
      end do

      ! for backward compatibility: signal by negative eqvfac
      if (resort_non_fixed) ieqvfac = -ieqvfac

      deallocate(vtx,topo,xlines,svertex,neqv,idx_eqv,rank)

      if (ntest.ge.100) then
        write(luout,*) 'returning:'
        write(luout,*) 'ieqvfac = ',ieqvfac
        write(luout,*) 'resort  = ',resort
        write(luout,*) 'vtx_reo : ',vtx_reo(1:nvtx)
      end if

      return
      end
