*----------------------------------------------------------------------*
      subroutine topo_set_unique(contr)
*----------------------------------------------------------------------*
*     set unique vtx, topo, xlines information
*     matthias, july 2010
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      type(contraction), intent(inout) ::
     &     contr

      integer ::
     &     nvtx, nj, idx, jdx
      integer(8) ::
     &     hash, ifac

      integer, pointer ::
     &     scr(:), svtx(:)

      integer, external ::
     &     njres_contr
      integer(8), external ::
     &     topo_hash

      nj = njres_contr(contr)
      nvtx = contr%nvtx
      contr%unique_set = .true.

      allocate(contr%vtx(nvtx),contr%topo(nvtx,nvtx),
     &         contr%xlines(nvtx,nj),scr(nvtx),svtx(nvtx))
      call pack_contr(scr,contr%vtx,contr%topo,contr%xlines,contr,nj)
      svtx = contr%svertex
      call topo_make_unique2(scr,contr%vtx,svtx,contr%topo,contr%xlines,
     &                       nvtx,nj)

      deallocate(scr,svtx)

      ! set hash value
      contr%hash = topo_hash(contr%vtx,contr%topo,contr%xlines,nvtx,nj)

      return
      end
