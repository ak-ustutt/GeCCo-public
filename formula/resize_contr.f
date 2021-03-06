*----------------------------------------------------------------------*
      subroutine resize_contr(contr,nvtx,narc,nxarc,nfac)
*----------------------------------------------------------------------*
*     check max dimensions of contr (mxvtx,mxarc,mxxarc,mxfac) and
*     resize if necessary
*     data-save version:
*     we first allocate new memory, deallocate the old one, and
*     finally redirect the pointer to the new array
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(inout) ::
     &     contr
      integer, intent(in) ::
     &     nvtx, narc, nxarc, nfac

      type(cntr_vtx), pointer ::
     &     vtx_new(:)
      type(cntr_arc), pointer ::
     &     arc_new(:),
     &     xarc_new(:)
      integer, pointer ::
     &     inf_new(:,:), joined_new(:,:), svertex_new(:)

      integer ::
     &     nsave

      if (contr%mxvtx.gt.0.and..not.associated(contr%vertex))
     &     call quit(1,'resize_contr','vertex pointer inconsistent')
      if (contr%mxvtx.gt.0.and..not.associated(contr%joined))
     &     call quit(1,'resize_contr','joined pointer inconsistent')
      if (contr%mxvtx.gt.0.and..not.associated(contr%svertex))
     &     call quit(1,'resize_contr','svertex pointer inconsistent')
      if (contr%mxarc.gt.0.and..not.associated(contr%arc))
     &     call quit(1,'resize_contr','arc pointer inconsistent')
      if (contr%mxxarc.gt.0.and..not.associated(contr%xarc))
     &     call quit(1,'resize_contr','xarc pointer inconsistent')
      if (contr%mxfac.gt.0.and..not.associated(contr%inffac))
     &     call quit(1,'resize_contr','inffac pointer inconsistent')

      if (contr%mxvtx.lt.nvtx) then
        allocate(vtx_new(nvtx),
     &       joined_new(0:nvtx,nvtx),svertex_new(nvtx))
        nsave = min(contr%mxvtx,contr%nvtx)
        if (nsave.gt.0) then
          vtx_new(1:nsave) = contr%vertex(1:nsave)
          svertex_new(1:nsave) = contr%svertex(1:nsave)
          joined_new(0:nsave,1:nsave) = contr%joined(0:nsave,1:nsave)
          joined_new(nsave+1:nvtx,1:nsave) = 0
        end if
c fix for newly introduced dagger flag -- ensure that it is always def.'d
        if (nsave+1.le.nvtx) vtx_new(nsave+1:nvtx)%dagger = .false.
        if (contr%mxvtx.gt.0) then
          deallocate(contr%vertex)
          deallocate(contr%joined)
          deallocate(contr%svertex)
        end if
        contr%vertex => vtx_new
        contr%joined => joined_new
        contr%svertex => svertex_new
        contr%mxvtx = nvtx
        nullify(vtx_new,joined_new,svertex_new)
      end if

      if (contr%mxarc.lt.narc) then
        allocate(arc_new(narc))
        nsave = min(contr%mxarc,contr%narc)
        if (nsave.gt.0)
     &       arc_new(1:nsave) = contr%arc(1:nsave)
        if (contr%mxarc.gt.0) deallocate(contr%arc)
        contr%arc => arc_new
        contr%mxarc = narc
        nullify(arc_new)
      end if

      if (contr%mxxarc.lt.nxarc) then
        allocate(xarc_new(nxarc))
        nsave = min(contr%mxxarc,contr%nxarc)
        if (nsave.gt.0)
     &       xarc_new(1:nsave) = contr%xarc(1:nsave)
        if (contr%mxxarc.gt.0) deallocate(contr%xarc)
        contr%xarc => xarc_new
        contr%mxxarc = nxarc
        nullify(xarc_new)
      end if

      if (contr%mxfac.lt.nfac) then
        allocate(inf_new(ld_inffac,nfac))
        nsave = min(contr%mxfac,contr%nfac)
        if (nsave.gt.0)
     &       inf_new(1:ld_inffac,1:nsave) =
     &       contr%inffac(1:ld_inffac,1:nsave)
        if (contr%mxfac.gt.0) deallocate(contr%inffac)
        contr%inffac => inf_new
        contr%mxfac = nfac
        nullify(inf_new)
      end if

      ! unique representation will look different
      if (contr%unique_set) then
        deallocate(contr%vtx,contr%topo,contr%xlines)
        contr%unique_set = .false.
      end if
      contr%vtx => null()
      contr%topo => null()
      contr%xlines => null()

      return
      end
