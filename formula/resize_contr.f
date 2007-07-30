*----------------------------------------------------------------------*
      subroutine resize_contr(contr,nvtx,narc,nfac)
*----------------------------------------------------------------------*
*     check max dimensions of contr (mxvtx,mxarc,mxfac) and resize
*     if necessary
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
     &     nvtx, narc, nfac

      type(cntr_vtx), pointer ::
     &     vtx_new(:)
      type(cntr_arc), pointer ::
     &     arc_new(:)
      integer, pointer ::
     &     inf_new(:,:)

      integer ::
     &     nsave

      if (contr%mxvtx.gt.0.and..not.associated(contr%vertex))
     &     call quit(1,'resize_contr','vertex pointer inconsistent')
      if (contr%mxarc.gt.0.and..not.associated(contr%arc))
     &     call quit(1,'resize_contr','arc pointer inconsistent')
      if (contr%mxfac.gt.0.and..not.associated(contr%inffac))
     &     call quit(1,'resize_contr','inffac pointer inconsistent')
      if (contr%mxvtx.lt.nvtx) then
        allocate(vtx_new(nvtx))
        nsave = min(contr%mxvtx,contr%nvtx)
        if (nsave.gt.0)
     &       vtx_new(1:nsave) = contr%vertex(1:nsave)
        if (contr%mxvtx.gt.0) deallocate(contr%vertex)
        contr%vertex => vtx_new
        contr%mxvtx = nvtx
      end if
      
      if (contr%mxarc.lt.narc) then
        allocate(arc_new(narc))
        nsave = min(contr%mxarc,contr%narc)
        if (nsave.gt.0)
     &       arc_new(1:nsave) = contr%arc(1:nsave)
        if (contr%mxarc.gt.0) deallocate(contr%arc)
        contr%arc => arc_new
        contr%mxarc = narc
      end if

      if (contr%mxfac.lt.nfac) then
        allocate(inf_new(ld_inffac,nfac))
        nsave = min(contr%mxfac,contr%nfac)
        if (nsave.gt.0)
     &       inf_new(1:ld_inffac,1:nsave) =
     &       contr%inffac(1:ld_inffac,1:nsave)
        if (contr%mxfac.gt.0) deallocate(contr%inffac)
        contr%inffac => inf_new
        contr%mxarc = nfac
      end if

      return
      end
