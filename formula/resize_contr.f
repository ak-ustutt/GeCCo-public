*----------------------------------------------------------------------*
      subroutine resize_contr(contr,nvtx,narc,nfac)
*----------------------------------------------------------------------*
*     check max dimensions of contr (mxvtx,mxarc,mxfac) and resize
*     if necessary
*     currently, all information previously contained on vertex, arc
*     and inffac will be lost !
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(inout) ::
     &     contr
      integer, intent(in) ::
     &     nvtx, narc, nfac

      if (contr%mxvtx.lt.nvtx) then
        if (contr%mxvtx.gt.0) deallocate(contr%vertex)
        allocate(contr%vertex(nvtx))
        contr%mxvtx = nvtx
      end if
      
      if (contr%mxarc.lt.narc) then
        if (contr%mxarc.gt.0) deallocate(contr%arc)
        allocate(contr%arc(narc))
        contr%mxarc = narc
      end if

      if (contr%mxfac.lt.nfac) then
        if (contr%mxfac.gt.0) deallocate(contr%inffac)
        allocate(contr%inffac(4,nfac))
        contr%mxfac = nfac
      end if

      return
      end
