*----------------------------------------------------------------------*
      subroutine rank_vtxs(rank,vtx,topo,xlines,nvtx,nj)
*----------------------------------------------------------------------*
*     rank vertices according to
*       a) their vtx-number
*       b) their xlines (with increasing nj)
*       c) their connection to other vertices
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nvtx, nj
      integer(8), intent(in) ::
     &     vtx(nvtx), topo(nvtx,nvtx), xlines(nvtx,nj)
      integer, intent(out) ::
     &     rank(nvtx)

      integer(8) ::
     &     vtx_scr(nvtx), topo_scr(nvtx,nvtx), xlines_scr(nvtx,nj)
      integer ::
     &     ivtx,
     &     ireo(nvtx)

      vtx_scr = vtx
      topo_scr = topo
      xlines_scr = xlines
      call topo_make_unique(ireo,vtx_scr,topo_scr,xlines_scr,nvtx,nj)

      do ivtx = 1, nvtx
        rank(ireo(ivtx)) = ivtx
      end do

      if (ntest.ge.100) then
        write(luout,*) 'ranking of vertices:'
        write(luout,'(1x,10i4)') rank(1:nvtx)
      end if

      return
      end
