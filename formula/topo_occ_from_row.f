*----------------------------------------------------------------------*
      subroutine topo_occ_from_row(occ,irow,iavoid,topo,xlines,
     &                             nvtx,nj_res)
*----------------------------------------------------------------------*
*     extract the occupation from row irow in topo and xlines
*     the entry iavoid is omitted
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'

      integer, intent(out) ::
     &     occ(ngastp,2)
      integer, intent(in) ::
     &     irow, iavoid, nvtx, nj_res
      integer(8), intent(in) ::
     &     topo(nvtx,nvtx), xlines(nvtx,nj_res)

      integer(8) ::
     &     base
      integer ::
     &     kvtx, nidx
      integer ::
     &     occ_scr(ngastp,2)

      integer, external ::
     &     int8_expand


      base = pack_base
      occ = 0
      do kvtx = 1, nvtx
        if (kvtx.eq.iavoid) cycle
        if (topo(irow,kvtx).eq.0) cycle
        occ_scr = 0
        nidx = int8_expand(topo(irow,kvtx),base,occ_scr)
        if (nidx.gt.ngastp*2)
     &       call quit(1,'topo_occ_from_row','range 2')
        occ = occ + occ_scr
      end do
      do kvtx = 1, nj_res
        if (xlines(irow,kvtx).eq.0) cycle
        occ_scr = 0
        nidx = int8_expand(xlines(irow,kvtx),base,occ_scr)
        if (nidx.gt.ngastp*2)
     &       call quit(1,'topo_occ_from_row','range 3')
        occ = occ + occ_scr
      end do

      return
      end
