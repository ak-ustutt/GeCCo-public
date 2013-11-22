*----------------------------------------------------------------------*
      subroutine prt_contr_p(lulog,svtx,vtx,topo,xlines,nvtx,nj)
*----------------------------------------------------------------------*
*     print contraction in topo-matrix form
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     lulog, nvtx, nj

      integer ::
     &     svtx(nvtx)
      integer(8), intent(in) ::
     &     vtx(nvtx), topo(nvtx,nvtx),
     &     xlines(nvtx,nj)

      integer ::
     &     ivtx
      character(32) ::
     &     fmt

      write(fmt,'("(x,i2,i5,""|"",",i2,"i9.8,""|"",",i1,"i9.8",")")')
     &     nvtx, nj
      if (nvtx.gt.0.and.svtx(1).gt.0) then
        do ivtx = 1, nvtx
          write(lulog,fmt) svtx(ivtx), vtx(ivtx),
     &         topo(ivtx,1:nvtx), xlines(ivtx,1:nj)
        end do
      else
        do ivtx = 1, nvtx
          write(lulog,fmt) 0, vtx(ivtx),
     &         topo(ivtx,1:nvtx), xlines(ivtx,1:nj)
        end do
      end if

      return
      end
