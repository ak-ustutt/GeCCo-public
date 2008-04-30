      subroutine prt_contr_p(luout,svtx,vtx,topo,xlines,nvtx,nj)

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     luout, nvtx, nj

      integer ::
     &     svtx(nvtx)
      integer(8), intent(in) ::
     &     vtx(nvtx), topo(nvtx,nvtx),
     &     xlines(nvtx,nj)

      integer ::
     &     ivtx
      character(32) ::
     &     fmt

      write(fmt,'("(x,i2,i5,""|"",",i1,"i9.8,""|"",",i1,"i9.8",")")')
     &     nvtx, nj
      do ivtx = 1, nvtx
        write(luout,fmt) svtx(ivtx), vtx(ivtx),
     &       topo(ivtx,1:nvtx), xlines(ivtx,1:nj)
      end do

      return
      end
