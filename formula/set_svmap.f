      subroutine set_svmap(svmap,occ_ol,occ_res,nvtx,nj_res,ierr)
      ! try to find a unique mapping of open lines per vtx to 
      ! result occupation
      ! primitive version 

      implicit none
      include 'opdim.h'
      include 'ifc_operators.h'

      integer, intent(in) ::
     &     nvtx, nj_res,
     &     occ_ol(ngastp,2,nvtx),
     &     occ_res(ngastp,2,nj_res)
      integer, intent(out) ::
     &     svmap(nvtx), ierr

      integer ::
     &     isvtx, ivtx
      integer ::
     &     occ_scr(ngastp,2,nj_res)

c dbg
c      call wrt_occ_n(6,occ_ol,nvtx)
c      call wrt_occ_n(6,occ_res,nj_res)
c dbg
      ierr = 0
      svmap(1:nvtx) = 0
      occ_scr = 0
      isvtx = 1
      do ivtx = 1, nvtx
c dbg
c        print *,'ivtx, isvtx: ',ivtx,isvtx
c        call wrt_occ_n(6,occ_scr,nj_res)
c        call wrt_occ_n(6,occ_res,nj_res)
c dbg
        if (iocc_equal(occ_scr(1:ngastp,1:2,isvtx),.false.,
     &                 occ_res(1:ngastp,1:2,isvtx),.false.))
     &       isvtx = isvtx+1
        if (iocc_zero(occ_ol(1:ngastp,1:2,ivtx))) cycle
        if (isvtx.gt.nj_res) call quit(1,'set_svmap','error 1')
        if (isvtx.gt.nj_res) ierr = 1!call quit(1,'set_svmap','error 1')
        if (isvtx.gt.nj_res) return

        svmap(ivtx) = isvtx
        occ_scr(1:ngastp,1:2,isvtx) =
     &       occ_scr(1:ngastp,1:2,isvtx) +
     &       occ_ol (1:ngastp,1:2,ivtx)

        if (.not.iocc_bound('<=',occ_scr(1:ngastp,1:2,isvtx),.false.,
     &                           occ_res(1:ngastp,1:2,isvtx),.false.))
     &       call quit(1,'set_svmap','error 2')
      end do

c dbg
c      print *,'svmap set to: ',svmap
c dbg
      return
      end
