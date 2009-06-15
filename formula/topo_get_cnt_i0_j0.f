*----------------------------------------------------------------------*
      subroutine topo_get_cnt_i0_j0(
     &     occ_cnt,occ_i0,occ_j0,n_enclosed,
     &     ivtx,jvtx,topo,xlines,nvtx,nj_res)
*----------------------------------------------------------------------*
*     extract occupations of CNT, I0 and J0 from topo and xlines
*     also get length of enclosed string (passive vertices)
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'

      integer, intent(out) ::
     &     occ_cnt(ngastp,2), occ_i0(ngastp,2), occ_j0(ngastp,2),
     &     n_enclosed
      integer, intent(in) ::
     &     ivtx, jvtx, nvtx, nj_res
      integer(8), intent(in) ::
     &     topo(nvtx,nvtx), xlines(nvtx,nj_res)
        
      integer(8) ::
     &     base
      integer ::
     &     idx, nidx
      integer ::
     &     occ_scr(ngastp,2)

      integer, external ::
     &     int8_expand
  
      base = pack_base
      if (ivtx.ge.jvtx)
     &     call quit(1,'topo_get_cnt_i0_j0','ivtx.ge.jvtx??')

      ! get CNT
      occ_cnt = 0
      nidx = int8_expand(topo(ivtx,jvtx),base,occ_cnt)
      if (nidx.gt.ngastp*2)
     &     call quit(1,'topo_get_cnt_i0_j0','range 1')

      ! assemble I0
      call topo_occ_from_row(occ_i0,ivtx,jvtx,topo,xlines,nvtx,nj_res)

      ! assemble J0
      call topo_occ_from_row(occ_j0,jvtx,ivtx,topo,xlines,nvtx,nj_res)

      ! count enclosed CA
      n_enclosed = 0
      do idx = ivtx+1, jvtx-1
        call topo_occ_from_row(occ_scr,idx,0,topo,xlines,nvtx,nj_res)
        n_enclosed = n_enclosed + sum(occ_scr(1:ngastp,1:2))
      end do

      return

      end
