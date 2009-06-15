*----------------------------------------------------------------------*
      subroutine topo_contract(sign_cnt,
     &     topo,xlines,nvtx,nj_res,
     &     svertex,isvtx1,isvtx2,
     &     list,nlist)
*----------------------------------------------------------------------*
*     carry out contraction: remove all arcs from the list
*       (i.e. zero the respective entries in topo)
*       at the same time, get the sign associated with this
*       event
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     nvtx, nj_res, nlist,
     &     list(2,nlist), svertex(nvtx), isvtx1, isvtx2
      integer, intent(out) ::
     &     sign_cnt
      integer(8), intent(inout) ::
     &     topo(nvtx,nvtx), xlines(nvtx,nj_res)

      integer ::
     &     idx, ivtx, jvtx, n_enclosed

      integer ::
     &     occ_cnt(ngastp,2), occ_i0(ngastp,2), occ_j0(ngastp,2)

      integer, external ::
     &     sign_hpvx, sign_contr

      sign_cnt = 1
      do idx = 1, nlist
        ! get the two vertices to be contracted
        ivtx = list(1,idx)
        jvtx = list(2,idx)

        ! turn topo and xline information into
        !  - CNT and I0,J0 occupations at contracted vertices
        !  - length of enclosed string
        call topo_get_cnt_i0_j0(occ_cnt,occ_i0,occ_j0,n_enclosed,
     &                          ivtx,jvtx,
     &                          topo,xlines,nvtx,nj_res)

c dbg
c        print *,'idx = ',idx,ivtx,jvtx
c        print *,'CNT,I0,J0:'
c        call wrt_occ(6,occ_cnt)
c        call wrt_occ(6,occ_i0)
c        call wrt_occ(6,occ_j0)
c        print *,'nenclosed = ',n_enclosed
c dbg
        ! get the sign for
        ! {IC IA} -> {KC I0C I0A KA}
        sign_cnt = sign_cnt*sign_hpvx(1,occ_cnt,.false.,occ_i0,.false.)
c dbg
c        print *,'after hpvx I ',sign_cnt
c dbg
        ! {JC JA} -> {KA+ J0C J0A KC+}
        sign_cnt = sign_cnt*sign_hpvx(1,occ_cnt,.true. ,occ_j0,.false.)
c dbg
c        print *,'after hpvx J ',sign_cnt
c dbg
        ! approach of KC and KC+, and KA and KA+
        sign_cnt = sign_cnt*
     &       sign_contr(occ_cnt,occ_i0,occ_j0,n_enclosed,.false.)
c dbg
c        print *,'after sign_contr ',sign_cnt
c dbg

        ! remove contraction from topo
        topo(ivtx,jvtx) = 0
        topo(jvtx,ivtx) = 0
      end do

      return
      end 
