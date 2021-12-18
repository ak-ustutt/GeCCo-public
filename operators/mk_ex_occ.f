*----------------------------------------------------------------------*
      subroutine mk_ex_occ(occ_ex,occ_op,occ_cnt,
     &                     mmap,dagger,nj)
*----------------------------------------------------------------------*
*     generate passive ("external") index occupation for operator
*     contracted by occ_cnt, with merging information on mmap
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'

      logical, intent(in) ::
     &     dagger
      integer, intent(in) ::
     &     nj, occ_op(ngastp,2,nj), occ_cnt(ngastp,2,*),
     &     mmap(*)
      integer, intent(out) ::
     &     occ_ex(ngastp,2,nj)

      integer ::
     &     ij, ivtx, idx, jdx, nvtx1, nvtx2,
     &     ivtx_op1, ivtx_op2, ivtx_cnt,
     &     ioff_map1, ioff_map2
      logical ::
     &     error

      occ_ex = occ_op
      error = .false.

c dbg
      print *,'occ_op: '
      call wrt_occ_n(6,occ_op,nj)
c dbg

      idx = 0
      do ij = 1, nj
        nvtx1 = mmap(idx+1)
        error = error.or.(.not.dagger.and.nvtx1.ne.1)
        ivtx_op1 = mmap(idx+2)
        ioff_map1 = idx+1
c dbg
        print *,'nvtx = ',nvtx1
        print *,'   A-> ',mmap(idx+2:idx+1+nvtx1)
c dbg
        idx = idx + nvtx1 + 1
        nvtx2 = mmap(idx+1)
        error = error.or.(dagger.and.nvtx2.ne.1)
        ivtx_op2 = mmap(idx+2)
        ioff_map2 = idx+1
c dbg
        print *,'nvtx = ',nvtx2
        print *,'   B-> ',mmap(idx+2:idx+1+nvtx2)
c dbg
        if (.not.dagger) then
          do jdx = 1, nvtx2
            ivtx_cnt = mmap(ioff_map2+jdx)
            occ_ex(1:ngastp,1:2,ivtx_op1) =
     &           occ_ex(1:ngastp,1:2,ivtx_op1) -
     &           occ_cnt(1:ngastp,1:2,ivtx_cnt)
          end do
        else
          do jdx = 1, nvtx1
            ivtx_cnt = mmap(ioff_map1+jdx)
            occ_ex(1:ngastp,1,ivtx_op2) =
     &           occ_ex(1:ngastp,1,ivtx_op2) -
     &           occ_cnt(1:ngastp,2,ivtx_cnt)
            occ_ex(1:ngastp,2,ivtx_op2) =
     &           occ_ex(1:ngastp,2,ivtx_op2) -
     &           occ_cnt(1:ngastp,1,ivtx_cnt)
c dbg
      print *,'occ_ex: jdx = ',jdx
      call wrt_occ_n(6,occ_ex,nj)
c dbg
          end do
        end if
        idx = idx + nvtx2 + 1
      end do
c dbg
      print *,'occ_ex: '
      call wrt_occ_n(6,occ_ex,nj)
c dbg

      if (error) call quit(1,'mk_ex_occ','strange case')

      return
      end
