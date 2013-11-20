*----------------------------------------------------------------------*
      subroutine standard_order(vtx_reo,rank,
     &                          svertex,topo,xlines,nvtx,nj)
*----------------------------------------------------------------------*
*     given a ranking of vertices, find the allowed ordering of vertices
*     which has the lowest number = sum_i^n rank(i) * base^(n-i)
*     i.e. which comes closest to the ordering proposed by rank
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nvtx, nj, rank(nvtx), svertex(nvtx)
      integer, intent(out) ::
     &     vtx_reo(nvtx)
      integer(8), intent(in) ::
     &     topo(nvtx,nvtx), xlines(nvtx,nj)

      integer ::
     &     ivtx, jvtx, ivtx_r, jvtx_r, idx, ihlp

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'enforcing standard_order')
        write(lulog,*) 'ranking: '
        write(lulog,'(1x,1i4)') rank(1:nvtx)
      end if

      do ivtx = 1, nvtx
        vtx_reo(ivtx) = ivtx
      end do

      do jvtx = 2, nvtx
        jvtx_r = vtx_reo(jvtx)
        ivtx = jvtx-1
        ivtx_r = vtx_reo(ivtx)
c dbg
c        print *,'present: ',jvtx,' rank = ',rank(jvtx_r)
c dbg
        ! maximum "upper" position to that we may shift vertex(jvtx):
        do while(ivtx.gt.0)
          if (.not.may_commute(jvtx_r,ivtx_r)) exit
          ivtx = ivtx-1
          ivtx_r = vtx_reo(ivtx)
        end do
        ivtx = ivtx+1 ! actual "upper" position is ivtx+1
        ivtx_r = vtx_reo(ivtx)
c dbg
c        print *,'upper = ',ivtx
c dbg

        ! now find the uppermost postion such that vertex(jvtx) has lower
        ! rank than the succeding vertex
        do while(ivtx.lt.jvtx .and. rank(jvtx_r).ge.rank(ivtx_r))
          ivtx = ivtx+1
          ivtx_r = vtx_reo(ivtx)
        end do

        ! insert jvtx at position ivtx        
c dbg
c        print *,'shifting (1): ',ivtx,jvtx,
c     &       ' rank: ',rank(ivtx_r),rank(jvtx_r)
c dbg
        call shift_reo(vtx_reo,jvtx,ivtx,nvtx)

c        ihlp = vtx_reo(jvtx)
c        do idx = jvtx, ivtx+1, -1
c          vtx_reo(idx) = vtx_reo(idx-1)
c        end do
c        vtx_reo(ivtx) = ihlp
c dbg
c        print *,'new reo: ',vtx_reo(1:nvtx)
c        print *,'new rank:',rank(vtx_reo(1:nvtx))
c dbg

      end do
      ! do a second round with reversed search direction
      ! helps agains "blocking" of vertices, but not yet the
      ! perfect solution
      do jvtx = nvtx-1, 1, -1
        jvtx_r = vtx_reo(jvtx)
        ivtx = jvtx+1
        ivtx_r = vtx_reo(ivtx)
        ! if we shift ivtx_r will end up at position jvtx_r
        ! so we do that only if this lowers to total rank of
        ! the contraction
        if (rank(jvtx_r).lt.rank(ivtx_r)) cycle
c dbg
c        print *,'present: ',jvtx,' rank = ',rank(jvtx_r)
c dbg
        ! maximum "lower" position to that we may shift vertex(jvtx):
        do while(ivtx.le.nvtx)
          if (.not.may_commute(jvtx_r,ivtx_r)) exit
          ivtx = ivtx+1
          ivtx_r = vtx_reo(ivtx)
        end do
        ivtx = ivtx-1 ! actual "lower" position is ivtx-1
        ivtx_r = vtx_reo(ivtx)
c dbg
c        print *,'lower = ',ivtx
c dbg

        ! now find the lowermost postion such that vertex(jvtx) has higher
        ! rank than the succeding vertex
        do while(ivtx.gt.jvtx .and. rank(jvtx_r).lt.rank(ivtx_r))
          ivtx = ivtx-1
          ivtx_r = vtx_reo(ivtx)
        end do

        ! insert jvtx at position ivtx        
c dbg
c        print *,'shifting (2): ',ivtx,jvtx,
c     &       ' rank: ',rank(ivtx_r),rank(jvtx_r)
c dbg
        call shift_reo(vtx_reo,jvtx,ivtx,nvtx)

c dbg
c        print *,'new reo: ',vtx_reo(1:nvtx)
c        print *,'new rank:',rank(vtx_reo(1:nvtx))
c dbg
          
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'best allowed ranking: '
        write(lulog,'(1x,1i4)') rank(vtx_reo(1:nvtx))
      end if

      return

      contains

      logical function may_commute(vtx1,vtx2)
      integer, intent(in) ::
     &     vtx1,vtx2
      integer ::
     &     ij
      logical ::
     &     xl1_scr, xl2_scr, xl1_t, xl2_t, xl12_same

      may_commute = topo(vtx1,vtx2).eq.0
      may_commute = may_commute.and.svertex(vtx1).ne.svertex(vtx2)

      if (.not.may_commute) return

      xl1_t = .false.
      xl2_t = .false.
      xl12_same = .true.
      do ij = 1, nj
        xl1_scr = xlines(vtx1,ij).gt.0
        xl2_scr = xlines(vtx2,ij).gt.0
        xl1_t = xl1_t.or.xl1_scr
        xl2_t = xl2_t.or.xl2_scr
        xl12_same = xl12_same.and.(xl1_scr.eqv.xl2_scr)
      end do

      may_commute = .not.xl1_t.or..not.xl2_t.or.xl12_same

      return
      end function

      subroutine shift_reo(reo,from,to,n)

      integer, intent(in) ::
     &     from, to, n
      integer, intent(inout) ::
     &     reo(n)
      integer ::
     &     perm(n), reo_scr(n)

      integer ::
     &     i, inc

      inc = +1
      if (from.gt.to) inc = -1
      do i = 1, n
        perm(i) = i
      end do
      do i = from, to-inc, inc
        perm(i) = i+inc
      end do
      perm(to) = from

      do i = 1, n
        reo_scr(i) = reo(perm(i))
      end do

      reo = reo_scr
      return
      
      end subroutine
      
      end
