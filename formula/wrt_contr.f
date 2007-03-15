
*----------------------------------------------------------------------*
      subroutine wrt_contr(lu,contr)
*----------------------------------------------------------------------*
*     write formula on unit lu in condensed form (assuming that 
*     usually word numbers (<32768) are sufficient)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     lu
      type(contraction), intent(in) ::
     &     contr

      integer, parameter ::
     &     lbuf = 1024
      integer(2) ::
     &     buffer(lbuf)
      integer ::
     &     idx, ierr, ii, ica, igastp
      
      ierr = 0
      buffer(1) = contr%iblk_res
      buffer(2) = contr%nvtx
      buffer(3) = contr%narc
      buffer(4) = contr%nfac
      if (buffer(1).ne.contr%iblk_res) ierr = ierr+1
      if (buffer(2).ne.contr%nvtx) ierr = ierr+1
      if (buffer(3).ne.contr%narc) ierr = ierr+1
      if (buffer(4).ne.contr%nfac) ierr = ierr+1

      idx = 4
      if (ierr.gt.0) goto 101
      if (idx+contr%nvtx*2.gt.lbuf) goto 103

      do ii = 1, contr%nvtx
        buffer(idx+1) = contr%vertex(ii)%idx_op
        if (buffer(idx+1).ne.contr%vertex(ii)%idx_op) ierr = ierr+1
        buffer(idx+2) = contr%vertex(ii)%iblk_op
        if (buffer(idx+2).ne.contr%vertex(ii)%iblk_op) ierr = ierr+1
        idx = idx+2
      end do
      if (ierr.gt.0) goto 102

      if (idx+contr%narc*8.gt.lbuf) goto 103
      do ii = 1, contr%narc
        buffer(idx+1) = contr%arc(ii)%link(1)
        buffer(idx+2) = contr%arc(ii)%link(2)
c new
        idx = idx+2
        do ica = 1, 2
          do igastp = 1, ngastp
            idx = idx+1
            buffer(idx) = contr%arc(ii)%occ_cnt(igastp,ica)
          end do
        end do
c end new
c        buffer(idx+3) = contr%arc(ii)%occ_cnt(1,1)
c        buffer(idx+4) = contr%arc(ii)%occ_cnt(2,1)
c        buffer(idx+5) = contr%arc(ii)%occ_cnt(3,1)
c        buffer(idx+6) = contr%arc(ii)%occ_cnt(1,2)
c        buffer(idx+7) = contr%arc(ii)%occ_cnt(2,2)
c        buffer(idx+8) = contr%arc(ii)%occ_cnt(3,2)
c        idx = idx+8
      end do

      if (idx+contr%nfac*3.gt.lbuf) goto 103
      do ii = 1, contr%nfac
        buffer(idx+1) = contr%inffac(1,ii)
        buffer(idx+2) = contr%inffac(2,ii)
        buffer(idx+3) = contr%inffac(3,ii)
        buffer(idx+4) = contr%inffac(4,ii)
        idx = idx+4
      end do

      write(lu) contr%fac,idx,buffer(1:idx)

      if (ntest.ge.100) then
        write(luout,*) 'wrote ',8+4+idx,' bytes'
      end if

      return

      ! error handling
 101  write(luout,*) 'wrt_contr: too large nvtx, narc'
      goto 200
 102  write(luout,*)
     &     'wrt_contr: too large numbers in vertex description'
      goto 200
 103  write(luout,*) 'wrt_contr: insufficient buffer size'

 200  call quit(1,'wrt_contr','error writing contraction')

      end
