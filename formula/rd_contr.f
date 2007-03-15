
*----------------------------------------------------------------------*
      logical function rd_contr(lu,contr,idx_res)
*----------------------------------------------------------------------*
*     read next contraction from lu and expand to long form on contr
*     return .false. if EOF
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lu, idx_res
      type(contraction), intent(out) ::
     &     contr

      integer, parameter ::
     &     lbuf = 1024
      integer(2) ::
     &     buffer(lbuf)
      integer ::
     &     idx, ierr, ii, ica, igastp


      rd_contr = .true.
      
      read(lu,end=100) contr%fac,idx,buffer(1:idx)

      contr%idx_res = idx_res
      contr%iblk_res = buffer(1)
      contr%nvtx = buffer(2)
      contr%narc = buffer(3) 
      contr%nfac = buffer(4) 

      idx = 4

      do ii = 1, contr%nvtx
        contr%vertex(ii)%idx_op  = buffer(idx+1)
        contr%vertex(ii)%iblk_op = buffer(idx+2)
        idx = idx+2
      end do

      do ii = 1, contr%narc
        contr%arc(ii)%link(1)      = buffer(idx+1)
        contr%arc(ii)%link(2)      = buffer(idx+2)
c new
        idx = idx+2
        do ica = 1,2
          do igastp = 1, ngastp
            idx = idx+1
            contr%arc(ii)%occ_cnt(igastp,ica) = buffer(idx)
          end do
        end do
c end new
c        contr%arc(ii)%occ_cnt(1,1) = buffer(idx+3)
c        contr%arc(ii)%occ_cnt(2,1) = buffer(idx+4)
c        contr%arc(ii)%occ_cnt(3,1) = buffer(idx+5)
c        contr%arc(ii)%occ_cnt(1,2) = buffer(idx+6)
c        contr%arc(ii)%occ_cnt(2,2) = buffer(idx+7)
c        contr%arc(ii)%occ_cnt(3,2) = buffer(idx+8)
c        idx = idx+8
      end do

      do ii = 1, contr%nfac
        contr%inffac(1,ii) = buffer(idx+1)
        contr%inffac(2,ii) = buffer(idx+2)
        contr%inffac(3,ii) = buffer(idx+3)
        contr%inffac(4,ii) = buffer(idx+4)
        idx = idx+4
      end do

      return

      ! EOF encountered:
 100  rd_contr = .false.
      return

      end
