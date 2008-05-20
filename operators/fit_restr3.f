*----------------------------------------------------------------------*
      subroutine fit_restr3(irst_ex,irst_cnt,
     &     iocc_op,iocc_ex,iocc_cnt,
     &     irst_op,merge_map,
     &     njoined_op,njoined_cnt,
     &     ihpvgas,ngas,nspin)
*----------------------------------------------------------------------*
*     one more dummy routine to deal with restrictions
*     given the occupations and restrictions of the operator 
*     and the occupations of external and the contraction part
*     decide on the contractions applied to external and contraction
*     part
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'stdunit.h'

      integer, intent(in) ::
     &     njoined_op, njoined_cnt, ngas, nspin,
     &     ihpvgas(ngas,nspin),
     &     merge_map(*),
     &     iocc_op(ngastp,2,njoined_op),
     &     iocc_ex(ngastp,2,njoined_op),     
     &     iocc_cnt(ngastp,2,njoined_cnt),
     &     irst_op(2,ngas,2,2,nspin,njoined_op)
      integer, intent(out) ::
     &     irst_cnt(2,ngas,2,2,nspin,njoined_cnt),
     &     irst_ex (2,ngas,2,2,nspin,njoined_op )

      integer ::
     &     idx_base, idx, iop, iex, icnt, nex, ncnt

c dbg
c      print *,'OP'
c      call wrt_occ_n(luout,iocc_op,njoined_op)
c      print *,'EX'
c      call wrt_occ_n(luout,iocc_ex,njoined_op)
c      print *,'CNT'
c      call wrt_occ_n(luout,iocc_cnt,njoined_cnt)
c dbg
      ! loop over merge map info
      idx_base = 1
      ! for each vertex of the operator ...
      do iop = 1, njoined_op
        ! ... we have the number ...
        ncnt = merge_map(idx_base)
        if (ncnt.gt.0) then
          do idx = 1, ncnt
            ! ... and the indices of the correspondig cnt vertices
            icnt = merge_map(idx_base+idx)
            ! still call for adaptation to nspin>1:
c dbg
c            print *,'OP, CNT: ',iop,icnt
c dbg
            call fit_restr(irst_cnt(1,1,1,1,1,icnt),
     &                     iocc_cnt(1,1,icnt),
     &                     irst_op (1,1,1,1,1,iop),
     &                     ihpvgas,ngas)
          end do
        end if
        idx_base = idx_base+ncnt+1
        ! ... and likewise for ex ...
        nex = merge_map(idx_base)
        if (nex.gt.0) then
          do idx = 1, nex
            iex = merge_map(idx_base+idx)
c dbg
c            print *,'OP, EX: ',iop,iex
c dbg
            ! still call for adaptation to nspin>1:
            call fit_restr(irst_ex(1,1,1,1,1,iex),
     &                     iocc_ex(1,1,iex),
     &                     irst_op(1,1,1,1,1,iop),
     &                     ihpvgas,ngas)
          end do
        end if
        idx_base = idx_base+nex+1
      end do
      
      return
      end
