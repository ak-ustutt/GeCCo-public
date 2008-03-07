*----------------------------------------------------------------------*
      logical function restr_cmp(irst_op,irst_str,ica,igastp,
     &     ihpvgas,ngas,nspin)
*----------------------------------------------------------------------*
*
*     compare the part of restriction on operator block given by
*     ica and igastp to restriction on string
*
*     return true if equal, i.e. if this string describes the given
*     part of the operator
*
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     ngas, nspin, igastp, ica,
     &     irst_op(2,ngas,2,2,nspin),
     &     irst_str(2,ngas,2,nspin), ihpvgas(ngas,nspin)

      logical ::
     &     same, inc_j
      integer ::
     &     igas, jgas, ispin

      jgas = 1 ! counter within irst_str
      same = .true.
      cmp_loop: do igas = 1, ngas
        do ispin = 1, nspin
          inc_j = .false.
          if (ihpvgas(igas,ispin).ne.igastp) cycle
c dbg
c          print *,'igas, jgas, ispin: ',igas, jgas, ispin
c          print '(x,a,4i4)',
c     &         'comparing: ',irst_op(1:2,igas,ica,1:2,ispin)
c          print '(x,a,4i4)',
c     &         '       to: ',irst_str(1:2,jgas,1:2,ispin)
c dbg          
          inc_j = .true.
          same = same.and.
     &       (irst_op(1,igas,ica,1,ispin).eq.
     &       irst_str(1,jgas,1,ispin)).and.
     &       (irst_op(2,igas,ica,1,ispin).eq.
     &       irst_str(2,jgas,1,ispin)) .and.
     &       (irst_op(1,igas,ica,2,ispin).eq.
     &       irst_str(1,jgas,2,ispin)).and.
     &       (irst_op(2,igas,ica,2,ispin).eq.
     &       irst_str(2,jgas,2,ispin))
          if (.not.same) exit cmp_loop
        end do
        if (inc_j) jgas = jgas+1
      end do cmp_loop

      restr_cmp = same
c dbg
c      print *,'result: ',same
c dbg

      return
      end

