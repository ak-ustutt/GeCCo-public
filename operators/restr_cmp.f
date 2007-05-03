*----------------------------------------------------------------------*
      logical function restr_cmp(irst_op,irst_str,ica,igastp,
     &     ihpvgas,ngas)
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
     &     ngas, igastp, ica,
     &     irst_op(2,ngas,2,2), irst_str(2,ngas,2), ihpvgas(ngas)

      logical ::
     &     same
      integer ::
     &     igas, jgas

      jgas = 1 ! counter within irst_str
      same = .true.
      cmp_loop: do igas = 1, ngas
        if (ihpvgas(igas).ne.igastp) cycle cmp_loop
        same = same.and.
     &       (irst_op(1,igas,ica,1).eq.
     &       irst_str(1,jgas,1)).and.
     &       (irst_op(2,igas,ica,1).eq.
     &       irst_str(2,jgas,1)) .and.
     &       (irst_op(1,igas,ica,2).eq.
     &       irst_str(1,jgas,2)).and.
     &       (irst_op(2,igas,ica,2).eq.
     &       irst_str(2,jgas,2))
        if (.not.same) exit cmp_loop
        jgas = jgas+1
      end do cmp_loop

      restr_cmp = same

      return
      end

