*----------------------------------------------------------------------*
      subroutine copy_contr(contr_src,contr_tgt)
*----------------------------------------------------------------------*
*     copy all information from src to tgt
*     the target is resized, if necessary     
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr_src
      type(contraction), intent(inout) ::
     &     contr_tgt

      integer ::
     &     nvtx, narc, nfac

      nvtx = contr_src%nvtx
      narc = contr_src%narc
      nfac = contr_src%nfac
      call resize_contr(contr_tgt,nvtx,narc,nfac)

      ! copy all elements
      ! contr_tgt = contr_src does not do the proper work:
      ! it would copy the pointer adresses for vertex etc.
      ! rather than copying the entries
      contr_tgt%idx_res = contr_src%idx_res
      contr_tgt%iblk_res = contr_src%iblk_res
      contr_tgt%fac = contr_src%fac
      contr_tgt%nvtx = contr_src%nvtx
      contr_tgt%narc = contr_src%narc
      contr_tgt%nfac = contr_src%nfac

      if (nvtx.gt.0) then
        contr_tgt%vertex(1:nvtx) = contr_src%vertex(1:nvtx)
      end if

      if (narc.gt.0) then
        contr_tgt%arc(1:narc) = contr_src%arc(1:narc)
      end if

      if (nfac.gt.0) then
        contr_tgt%inffac(1:ld_inffac,1:nfac) =
     &       contr_src%inffac(1:ld_inffac,1:nfac)
      end if

      return
      end
