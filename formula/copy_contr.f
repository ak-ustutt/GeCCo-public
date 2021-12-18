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
     &     nvtx, narc, nxarc, nfac

      nvtx = contr_src%nvtx
      narc = contr_src%narc
      nxarc = contr_src%nxarc
      nfac = contr_src%nfac
      call resize_contr(contr_tgt,nvtx,narc,nxarc,nfac)

      ! copy all elements
      ! contr_tgt = contr_src does not do the proper work:
      ! it would copy the pointer adresses for vertex etc.
      ! rather than copying the entries
      contr_tgt%idx_res = contr_src%idx_res
      contr_tgt%iblk_res = contr_src%iblk_res
      contr_tgt%dagger = contr_src%dagger
      contr_tgt%fac = contr_src%fac
      contr_tgt%nvtx = contr_src%nvtx
      contr_tgt%nsupvtx = contr_src%nsupvtx
      contr_tgt%narc = contr_src%narc
      contr_tgt%nxarc = contr_src%nxarc
      contr_tgt%nfac = contr_src%nfac

      if (nvtx.gt.0) then
        contr_tgt%vertex(1:nvtx) = contr_src%vertex(1:nvtx)
        contr_tgt%svertex(1:nvtx) = contr_src%svertex(1:nvtx)
        contr_tgt%joined(0:nvtx,1:nvtx) =contr_src%joined(0:nvtx,1:nvtx)
      end if

      if (narc.gt.0) then
        contr_tgt%arc(1:narc) = contr_src%arc(1:narc)
      end if

      if (nxarc.gt.0) then
        contr_tgt%xarc(1:nxarc) = contr_src%xarc(1:nxarc)
      end if

      if (nfac.gt.0) then
        contr_tgt%inffac(1:ld_inffac,1:nfac) =
     &       contr_src%inffac(1:ld_inffac,1:nfac)
      end if

      if (contr_src%index_info) then
        contr_tgt%index_info = .true.
        if (associated(contr_tgt%contr_string))
     &       deallocate(contr_tgt%contr_string)
        if (associated(contr_tgt%result_string))
     &       deallocate(contr_tgt%result_string)
        allocate(contr_tgt%contr_string(contr_src%nidx))
        if (contr_src%nxidx.gt.0)
     &       allocate(contr_tgt%result_string(contr_src%nxidx))
        contr_tgt%nidx = contr_src%nidx
        contr_tgt%contr_string = contr_src%contr_string
        contr_tgt%nxidx = contr_src%nxidx
        if (contr_src%nxidx.gt.0)
     &       contr_tgt%result_string = contr_src%result_string
      end if
      contr_tgt%total_sign = contr_src%total_sign
      contr_tgt%eqvl_fact = contr_src%eqvl_fact
      
      return
      end
