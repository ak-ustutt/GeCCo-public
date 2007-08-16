*----------------------------------------------------------------------*
      integer function idx_msgmdst2(
     &     iblk,idxmsa_blk,gama_blk,
     &     occ_c,idxms_c,gam_c,nc,
     &     occ_a,idxms_a,gam_a,na,
     &     dagger,op,nsym)
*----------------------------------------------------------------------*
*     version for condensed representation
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      integer, intent(in) ::
     &     iblk, idxmsa_blk, gama_blk,
     &     nc, na, nsym,
     &     occ_c(nc), idxms_c(nc), gam_c(nc),
     &     occ_a(na), idxms_a(na), gam_a(na)
      
      type(operator), intent(in) ::
     &     op

      logical, intent(in) ::
     &     dagger

      integer ::
     &     mgdid, idx, idx_end

      integer, pointer ::
     &     didarr(:,:,:)      

      integer, external ::
     &     msgmdid2, ielsum

      if (.not.dagger) then
        ! calculate integer-valued ID of (ms,Gamma) distribution
        mgdid = msgmdid2(occ_c,idxms_c,gam_c,nc,
     &                   occ_a,idxms_a,gam_a,na,nsym)
      else
        mgdid = msgmdid2(occ_a,idxms_a,gam_a,na,
     &                   occ_c,idxms_c,gam_c,nc,nsym)
      end if

      didarr => op%off_op_gmox(iblk)%did

      idx_msgmdst2 = -1

      idx_end = op%off_op_gmox(iblk)%ndis(gama_blk,idxmsa_blk)
c dbg
c      print *,'-->',iocc_cls,igamt,idxms
c dbg
      do idx = 1, idx_end
        if (didarr(idx,gama_blk,idxmsa_blk).eq.mgdid) then
          idx_msgmdst2 = idx
          exit
        end if
      end do

      return
      end
