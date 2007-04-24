*----------------------------------------------------------------------*
      integer function idx_msgmdst(iocc_cls,mst,igamt,
     &     msd,gmd,dagger,op,nsym)
*----------------------------------------------------------------------*
*     note that mst, igamt are the ms_a and igam_a of the current
*     operator block
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      integer, intent(in) ::
     &     iocc_cls, mst, igamt, nsym, msd(ngastp,2), gmd(ngastp,2)
      
      type(operator), intent(in) ::
     &     op

      logical, intent(in) ::
     &     dagger

      integer ::
     &     mgdid, idx, idxms
      integer ::
     &     iocc(ngastp,2), msdd(ngastp,2), gmdd(ngastp,2)

      integer, pointer ::
     &     didarr(:,:,:)      

      integer, external ::
     &     msgmdid, ielsum

      ! get occupation
      iocc(1:ngastp,1:2) = op%ihpvca_occ(1:ngastp,1:2,iocc_cls) 

      if (.not.dagger) then
        ! calculate integer-valued ID of (ms,Gamma) distribution
        mgdid = msgmdid(iocc,msd,gmd,nsym)
      else
        msdd(1:ngastp,1) = msd(1:ngastp,2)
        msdd(1:ngastp,2) = msd(1:ngastp,1)
        gmdd(1:ngastp,1) = gmd(1:ngastp,2)
        gmdd(1:ngastp,2) = gmd(1:ngastp,1)
        mgdid = msgmdid(iocc,msdd,gmdd,nsym)
      end if

      didarr => op%off_op_gmox(iocc_cls)%did

      idxms = (ielsum(iocc,ngastp)-mst)/2+1

      idx_msgmdst = -1
      do idx = 1, op%off_op_gmox(iocc_cls)%ndis(igamt,idxms)
        if (didarr(idx,igamt,idxms).eq.mgdid) then
          idx_msgmdst = idx
          exit
        end if
      end do

      return
      end
