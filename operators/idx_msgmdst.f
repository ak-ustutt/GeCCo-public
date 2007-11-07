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

      logical ::
     &     dag_temp

      integer, external ::
     &     msgmdid, ielsum

      ! get occupation
      iocc(1:ngastp,1:2) = op%ihpvca_occ(1:ngastp,1:2,iocc_cls) 

c bodge
      if(op%njoined.gt.1)then
        iocc(1:ngastp,1:2)=iocc(1:ngastp,1:2)+
     &       op%ihpvca_occ(1:ngastp,1:2,2)
      endif
      dag_temp = dagger
      if(op%dagger) dag_temp = .not.dag_temp
c bodge

c      if (.not.dagger) then
      if (.not.dag_temp) then
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
c dbg
c      print *,'-->',iocc_cls,igamt,idxms
c dbg
      do idx = 1, op%off_op_gmox(iocc_cls)%ndis(igamt,idxms)
        if (didarr(idx,igamt,idxms).eq.mgdid) then
          idx_msgmdst = idx
          exit
        end if
      end do

      return
      end
