*----------------------------------------------------------------------*
      integer function idx_msgmdst(iocc_cls,mst,igamt,
     &     msd,gmd,dagger,mel,nsym)
*----------------------------------------------------------------------*
*     note that mst, igamt are the ms_a and igam_a of the current
*     ME-list block
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_me_list.h'

      integer, intent(in) ::
     &     iocc_cls, mst, igamt, nsym, msd(ngastp,2), gmd(ngastp,2)
      
      type(me_list), intent(in) ::
     &     mel

      logical, intent(in) ::
     &     dagger

      integer ::
     &     mgdid, idx, idxms
      integer ::
     &     iocc(ngastp,2), msdd(ngastp,2), gmdd(ngastp,2)

      type(operator), pointer ::
     &     op

      integer, pointer ::
     &     didarr(:,:,:) 

      logical ::
     &     dag_temp

      integer, external ::
     &     msgmdid, ielsum

      op => mel%op

      ! get occupation
      iocc(1:ngastp,1:2) = op%ihpvca_occ(1:ngastp,1:2,iocc_cls) 

c bodge
      if(op%njoined.gt.1)then
        iocc(1:ngastp,1:2)=iocc(1:ngastp,1:2)+
     &       op%ihpvca_occ(1:ngastp,1:2,2)
        ! allow this only for call from idx42str2 (else we have idx_msgmdst2())
        ! check this way:
        if (op%njoined.gt.2.or.ielsum(iocc,ngastp*2).ne.4)
     &       call quit(1,'idx_msgmdst','you are definitely abusing me!')

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

      didarr => mel%off_op_gmox(iocc_cls)%did

      idxms = (ielsum(iocc,ngastp)-mst)/2+1

      idx_msgmdst = -1
c dbg
c      print *,'dagger = ',dag_temp
c      print *,'-->',iocc_cls,igamt,idxms
c      print *,'   ',mgdid,mgdid2
c      print *,'   ',
c     &     didarr(1:mel%off_op_gmox(iocc_cls)%ndis(igamt,idxms),
c     &     igamt,idxms)
c dbg
      do idx = 1, mel%off_op_gmox(iocc_cls)%ndis(igamt,idxms)
        if (didarr(idx,igamt,idxms).eq.mgdid) then
          idx_msgmdst = idx
          exit
        end if
      end do

      return
      end
