*----------------------------------------------------------------------*
      integer function idx_msgmdst2(
     &     iblk,idxmsa_blk,gama_blk,
     &     occ_c,idxms_c,gam_c,nc,
     &     occ_a,idxms_a,gam_a,na,
     &     dagger,mel,nsym)
*----------------------------------------------------------------------*
*     version for condensed representation
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'stdunit.h'

      integer, intent(in) ::
     &     iblk, idxmsa_blk, gama_blk,
     &     nc, na, nsym,
     &     occ_c(nc), idxms_c(nc), gam_c(nc),
     &     occ_a(na), idxms_a(na), gam_a(na)
      
      type(me_list), intent(in) ::
     &     mel

      logical, intent(in) ::
     &     dagger

      integer ::
     &     mgdid, idx, idx_end

      integer, pointer ::
     &     didarr(:,:,:)      

      logical ::
     &     dag_temp

      integer, external ::
     &     msgmdid2, ielsum

c GWR
c DANGER DANGER .....
      dag_temp = dagger
      if(mel%op%dagger) dag_temp = .not.dag_temp
c GWR

      if (.not.dag_temp) then
        ! calculate integer-valued ID of (ms,Gamma) distribution
        mgdid = msgmdid2(occ_c,idxms_c,gam_c,nc,
     &                   occ_a,idxms_a,gam_a,na,nsym)
      else
        mgdid = msgmdid2(occ_a,idxms_a,gam_a,na,
     &                   occ_c,idxms_c,gam_c,nc,nsym)
      end if

      didarr => mel%off_op_gmox(iblk)%did

      idx_msgmdst2 = -1

      idx_end = mel%off_op_gmox(iblk)%ndis(gama_blk,idxmsa_blk)
c dbg
c      print *,'-->',iocc_cls,igamt,idxms
c dbg
      do idx = 1, idx_end
        if (didarr(idx,gama_blk,idxmsa_blk).eq.mgdid) then
          idx_msgmdst2 = idx
          exit
        end if
      end do

      if (idx_msgmdst2.eq.-1) then
c dbg
c        print *,'nsym = ',nsym
c dbg
        write(luout,*) 'list: ',trim(mel%label)
        write(luout,*) 'op:   ',trim(mel%op%name)
        write(luout,*) 'gamma(A),idxms(A):',gama_blk,idxmsa_blk
        write(luout,*) 'occ_c:   ',occ_c(1:nc)
        write(luout,*) 'idxms_c: ',idxms_c(1:nc)
        write(luout,*) 'gam_c:   ',gam_c(1:nc)
        write(luout,*) 'occ_a:   ',occ_a(1:na)
        write(luout,*) 'idxms_a: ',idxms_a(1:na)
        write(luout,*) 'gam_a:   ',gam_a(1:na)
        write(luout,*) 'mgdid: ',mgdid
        write(luout,*) 'didarr:',
     &       didarr(1:idx_end,gama_blk,idxmsa_blk)
        call quit(1,'idx_msgmdst2','distribution not found')
      end if

      return
      end
