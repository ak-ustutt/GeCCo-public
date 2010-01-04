*----------------------------------------------------------------------*
      integer function idx_msgmdst2(
     &     iblk,idxmsa_blk,gama_blk,
     &     occ_c,idxms_c,gam_c,nc,
     &     occ_a,idxms_a,gam_a,na,
     &     dagger,tra_map_c,tra_map_a,mel,nsym)
*----------------------------------------------------------------------*
*     return the distribution number for given block, MS(A), GAMMA(A)
*     version for condensed representation
*
*     FIX for transposed operators with multiple vertices
*     tra_map_c/tra_map_a contain the necessary reordering arrays
*     for the time being, we allow an entry of "-1" to indicate that
*     the maps have not been set
*     (they can be obtained from set_dis_tra_map())
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
     &     occ_a(na), idxms_a(na), gam_a(na),
     &     tra_map_c(nc), tra_map_a(na)
      
      type(me_list), intent(in) ::
     &     mel

      logical, intent(in) ::
     &     dagger

      integer ::
     &     mgdid, idx, idx_end

      integer, pointer ::
     &     didarr(:,:,:)      
      integer ::
     &     occ_c_tra(nc), gam_c_tra(nc), idxms_c_tra(nc),
     &     occ_a_tra(na), gam_a_tra(na), idxms_a_tra(na)

      logical ::
     &     dag_temp

      integer, external ::
     &     msgmdid2, ielsum

      dag_temp = dagger
      if(mel%op%dagger) !dag_temp = .not.dag_temp
     &     call quit(1,'idx_msgmdst2',
     &     'op%dagger should not be set to .true. any more....')

      if (.not.dag_temp) then
        ! calculate integer-valued ID of (ms,Gamma) distribution
        mgdid = msgmdid2(occ_c,idxms_c,gam_c,nc,
     &                   occ_a,idxms_a,gam_a,na,nsym)
      else
        if (na.le.1.and.nc.le.1) then
          mgdid = msgmdid2(occ_a,idxms_a,gam_a,na,
     &                     occ_c,idxms_c,gam_c,nc,nsym)
        else
          if (nc.gt.0.and.tra_map_c(1).lt.0
     &        .or.na.gt.0.and.tra_map_a(1).lt.0)
     &         call quit(1,'idx_msgmdst2','need to set tra_map_c/a')
          occ_c_tra = occ_c(tra_map_c(1:nc))
          occ_a_tra = occ_a(tra_map_a(1:na))
          idxms_c_tra = idxms_c(tra_map_c(1:nc))
          idxms_a_tra = idxms_a(tra_map_a(1:na))
          gam_c_tra = gam_c(tra_map_c(1:nc))
          gam_a_tra = gam_a(tra_map_a(1:na))
          mgdid = msgmdid2(occ_a_tra,idxms_a_tra,gam_a_tra,na,
     &                     occ_c_tra,idxms_c_tra,gam_c_tra,nc,nsym)
        end if
      end if

      didarr => mel%off_op_gmox(iblk)%did

      idx_msgmdst2 = -1

      idx_end = mel%off_op_gmox(iblk)%ndis(gama_blk,idxmsa_blk)
c dbg
c      print *,'-->',gama_blk,idxmsa_blk,idx_end
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
        write(luout,*) 'nc,na    ',nc,na
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

c dbg
c        write(luout,*) 'list: ',trim(mel%label)
c        write(luout,*) 'op:   ',trim(mel%op%name)
c        write(luout,*) 'gamma(A),idxms(A):',gama_blk,idxmsa_blk
c        write(luout,*) 'nc,na    ',nc,na
c        write(luout,*) 'occ_c:   ',occ_c(1:nc)
c        write(luout,*) 'idxms_c: ',idxms_c(1:nc)
c        write(luout,*) 'gam_c:   ',gam_c(1:nc)
c        write(luout,*) 'occ_a:   ',occ_a(1:na)
c        write(luout,*) 'idxms_a: ',idxms_a(1:na)
c        write(luout,*) 'gam_a:   ',gam_a(1:na)
c        write(luout,*) 'mgdid: ',mgdid
c        write(luout,*) 'didarr:',
c     &       didarr(1:idx_end,gama_blk,idxmsa_blk)
c dbg

      return
      end
