*----------------------------------------------------------------------*
      integer function idx_msgmdst2(err,
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
     &     dagger, err

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
c dbg
c          write(lulog,*) 'idxms_c_tra: ',idxms_c_tra(1:nc)
c          write(lulog,*) 'idxms_a_tra: ',idxms_a_tra(1:na)
c dbg
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
c          write(lulog,*) 'mgdid: ', mgdid, ' for ',gama_blk,idxmsa_blk 
c          write(lulog,*) 'didarr:',
c     &       didarr(1:idx_end,gama_blk,idxmsa_blk)
c
c dbg
c dbg
c      print *,'-->',gama_blk,idxmsa_blk,idx_end
c dbg
      do idx = 1, idx_end
        if (didarr(idx,gama_blk,idxmsa_blk).eq.mgdid) then
          idx_msgmdst2 = idx
          exit
        end if
      end do

      if (idx_msgmdst2.eq.-1.and.err) then
c dbg
c        print *,'nsym = ',nsym
c dbg
        write(lulog,*) 'list: ',trim(mel%label)
        write(lulog,*) 'op:   ',trim(mel%op%name)
        write(lulog,*) 'gamma(A),idxms(A):',gama_blk,idxmsa_blk
        write(lulog,*) 'nc,na    ',nc,na
        write(lulog,*) 'occ_c:   ',occ_c(1:nc)
        write(lulog,*) 'idxms_c: ',idxms_c(1:nc)
        write(lulog,*) 'gam_c:   ',gam_c(1:nc)
        write(lulog,*) 'occ_a:   ',occ_a(1:na)
        write(lulog,*) 'idxms_a: ',idxms_a(1:na)
        write(lulog,*) 'gam_a:   ',gam_a(1:na)
        write(lulog,*) 'mgdid: ',mgdid
        write(lulog,*) 'didarr:',
     &       didarr(1:idx_end,gama_blk,idxmsa_blk)
        call quit(1,'idx_msgmdst2','distribution not found')
      end if

c dbg
c        write(lulog,*) 'list: ',trim(mel%label)
c        write(lulog,*) 'op:   ',trim(mel%op%name)
c        write(lulog,*) 'gamma(A),idxms(A):',gama_blk,idxmsa_blk
c        write(lulog,*) 'nc,na    ',nc,na
c        write(lulog,*) 'occ_c:   ',occ_c(1:nc)
c        write(lulog,*) 'idxms_c: ',idxms_c(1:nc)
c        write(lulog,*) 'gam_c:   ',gam_c(1:nc)
c        write(lulog,*) 'occ_a:   ',occ_a(1:na)
c        write(lulog,*) 'idxms_a: ',idxms_a(1:na)
c        write(lulog,*) 'gam_a:   ',gam_a(1:na)
c        write(lulog,*) 'mgdid: ',mgdid
c        write(lulog,*) 'didarr:',
c     &       didarr(1:idx_end,gama_blk,idxmsa_blk)
c dbg

      return
      end
