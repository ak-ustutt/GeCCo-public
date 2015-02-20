*------------------------------------------------------------*
      subroutine copy_mel(label_res,label_in,trnsps,fac,
     &     op_info,orb_info,str_info)
*------------------------------------------------------------*
*     copy mel from label_in to label_out.
*     this routine allows to copy mel of one operator 
*     to that of its adjoint operator.
*------------------------------------------------------------*

      
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
       
      character(*), intent(in) ::
     &     label_in
      character(*), intent(inout) ::
     &     label_res
      logical, intent(in) ::
     &     trnsps 
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      real(8), intent(in) ::
     &     fac

      type(me_list), pointer ::
     &     me_res, me_current

      integer ::
     &     idx_res, idx_cur, nblk, iblk, iblkoff, njoined,
     &     iblk_cur

      logical ::
     &     open_close_res, copy_only, open_close_cur

      real(8) ::
     &     xnorm2

      integer, pointer  ::
     &     occ(:,:,:)

       integer, external ::
     &     idx_mel_list, iblk_occ


      copy_only = .True.


      idx_res = idx_mel_list(label_res,op_info)

      if (idx_res.lt.0) then
        write(lulog,*) '"',trim(label_res),'"'
        write(lulog,*) idx_res
        call quit(1,'copy_mel','label not on list')
      end if

      me_res => op_info%mel_arr(idx_res)%mel
      if (.not.associated(me_res%fhand))
     &     call quit(1,'copy_mel','no file handle defined for '//
     &                  trim(me_res%label))
      open_close_res = me_res%fhand%unit.le.0

      idx_cur = idx_mel_list(label_in,op_info)

      if (idx_cur.lt.0) then
        write(lulog,*) '"',trim(label_in),'"'
        write(lulog,*) idx_cur
        call quit(1,'copy_mel','label not on list')
      end if

      me_current => op_info%mel_arr(idx_cur)%mel
      if (.not.associated(me_current%fhand))
     &     call quit(1,'copy_mel','no file handle defined for '//
     &                  trim(me_current%label))
      open_close_cur = me_current%fhand%unit.le.0

      if(open_close_res)then
        call file_open(me_res%fhand)
      endif

      call zeroop(me_res)

      nblk = me_res%op%n_occ_cls
      njoined = me_res%op%njoined 


      if(open_close_cur)then
        call file_open(me_current%fhand)
      endif

      do iblk = 1, nblk
        iblkoff = (iblk-1)*njoined
        occ => me_res%op%ihpvca_occ(1:ngastp,1:2,
     &                               iblkoff+1:iblkoff+njoined)

        iblk_cur = iblk_occ(occ,.true.,me_current%op,
     &                        me_res%op%blk_version(iblk))
        
        
        if (iblk_cur.lt.1) then
          call wrt_occ_n(lulog,occ,njoined)
          call quit(1,'copy_mel',
     &         'block not found: '//trim(me_current%op%name))
        end if

        if (trnsps) then

          call add_opblk_transp(xnorm2,1,fac,me_current,me_res,
     &         .false.,.true.,   
     &         iblk_cur,iblk,op_info,str_info,orb_info,copy_only)

        else

          call add_opblk(xnorm2,1,fac,me_current,me_res,
     &         iblk_cur,iblk,orb_info,copy_only)

        endif
      end do

      if (open_close_res)
     &     call file_close_keep(me_res%fhand)

      if (open_close_cur)
     &     call file_close_keep(me_current%fhand)

      if (ntest.ge.100) then
        write(lulog,*) 'dump of result list:'
        call wrt_mel_file(lulog,5,
     &       me_current,
     &       1,me_res%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'dump of result list:'
        call wrt_mel_file(lulog,5,
     &       me_res,
     &       1,me_res%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      return
      end
