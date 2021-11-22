      subroutine add_op(label_res,fac,label_sum,nsum,
     &     op_info,orb_info,str_info,replace)
*----------------------------------------------------------------------*
*     add nsum operator lists block by block
*     the block structure of operator "label_res" is decisive for the
*     blocks considered; they should be contained in all summands
*     replace = true : only replace blocks by existing blks of op. list
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000
      
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      integer, intent(in) ::
     &     nsum
      real(8), intent(in) ::
     &     fac(nsum)
      character(*), intent(in) ::
     &     label_res, label_sum(nsum)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      logical, intent(in) ::
     &     replace

      type(me_list), pointer ::
     &     me_res, me_current

      integer ::
     &     idx_res, idx_sum(nsum), isum,
     &     njoined, iblkoff, nblk, iblk, iblk_sum
      logical ::
     &     open_close_res, open_close_sum
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0, xnorm2
      integer, pointer  ::
     &     occ(:,:,:)

      integer, external ::
     &     idx_mel_list, iblk_occ

      call atim_csw(cpu0,sys0,wall0)

      if(ntest.ge.100)then
        write(lulog,*) '===================='
        write(lulog,*) ' Add up operators   '
        write(lulog,*) '===================='
        write(lulog,*) 'Result: ',trim(label_res)
        write(lulog,*) 'The factors and summands: '
        do isum = 1, nsum
          write(lulog,'(3x,f12.6,x,a)') fac(isum),
     &         trim(label_sum(isum))
        end do
      endif

      idx_res = idx_mel_list(label_res,op_info)

      if (idx_res.lt.0) then
        write(lulog,*) '"',trim(label_res),'"'
        write(lulog,*) idx_res
        call quit(1,'add_op','label not on list')
      end if

      ! Point to the relevant operators and their associated files.
      me_res => op_info%mel_arr(idx_res)%mel
      if (.not.associated(me_res%fhand))
     &     call quit(1,'add_op','no file handle defined for '//
     &                  trim(me_res%label))
      open_close_res = me_res%fhand%unit.le.0
     &    .and..not.(me_res%op%n_occ_cls.eq.1.and.me_res%fhand%buffered)

      do isum = 1, nsum
        idx_sum(isum) = idx_mel_list(label_sum(isum),op_info)
        if (idx_sum(isum).lt.0) then
          write(lulog,*) '"',trim(label_sum(isum)),'"'
          write(lulog,*) idx_sum
          call quit(1,'add_op','label not on list')
        end if
        me_current => op_info%mel_arr(idx_sum(isum))%mel
        if (.not.associated(me_current%fhand))
     &       call quit(1,'add_op','no file handle defined for '//
     &                  trim(me_current%label))
      end do

      if(open_close_res)then
        call file_open(me_res%fhand)
      endif

      ! zero result file
      if (.not.replace) call zeroop(me_res)

      nblk    = me_res%op%n_occ_cls
      njoined = me_res%op%njoined
      
      ! loop over summands
      do isum = 1, nsum
        me_current => op_info%mel_arr(idx_sum(isum))%mel
        if (me_current%op%njoined.ne.njoined)
     &       call quit(1,'add_op',
     &       'incompatible njoined: '//trim(me_current%op%name))
        open_close_sum = me_current%fhand%unit.le.0
     &        .and..not.(me_current%op%n_occ_cls.eq.1.and.
     &                   me_current%fhand%buffered)

        if(open_close_sum)then
          call file_open(me_current%fhand)
        endif

        ! loop over result blocks
        do iblk = 1, nblk
          iblkoff = (iblk-1)*njoined
          occ => me_res%op%ihpvca_occ(1:ngastp,1:2,
     &                               iblkoff+1:iblkoff+njoined)

          iblk_sum = iblk_occ(occ,.false.,me_current%op,
     &                        me_res%op%blk_version(iblk))

          if (iblk_sum.lt.1) then
            if (replace) cycle
            call wrt_occ_n(lulog,occ,njoined)
            call quit(1,'add_op',
     &           'block not found: '//trim(me_current%op%name))
          end if

          call add_opblk(xnorm2,1,fac(isum),me_current,me_res,
     &         iblk_sum,iblk,orb_info,replace)

        end do

        if (open_close_sum)
     &     call file_close_keep(me_current%fhand)

      end do

      if (open_close_res)
     &     call file_close_keep(me_res%fhand)

      if (ntest.ge.1000) then
        write(lulog,*) 'dump of result list:'
        call wrt_mel_file(lulog,5,
     &       me_res,
     &       1,me_res%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      call atim_csw(cpu,sys,wall)
      call prtim(lulog,'time for adding',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
