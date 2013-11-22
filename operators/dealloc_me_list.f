*----------------------------------------------------------------------*
      subroutine dealloc_me_list(mel)
*----------------------------------------------------------------------*
*     deallocate all me_list sub-arrays
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'
      include 'par_globalmarks.h'

      type(me_list), intent(inout) ::
     &     mel

      integer ::
     &     nblk, iblk, ifree

      call mem_pushmark()

      ifree = mem_gotomark(me_list_def)

      if (mel%op%n_occ_cls.le.0.or.mel%op%n_occ_cls.ge.1000) then
        write(lulog,*) 'n_occ_cls = ',mel%op%n_occ_cls
        call quit(1,'dealloc_me_list',
     &              'suspicious number of blocks (bug?)')
      end if

      ! no operator associated? nothing to do for us ...
      if (.not.associated(mel%op))
     &     return

      nblk = mel%op%n_occ_cls

      if (associated(mel%off_op_occ)) then
        deallocate(mel%off_op_occ)
        mel%off_op_occ => null()
      end if
      if (associated(mel%len_op_occ)) then
        deallocate(mel%len_op_occ)
        mel%len_op_occ => null()
      end if
      if (associated(mel%idx_graph)) then
        deallocate(mel%idx_graph)
        mel%idx_graph => null()
      end if
      if (associated(mel%ld_op_gmox)) then
        do iblk = 1, nblk
          deallocate(mel%ld_op_gmox(iblk)%d_gam_ms)
        end do
        deallocate(mel%ld_op_gmox)
        mel%ld_op_gmox => null()
      end if
      if (associated(mel%off_op_gmo)) then
        do iblk = 1, nblk
          deallocate(mel%off_op_gmo(iblk)%gam_ms)
        end do
        deallocate(mel%off_op_gmo)
        mel%off_op_gmo => null()
      end if
      if (associated(mel%len_op_gmo)) then
        do iblk = 1, nblk
          deallocate(mel%len_op_gmo(iblk)%gam_ms)
        end do
        deallocate(mel%len_op_gmo)
        mel%len_op_gmo => null()
      end if
      if (associated(mel%off_op_gmox)) then
        do iblk = 1, nblk
          deallocate(mel%off_op_gmox(iblk)%d_gam_ms)
          deallocate(mel%off_op_gmox(iblk)%did)
          deallocate(mel%off_op_gmox(iblk)%ndis)
        end do
        deallocate(mel%off_op_gmox)
        mel%off_op_gmox => null()
      end if
      if (associated(mel%len_op_gmox)) then
        do iblk = 1, nblk
          deallocate(mel%len_op_gmox(iblk)%d_gam_ms)
        end do
        deallocate(mel%len_op_gmox)
        mel%len_op_gmox => null()
        ! if noone has cheated, we should end up here
        ! could be made safer
        ifree = mem_dealloc(trim(mel%label)//'-1')
        ifree = mem_dealloc(trim(mel%label)//'-2')
      end if

      call mem_popmark()

      if (associated(mel%fhand)) then
        call detach_mel_file(mel,.true.)
      end if

      return
      end
