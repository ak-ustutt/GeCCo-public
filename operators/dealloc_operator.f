*----------------------------------------------------------------------*
      subroutine dealloc_operator(op)
*----------------------------------------------------------------------*
*     deallocate all operator sub-arrays
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'ifc_memman.h'
      include 'par_globalmarks.h'

      type(operator), intent(inout) ::
     &     op

      integer ::
     &     nblk, iblk, ifree

      call mem_pushmark()

      ifree = mem_gotomark(operator_def)

      if (op%n_occ_cls.le.0.or.op%n_occ_cls.ge.1000) then
        write(luout,*) 'n_occ_cls = ',op%n_occ_cls
        call quit(1,'dealloc_operator',
     &              'suspicious number of blocks (bug?)')
      end if

      nblk = op%n_occ_cls

      if (associated(op%ihpvca_occ)) then
        deallocate(op%ihpvca_occ)
        op%ihpvca_occ => null()
      end if
      if (associated(op%ica_occ)) then
        deallocate(op%ica_occ)
        op%ica_occ => null()
      end if
      if (associated(op%igasca_restr)) then
        deallocate(op%igasca_restr)
        op%igasca_restr => null()
      end if
      if (associated(op%formal_blk)) then
        deallocate(op%formal_blk)
        op%formal_blk => null()
        ! if no-one has cheated, we should end up here when everything
        ! from this section has been removed cleanly
        ! could be made safer
        ifree = mem_dealloc(trim(op%name))
      end if

      call mem_popmark()

      return
      end
