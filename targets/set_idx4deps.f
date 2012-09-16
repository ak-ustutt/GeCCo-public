*----------------------------------------------------------------------*
      subroutine set_idx4deps(tgt_info)
*----------------------------------------------------------------------*
*     scan through target list and set for all dependencies the actual
*     target indices
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'

      type(target_info), intent(inout) ::
     &     tgt_info

      type(target), pointer ::
     &     tgt
      integer ::
     &     ierr, itgt, jtgt, idx
      integer, external ::
     &     idx_target

      ierr = 0
      do itgt = 1, tgt_info%ntargets
        tgt => tgt_info%array(itgt)%tgt
        if (tgt%n_joined_with.gt.0) then
          if (associated(tgt%idx_joined_with))
     &         deallocate(tgt%idx_joined_with)
          allocate(tgt%idx_joined_with(tgt%n_joined_with))
          do idx = 1, tgt%n_joined_with
            jtgt = idx_target(tgt%joined_with(idx),tgt_info)
            if (jtgt.lt.1) then
              ierr = ierr+1
              write(luout,'(x,a," - ",a)')
     &             trim(tgt%name),trim(tgt%joined_with(idx))
            end if
            tgt%idx_joined_with(idx) = jtgt
          end do
        end if
        if (tgt%n_depends_on.gt.0) then
          if (associated(tgt%idx_depends_on))
     &         deallocate(tgt%idx_depends_on)
          allocate(tgt%idx_depends_on(tgt%n_depends_on))
          do idx = 1, tgt%n_depends_on
            jtgt = idx_target(tgt%depends_on(idx),tgt_info)
            if (jtgt.lt.1) then
              ierr = ierr+1
              write(luout,'(x,a," : """,a,"""")')
     &             trim(tgt%name),trim(tgt%depends_on(idx))
            end if
            tgt%idx_depends_on(idx) = jtgt
          end do
        end if

      end do

      if (ierr.gt.0) then
        call quit(1,'set_idx4deps','undefined dependencies (see above)')
      end if

      return
      end
