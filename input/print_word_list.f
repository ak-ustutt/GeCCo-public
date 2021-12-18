      subroutine print_word_list(lulog,wlist)

      implicit none

      include 'def_word_list.h'

      integer, intent(in) ::
     &     lulog
      type(word_list), intent(inout) ::
     &     wlist

      type(word_list_entry), pointer ::
     &     wl_pnt
      integer ::
     &     level
      logical ::
     &     just_moved_up

      ! not active?
      if (.not.associated(wlist%head)) then
        write(lulog,*) 'word list is empty'
        return
      end if

      wl_pnt => wlist%head

      level = 0
      just_moved_up = .false.
      ! loop over all entries
      do 

        if (.not.just_moved_up) then
          write(lulog,'(1x,i3,4a)') level,'"',trim(wl_pnt%word),
     &                   '", sep=',wl_pnt%sep
          write(lulog,'(9x,a,i4,a,i4)')               
     &                   ' line=',wl_pnt%line,
     &                   ' col=',wl_pnt%col
        else
          write(lulog,*) 'returned to level ',level
        end if

        ! level down?
        if (.not.just_moved_up.and.associated(wl_pnt%down)) then
          level = level+1
          wl_pnt => wl_pnt%down
        else if (associated(wl_pnt%next)) then
          just_moved_up = .false.
          wl_pnt => wl_pnt%next
        else if (associated(wl_pnt%up)) then
          level = level-1
          wl_pnt => wl_pnt%up
          just_moved_up = .true.
        else
          exit
        end if
      end do

      if (level/=0) then
        write(lulog,*) 'warning: level/=0 (level = ',level,')'
      end if

      end
