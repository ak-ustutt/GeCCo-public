      subroutine clean_word_list(wlist)

      implicit none

      include 'def_word_list.h'

      type(word_list), intent(inout) ::
     &     wlist

      type(word_list_entry), pointer ::
     &     wl_pnt, wl_next, wl_up
 
      logical ::
     &     just_went_up

      ! not active?
      if (.not.associated(wlist%head)) return

      wl_pnt => wlist%head

      just_went_up = .false.
      ! deallocate all entries
      do 
        ! follow list one level down, if present 
        !  (and if we not just returned from there)
        if (associated(wl_pnt%down).and..not.just_went_up) then
          wl_pnt => wl_pnt%down
          cycle
        end if
        just_went_up = .false.
        ! remember pointers for next or level up entry
        wl_next => wl_pnt%next
        wl_up   => wl_pnt%up
        ! deallocate present entry
        deallocate(wl_pnt)
        ! next entry present?
        if (.not.associated(wl_next)) then
          ! level up present? if no, we are done ...
          if (.not.associated(wl_up)) exit
          ! else we go up ...
          wl_pnt => wl_up
          just_went_up = .true.
          cycle
        end if
        ! ... or go to next entry
        wl_pnt => wl_next
      end do

      ! reset all pointers
      wlist%head => null()
      wlist%tail => null()
      wlist%current => null()

      end
