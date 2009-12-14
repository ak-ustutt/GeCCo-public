      subroutine new_word_list_entry(wlist)
      !
      ! insert new entry at position "current"
      ! (so we can also insert nodes)
      !

      implicit none

      include 'def_word_list.h'

      type(word_list), intent(inout) ::
     &     wlist

      type(word_list_entry), pointer ::
     &     wl_pnt

      if (.not.associated(wlist%head)) then
        ! not active? make first entry then

        allocate(wlist%head)
        wlist%tail => wlist%head
        wlist%current => wlist%head
        wl_pnt => wlist%current
        wl_pnt%next => null()
      
      else if (.not.associated(wlist%current%next)) then
        ! on last entry: extend
        
        wl_pnt => wlist%current
        allocate(wl_pnt%next)
        wl_pnt => wl_pnt%next
        wl_pnt%next => null()
        wlist%current => wl_pnt
        wlist%tail => wl_pnt

      else

        wl_pnt => wlist%current%next
        wlist%current%next => null()
        allocate(wlist%current%next)
        wlist%current%next%next => wl_pnt
        wlist%current => wlist%current%next
        wl_pnt => wlist%current

      end if

      wl_pnt%word(1:maxlen_word) = ' '
      wl_pnt%sep = ' '

      end
