      subroutine new_word_list_entry(wlist,branch)
      !
      ! insert new entry at position "current"
      ! (so we can also insert nodes)
      !

      implicit none

      include 'def_word_list.h'

      type(word_list), intent(inout) ::
     &     wlist
      logical, intent(in) ::
     &     branch

      type(word_list_entry), pointer ::
     &     wl_pnt

      if (.not.associated(wlist%head)) then
        ! not active? make first entry then

        ! ignore "branch"

        allocate(wlist%head)
        wlist%tail => wlist%head
        wlist%current => wlist%head
        wl_pnt => wlist%current
        wl_pnt%next => null()
        wl_pnt%up   => null()
        wl_pnt%down => null()
      
      else if (.not.branch) then
        if (.not.associated(wlist%current%next)) then
          ! on last entry: extend
        
          wl_pnt => wlist%current
          allocate(wl_pnt%next)
          wl_pnt%next%up => wl_pnt%up
          wl_pnt => wl_pnt%next
          wl_pnt%next => null()
          wl_pnt%down => null()
          wlist%current => wl_pnt
          wlist%tail => wl_pnt

        else
          ! insert

          wl_pnt => wlist%current%next
          wlist%current%next => null()
          allocate(wlist%current%next)
          wlist%current%next%up => wlist%current%up
          wlist%current%next%down => null()
          wlist%current%next%next => wl_pnt
          wlist%current => wlist%current%next
          wl_pnt => wlist%current

        end if

      else ! branch
   
          wl_pnt => wlist%current
          allocate(wl_pnt%down)
          wl_pnt%down%up => wl_pnt
          wl_pnt => wl_pnt%down
          wl_pnt%next => null()
          wl_pnt%down => null()
          wlist%current => wl_pnt
          wlist%tail => wl_pnt

      end if

      wl_pnt%word(1:maxlen_word) = ' '
      wl_pnt%sep = ' '
      wl_pnt%line = -1
      wl_pnt%col  = -1

      end
