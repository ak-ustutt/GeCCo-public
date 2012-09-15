      logical function advance_word_list_entry(wlist,direction)
      !
      ! advance the "current" pointer in direction:
      !    ' '    same level
      !    'u'    level up
      !    'd'    level down
      ! return .false., if this is not possible
      !

      implicit none

      include 'def_word_list.h'

      type(word_list), intent(inout) ::
     &     wlist
      character, intent(in) ::
     &     direction

      type(word_list_entry), pointer ::
     &     wl_pnt

      advance_word_list_entry = .false.
      if (.not.associated(wlist%head)) return

      select case (direction)
      case(' ')
        if (associated(wlist%current%next)) then
          wlist%current => wlist%current%next
          advance_word_list_entry = .true.
        end if   
      case('d')
        if (associated(wlist%current%down)) then
          wlist%current => wlist%current%down
          advance_word_list_entry = .true.
        end if   
      case('u')
        if (associated(wlist%current%up)) then
          wlist%current => wlist%current%up
          advance_word_list_entry = .true.
        end if   
      case default
         call quit(1,'advance_word_list','illegal direction "'//
     &               direction//'"')
      end select

      return

      end
