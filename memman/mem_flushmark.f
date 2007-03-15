      integer function mem_flushmark(name)
      use memman
      implicit none
      
      character, intent(in), optional ::
     &     name*(*)

      if (present(name)) then
        mem_flushmark = memman_remsection(name)
      else
        mem_flushmark = memman_remsection()
      end if

      return
      end
