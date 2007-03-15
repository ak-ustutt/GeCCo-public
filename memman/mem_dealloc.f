      integer function mem_dealloc(name)
      use memman
      implicit none
      
      character, intent(in), optional ::
     &     name*(*)

      if (present(name)) then
        mem_dealloc = memman_dealloc(name)
      else
        mem_dealloc = memman_dealloc()
      end if

      return
      end
