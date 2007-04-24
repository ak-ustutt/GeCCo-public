*----------------------------------------------------------------------*
      integer function mem_flushmark(name)
*----------------------------------------------------------------------*
*     flush last section, or section with name <name>
*     all memory allocated in that section is freed
*----------------------------------------------------------------------*
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
