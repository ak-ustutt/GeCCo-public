*----------------------------------------------------------------------*
      integer function mem_gotomark(name)
*----------------------------------------------------------------------*
*     go to a section with the name <name>
*----------------------------------------------------------------------*
      use memman
      implicit none
      
      character, intent(in) ::
     &     name*(*)

      mem_gotomark = memman_set_cursection(name)

      return
      end
