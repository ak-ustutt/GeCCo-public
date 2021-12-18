*----------------------------------------------------------------------*
      integer function mem_setmark(name)
*----------------------------------------------------------------------*
*     create a new section
*----------------------------------------------------------------------*
      use memman
      implicit none
      
      character, intent(in) ::
     &     name*(*)

      mem_setmark = memman_addsection(name)

      return
      end
