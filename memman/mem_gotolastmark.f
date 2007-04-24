*----------------------------------------------------------------------*
      integer function mem_gotolastmark()
*----------------------------------------------------------------------*
*     go to last memory section
*----------------------------------------------------------------------*
      use memman
      implicit none
      
      mem_gotolastmark = memman_set_cursection()

      return
      end
