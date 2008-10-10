*----------------------------------------------------------------------*
      integer function mem_free_words()
*----------------------------------------------------------------------*
*     return number of free words
*----------------------------------------------------------------------*
      use memman
      implicit none
      
      mem_free_words = mem_free

      return
      end
