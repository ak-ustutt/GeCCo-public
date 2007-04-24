*----------------------------------------------------------------------*
      subroutine mem_popmark()
*----------------------------------------------------------------------*
*     pop current memory section from internal stack
*----------------------------------------------------------------------*
      use memman
      implicit none
      
      call memman_section_stack(-1)

      return
      end
