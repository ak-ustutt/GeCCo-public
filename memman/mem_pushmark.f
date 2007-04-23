*----------------------------------------------------------------------*
      subroutine mem_pushmark()
*----------------------------------------------------------------------*
*     push current memory section to internal stack
*----------------------------------------------------------------------*
      use memman
      implicit none
      
      call memman_section_stack(+1)

      return
      end
