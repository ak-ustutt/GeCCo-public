*----------------------------------------------------------------------*
*     wrapper routine for reg_show
*----------------------------------------------------------------------*
      subroutine show_keywords(lulog)
*----------------------------------------------------------------------*

      use parse_input2,only: reg_show
      implicit none

      integer, intent(in) ::
     &     lulog

      call reg_show(lulog)
      return
      end
