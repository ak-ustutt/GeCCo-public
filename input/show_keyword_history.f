*----------------------------------------------------------------------*
      subroutine show_keyword_history(lulog)
*----------------------------------------------------------------------*
*     wrapper routine
*----------------------------------------------------------------------*

      use parse_input
      implicit none

      integer, intent(in) ::
     &     lulog

      call keyword_list(lulog,keyword_history,show_args=.false.)

      return
      end
