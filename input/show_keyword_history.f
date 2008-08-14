*----------------------------------------------------------------------*
      subroutine show_keyword_history(luout)
*----------------------------------------------------------------------*
*     wrapper routine
*----------------------------------------------------------------------*

      use parse_input
      implicit none

      integer, intent(in) ::
     &     luout

      call keyword_list(luout,keyword_history,show_args=.false.)

      return
      end
