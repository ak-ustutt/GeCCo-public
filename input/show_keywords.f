*----------------------------------------------------------------------*
      subroutine show_keywords(luout)
*----------------------------------------------------------------------*
*     wrapper routine
*----------------------------------------------------------------------*

      use parse_input
      implicit none

      integer, intent(in) ::
     &     luout

      call keyword_list(luout,keyword_root,show_args=.true.)

      return
      end
