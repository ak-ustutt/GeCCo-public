*----------------------------------------------------------------------*
      subroutine show_keywords(lulog)
*----------------------------------------------------------------------*
*     wrapper routine
*----------------------------------------------------------------------*

      use parse_input
      implicit none

      integer, intent(in) ::
     &     lulog

      call keyword_list(lulog,keyword_root,show_args=.true.)

      return
      end
