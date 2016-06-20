*----------------------------------------------------------------------*
!>     wrapper routine for inp_show
*----------------------------------------------------------------------*
      subroutine show_input(lulog)
*----------------------------------------------------------------------*

      use parse_input,only: inp_show
      implicit none

      integer, intent(in) ::
     &     lulog

      call inp_show(lulog)
      return
      end
