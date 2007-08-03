*----------------------------------------------------------------------*
      subroutine cc_input()
*----------------------------------------------------------------------*
*     find out which kind of CC calculation is to run and set switches
*     in common /cc_route/
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'
      include 'cc_routes.h'

      integer ::
     &     icnt

      do_cc = .false.

      icnt = is_keyword_set('method.CC')

      if (icnt.gt.0) then

        do_cc = .true.

        ! defaults
        solve_tbar = .false.
        densities = 0

        icnt = is_keyword_set('calculate.CC_solve_tbar')
        solve_tbar = solve_tbar.or.icnt.gt.0
        icnt = is_keyword_set('calculate.properties')
        solve_tbar = solve_tbar.or.icnt.gt.0
        if (icnt.gt.0) densities = 1

        call get_argument_value('calculate.routes','simtraf',
     &       ival=ccsimtrf)

      end if

      return
      end

