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

      ! defaults
      solve_tbar = .false.
      densities = 0

      icnt = is_keyword_set('calculate.CC_solve_tbar')
c dbg
      print *,'icnt = ',icnt
c dbg
      solve_tbar = solve_tbar.or.icnt.gt.0
c      icnt = is_keyword_set('calculate.CC_solve_lambda')
c      solve_tbar = icnt.gt.0
      icnt = is_keyword_set('calculate.properties')
      solve_tbar = solve_tbar.or.icnt.gt.0
      if (icnt.gt.0) densities = 1

      call get_argument_value('calculate.routes','simtraf',
     &     ival=ccsimtrf)

      return
      end

