*----------------------------------------------------------------------*
      subroutine r12_input()
*----------------------------------------------------------------------*
*     find out which kind of R12 calculation is to run and set switches
*     on common /explicit/
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'
      include 'explicit.h'

      integer ::
     &     icnt, lev

      ! necessary defaults
      explicit = .false.

      icnt = is_keyword_set('method.R12')

      if (icnt.gt.0) then
        explicit=.true.
        call get_argument_value('method.R12','ansatz',ival=ansatze)
        if(ansatze.gt.3.or.ansatze.lt.1)then
          call quit(1,'do_calc',
     &         'Undefined R12 ansatz requested.')
        endif
        call get_argument_value('method.R12','triples',ival=trir12)
        r12_apprx = ' '
        call get_argument_value('method.R12','approx',str=r12_apprx)
      end if

      icnt = is_keyword_set('method.MP')
      mp2 = .false.
      if(icnt.gt.0)then
        call get_argument_value('method.MP','level',ival=lev)
        if(lev.eq.2) mp2 = .true.
      endif  

      return
      end
