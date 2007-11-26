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
     &     icnt

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
        call get_argument_value('method.R12','approx',str=r12_apprx)
      end if

      return
      end
