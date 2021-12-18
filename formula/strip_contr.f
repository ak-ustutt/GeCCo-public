*----------------------------------------------------------------------*
      subroutine strip_contr(contr,mode)
*----------------------------------------------------------------------*
*     mode = 1:  remove (x 0 [C]) and (0 x [C]) proto-arcs
*     mode = .... to come
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      type(contraction), intent(inout) ::
     &     contr
      integer, intent(in) ::
     &     mode

      integer ::
     &     narc, narc_new, iarc
      type(cntr_arc), pointer ::
     &     arc(:)

      narc = contr%narc
      arc => contr%arc
      
      select case (mode)
      case(1)
        narc_new = 0
        do iarc = 1, narc
          if (arc(iarc)%link(1).eq.0.or.
     &         arc(iarc)%link(2).eq.0) cycle
          narc_new = narc_new + 1
          if (narc_new.lt.iarc)
     &         arc(narc_new) = arc(iarc)
        end do
        contr%narc = narc_new
      case default
        call quit(1,'strip_contr','undefined mode')
      end select

      return
      end
