*----------------------------------------------------------------------*
      subroutine set_primitive_contr(contr)
*----------------------------------------------------------------------*
*     assuming a contraction for which nvtx is set and all arrays
*     have sufficient size: set super vertex info for primitive case
*     (i.e. each vertex is an operator)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      type(contraction), intent(inout) ::
     &     contr

      integer ::
     &     ivtx

      if (contr%nvtx.lt.0.or.contr%nvtx.gt.1000) then
        write(lulog,*) ' nvtx = ', contr%nvtx
        call quit(1,'set_primitive_contr',
     &       'suspicious value of nvtx (bug)?')
      end if

      contr%nsupvtx = contr%nvtx
      do ivtx = 1, contr%nvtx
        contr%svertex(ivtx) = ivtx
        contr%joined(0,ivtx) = 1
        contr%joined(1,ivtx) = ivtx
        contr%joined(2:contr%nvtx,ivtx) = 0
      end do

      return
      end
