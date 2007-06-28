*----------------------------------------------------------------------*
      subroutine init_contr(contr)
*----------------------------------------------------------------------*
*     initialize all pointers and set all mxvtx, mxarc, mxfac variables
*     to zero (describing the current size of allocated space for 
*     vertices, arcs, and factorization info)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(inout) ::
     &     contr

      contr%mxvtx = 0
      contr%mxarc = 0
      contr%mxfac = 0
      contr%nvtx = 0
      contr%narc = 0
      contr%nfac = 0
      nullify(contr%vertex)
      nullify(contr%arc)
      nullify(contr%inffac)
      
      return
      end
