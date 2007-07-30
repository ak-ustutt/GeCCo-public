*----------------------------------------------------------------------*
      subroutine dealloc_contr(contr)
*----------------------------------------------------------------------*
*     deallocate all subfields in contraction contr
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(inout) ::
     &     contr

      if (contr%mxvtx.gt.0) deallocate(contr%vertex)
      contr%vertex => null()
      contr%mxvtx = 0
      if (contr%mxarc.gt.0) deallocate(contr%arc)
      contr%arc => null()
      contr%mxarc = 0
      if (contr%mxfac.gt.0) deallocate(contr%inffac)
      contr%inffac => null()
      contr%mxfac = 0

      return
      end
