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

      if (contr%mxvtx.gt.0) then
        deallocate(contr%vertex)
        deallocate(contr%svertex)
        deallocate(contr%joined)
      end if
      contr%vertex => null()
      contr%svertex => null()
      contr%joined => null()
      contr%mxvtx = 0
      if (contr%mxarc.gt.0) deallocate(contr%arc)
      contr%arc => null()
      contr%mxarc = 0
      if (contr%mxxarc.gt.0) deallocate(contr%xarc)
      contr%xarc => null()
      contr%mxxarc = 0
      if (contr%mxfac.gt.0) deallocate(contr%inffac)
      contr%inffac => null()
      contr%mxfac = 0
      if (contr%unique_set)
     &     deallocate(contr%vtx,contr%topo,contr%xlines)
      contr%vtx => null()
      contr%topo => null()
      contr%xlines => null()
      if (contr%index_info)
     &     deallocate(contr%contr_string)
      if (associated(contr%result_string))
     &     deallocate(contr%result_string)
      contr%contr_string => null()
      contr%result_string => null()

      return
      end
