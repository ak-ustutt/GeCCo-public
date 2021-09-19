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
      contr%mxxarc = 0
      contr%mxfac = 0
      contr%nvtx = 0
      contr%narc = 0
      contr%nxarc = 0
      contr%nfac = 0
      contr%dagger = .false.
      contr%unique_set = .false.
      nullify(contr%vertex)
      nullify(contr%arc)
      nullify(contr%xarc)
      nullify(contr%joined)
      nullify(contr%svertex)
      nullify(contr%inffac)
      nullify(contr%vtx)
      nullify(contr%topo)
      nullify(contr%xlines)
      contr%index_info = .false.
      nullify(contr%contr_string)
      nullify(contr%result_string)
      contr%nidx = 0
      contr%nxidx = 0
      contr%total_sign = 0
      contr%eqvl_fact = 1d0
      
      return
      end
