*----------------------------------------------------------------------*
      subroutine canon_contr(contr,reo,vtx_reo)
*----------------------------------------------------------------------*
*     canonicalize contraction
*     a) apply reordering suggested by topo_contr() to vertices
*     b) for all all arcs: (left_op,right_op) in lexcically ascending
*        order (and we check that link(1) is left of link(2))
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      type(contraction), intent(inout) ::
     &     contr
      logical, intent(in) ::
     &     reo
      integer, intent(in) ::
     &     vtx_reo(*)

      integer ::
     &     nxarc, narc, nvtx, idx, jdx, ivtx, ival, jval
      integer, pointer ::
     &     vtx_oer(:)
      type(cntr_arc), pointer ::
     &     arc(:), xarc(:)
      type(cntr_arc) ::
     &     arc_sv
      type(cntr_vtx), pointer ::
     &     vtx_new(:)
      integer, pointer ::
     &     svertex_new(:)

      integer, external ::
     &     int_pack

      narc = contr%narc
      nxarc = contr%nxarc
      nvtx = contr%nvtx

      ! reorder vertices, if necessary
      if (reo) then
        allocate(vtx_new(contr%mxvtx)) ! use mxvtx!
        do ivtx = 1, nvtx
          vtx_new(ivtx) = contr%vertex(vtx_reo(ivtx))
        end do
        deallocate(contr%vertex)
        contr%vertex => vtx_new
        ! update svertex array
        allocate(svertex_new(contr%mxvtx))
        do ivtx = 1, nvtx
          svertex_new(ivtx) = contr%svertex(vtx_reo(ivtx))
        end do
        deallocate(contr%svertex)
        contr%svertex => svertex_new
        ! revert reordering info to old->new        
        allocate(vtx_oer(nvtx))
        do ivtx = 1, nvtx
          vtx_oer(vtx_reo(ivtx)) = ivtx
        end do
        ! update vertex references in arcs:
        arc => contr%arc
        do idx = 1, narc
          arc(idx)%link(1) = vtx_oer(arc(idx)%link(1))
          arc(idx)%link(2) = vtx_oer(arc(idx)%link(2))
        end do
        ! update vertex references in xarcs:
        xarc => contr%xarc
        do idx = 1, nxarc
          xarc(idx)%link(1) = vtx_oer(xarc(idx)%link(1))
          ! no update of link(2)!
        end do
        deallocate(vtx_oer)
      end if

      ! re-sort arcs, if necessary
      ! first round: check that link(1) < link(2)
      ! and remove zero contractions
      arc => contr%arc
      jdx = 0
      do idx = 1, narc
        if (arc(idx)%link(1).gt.arc(idx)%link(2)) then
          ivtx = arc(idx)%link(1)
          arc(idx)%link(1) = arc(idx)%link(2)
          arc(idx)%link(2) = ivtx
          arc(idx)%occ_cnt = iocc_dagger(arc(idx)%occ_cnt)
        end if
        if (iocc_nonzero(arc(idx)%occ_cnt)) then
          jdx = jdx+1
          if (jdx.lt.idx) arc(jdx) = arc(idx)
        end if
      end do
      narc= jdx
      contr%narc = narc

      ! second round: insertion sort
      if (narc.gt.1)
     &     call arc_sort(arc,narc,nvtx)

      ! sort xarcs as well
      if (contr%nxarc.gt.1)
     &     call arc_sort(contr%xarc,contr%nxarc,nvtx)

      ! update svertex info
      call update_svtx4contr(contr)

      return
      end
