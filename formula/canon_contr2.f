*----------------------------------------------------------------------*
      subroutine canon_contr2(contr,reo,vtx_reo)
*----------------------------------------------------------------------*
*     canonicalize contraction
*     a) apply reordering suggested by topo_contr() to vertices
*     b) for all all arcs: (left_op,right_op) in lexcically ascending
*        order (and we check that link(1) is left of link(2))
*     Debug version used for contractions where njoined=1 in all cases.
*     OBSOLETE
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
     &     narc, nvtx, idx, jdx, ivtx, ival, jval
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
      call quit(1,'canon_contr2','call to obsolete routine')

      narc = contr%narc
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
        ! update vertex references in arcs:
        ! revert reordering info to old->new        
        allocate(vtx_oer(nvtx))
        do ivtx = 1, nvtx
          vtx_oer(vtx_reo(ivtx)) = ivtx
        end do
        arc => contr%arc
        do idx = 1, narc
          arc(idx)%link(1) = vtx_oer(arc(idx)%link(1))
          arc(idx)%link(2) = vtx_oer(arc(idx)%link(2))
        end do
        ! Update the vertex references in xarcs.
        xarc => contr%xarc
        do idx = 1, contr%nxarc
          xarc(idx)%link(1) = vtx_oer(xarc(idx)%link(1))
          xarc(idx)%link(2) = xarc(idx)%link(1)
        enddo
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

c      do idx = 2, narc
c        ival = int_pack(arc(idx)%link,2,nvtx+1)
c        arc_sv = arc(idx)
c        jdx = idx-1
c        do while(jdx.gt.0)
c          jval = int_pack(arc(jdx)%link,2,nvtx+1)
c          if (jval.le.ival) exit
c          arc(jdx+1) = arc(jdx)
c          jdx = jdx-1
c        end do
c        arc(jdx+1) = arc_sv
c      end do

      ! sort xarcs as well
      if (contr%nxarc.gt.1)
     &     call arc_sort(contr%xarc,contr%nxarc,nvtx)

      ! update svertex info
      call update_svtx4contr(contr)

      return
      end
