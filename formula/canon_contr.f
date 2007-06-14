*----------------------------------------------------------------------*
      subroutine canon_contr(contr)
*----------------------------------------------------------------------*
*     canonicalize contraction
*     a) for all commuting vertices: (idx_op,iblk_op) in lexically 
*        ascending order (currently: not done)
*     b) for all all arcs: (left_op,right_op) in lexcically ascending
*        order (and we check that link(1) is left of link(2))
*     note: commuting vertices means that there is no contraction 
*     between the two vertices such that their sequence is arbitrary
*     commuting operators cannot be identified at this stage
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      type(contraction), intent(inout) ::
     &     contr

      integer ::
     &     narc, nvtx, idx, jdx, ivtx, ival, jval
      type(cntr_arc), pointer ::
     &     arc(:)
      type(cntr_arc) ::
     &     arc_sv

      integer, external ::
     &     int_pack

      narc = contr%narc
      nvtx = contr%nvtx
      
c      call quit(1,'canon_contr', 'baustelle')
c      ! loop over vertices, choosing the current one as root
c      ! for contractions with any of the ones that follow
c      do ivtx_root = 1, nvtx
c        ! set up list of vertices directly connecting to
c        ! current root
c        do jvtx = ivtx_root+1, nvtx
c          if (connec
c        ! test for canonical sequence
c        ! if OK, we go to next vertex as root
c        ! else:
c        ! analyze connectivity of effective vertices
c
c      end do
c
      ! re-sort arcs, if necessary
      ! first round: check that link(1) < link(2)
      arc => contr%arc
      do idx = 1, narc
        if (arc(idx)%link(1).gt.arc(idx)%link(2)) then
          ivtx = arc(idx)%link(1)
          arc(idx)%link(1) = arc(idx)%link(2)
          arc(idx)%link(2) = ivtx
          arc(idx)%occ_cnt = iocc_dagger(arc(idx)%occ_cnt)
        end if
      end do

      ! second round: insertion sort
      do idx = 2, narc
        ival = int_pack(arc(idx)%link,2,nvtx+1)
        arc_sv = arc(idx)
        jdx = idx-1
        do while(jdx.gt.0)
          jval = int_pack(arc(jdx)%link,2,nvtx+1)
          if (jval.le.ival) exit
          arc(jdx+1) = arc(idx)
          jdx = jdx-1
        end do
        arc(jdx+1) = arc_sv
      end do

      return
      end
