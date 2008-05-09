*----------------------------------------------------------------------*
      subroutine occvtx_from_arcs(mode,occ_vtx,contr,nj_res)
*----------------------------------------------------------------------*
*     set up occ_vtx array from arc and xarc info on contr
*     mode = 0:  with result info (on first nj_res occ's)
*     mode = 1:  vertex occ's only
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     nj_res, mode
      integer, intent(out) ::
     &     occ_vtx(ngastp,2,*)

      integer ::
     &     idxres, nvtx, narc, nxarc, njoined_res, iarc, ivtx1, ivtx2,
     &     nend

      type(cntr_arc), pointer ::
     &     arc(:), xarc(:)

      nvtx   = contr%nvtx
      narc   = contr%narc
      nxarc   = contr%nxarc
      idxres = contr%idx_res

      arc => contr%arc
      xarc => contr%xarc

      ! use xarcs to reconstruct the uncontracted part
      nend = nvtx
      if (mode.eq.0) nend = nvtx+nj_res
      occ_vtx(1:ngastp,1:2,1:nend) = 0
      do iarc = 1, nxarc
        ivtx1 = xarc(iarc)%link(1)
        if (mode.eq.0) ivtx1 = ivtx1+nj_res
        ivtx2 = xarc(iarc)%link(2)
        occ_vtx(1:ngastp,1:2,ivtx1) =
     &       occ_vtx(1:ngastp,1:2,ivtx1) +
     &       xarc(iarc)%occ_cnt
        if (mode.eq.0) then
          occ_vtx(1:ngastp,1:2,ivtx2) =
     &         occ_vtx(1:ngastp,1:2,ivtx2) +
     &         xarc(iarc)%occ_cnt
        end if
      end do

      ! use arcs to reconstruct the contracted part
      do iarc = 1, narc
        ivtx1 = arc(iarc)%link(1)
        if (mode.eq.0) ivtx1 = ivtx1+nj_res
        ivtx2 = arc(iarc)%link(2)
        if (mode.eq.0) ivtx2 = ivtx2+nj_res
        occ_vtx(1:ngastp,1:2,ivtx1) =
     &       occ_vtx(1:ngastp,1:2,ivtx1) +
     &       arc(iarc)%occ_cnt
        occ_vtx(1:ngastp,1:2,ivtx2) =
     &       occ_vtx(1:ngastp,1:2,ivtx2) +
     &       iocc_dagger(arc(iarc)%occ_cnt)
      end do

      return
      end
