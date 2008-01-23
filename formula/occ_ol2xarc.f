*----------------------------------------------------------------------*
      subroutine occ_ol2xarc(contr,occ_ol_vtx,svmap)
*----------------------------------------------------------------------*
*     given: remaining open lines per vertex (occ_ol_vtx)
*            a map telling us to which vertex of the result these
*            should contribute (svmap)
*     return: contr with xarc info telling the same
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      type(contraction), intent(inout) ::
     &     contr
      integer, intent(in) ::
     &     occ_ol_vtx(ngastp,2,contr%nvtx), svmap(contr%nvtx)

      integer ::
     &     nxarc, ixarc, nvtx, ivtx

      nvtx = contr%nvtx
      nxarc = 0
      do ivtx = 1, nvtx
        if (iocc_nonzero(occ_ol_vtx(1:ngastp,1:2,ivtx))) then
          nxarc = nxarc+1
        end if
      end do

      call resize_contr(contr,contr%nvtx,contr%narc,nxarc,0)
      contr%nxarc = nxarc

      ixarc = 0
      do ivtx = 1, nvtx
        if (iocc_nonzero(occ_ol_vtx(1:ngastp,1:2,ivtx))) then
          ixarc = ixarc+1
          if (svmap(ivtx).le.0)
     &         call quit(1,'occ_ol2xarc','inconsistent svmap')
          contr%xarc(ixarc)%link(1) = ivtx
          contr%xarc(ixarc)%link(2) = svmap(ivtx)
          contr%xarc(ixarc)%occ_cnt =
     &         occ_ol_vtx(1:ngastp,1:2,ivtx)

        end if
      end do

      return
      end
