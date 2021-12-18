*----------------------------------------------------------------------*
      subroutine contr_add_arc(contr,vtx1,vtx2,occ)
*----------------------------------------------------------------------*
*     add arc to contraction
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     vtx1, vtx2, occ(ngastp,2)
      type(contraction), intent(inout) ::
     &     contr

      call resize_contr(contr,0,contr%narc+1,0,0)

      contr%narc = contr%narc+1

      contr%arc(contr%narc)%link(1) = vtx1
      contr%arc(contr%narc)%link(2) = vtx2
      contr%arc(contr%narc)%occ_cnt(1:ngastp,1:2) = occ(1:ngastp,1:2)

      end
