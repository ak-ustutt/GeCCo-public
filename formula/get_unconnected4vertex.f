*----------------------------------------------------------------------*
      subroutine get_unconnected4vertex(iocc,ivtx,contr,op_info)
*----------------------------------------------------------------------*
*     return (in occupation form) the number of unconnected indices 
*     of vertex #ivtx in contraction contr
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'

      integer, intent(out) ::
     &     iocc(ngastp,2)
      integer, intent(in) ::
     &     ivtx
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     dag
      integer ::
     &     idx_op, iblk_op, iarc
      type(operator_array), pointer ::
     &     op_arr(:)
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(cntr_arc), pointer ::
     &     arc(:)

      ! get vertex occupation
      vertex => contr%vertex
      op_arr => op_info%op_arr
      idx_op = vertex(ivtx)%idx_op
      iblk_op = vertex(ivtx)%iblk_op
      iocc = op_arr(idx_op)%op%ihpvca_occ(1:ngastp,1:2,iblk_op)
      dag = op_arr(idx_op)%op%dagger
      if (dag) iocc = iocc_dagger(iocc)
      
      ! loop over all arcs and remove contractions
      arc => contr%arc
      do iarc = 1, contr%narc
        ! for proto-contractions: ignore certain arcs
        if (arc(iarc)%occ_cnt(1,1).lt.0) cycle
        if (arc(iarc)%link(1).eq.ivtx) then
          iocc = iocc - arc(iarc)%occ_cnt
        else if (arc(iarc)%link(2).eq.ivtx) then
          iocc = iocc - iocc_dagger(arc(iarc)%occ_cnt)
        end if
      end do

      return
      end
