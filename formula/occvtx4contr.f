*----------------------------------------------------------------------*
      subroutine occvtx4contr(occ_vtx,contr,op_info)
*----------------------------------------------------------------------*
*     set occupations of operator blocks corresponding to result
*     and vertices of contraction contr
*     first occupation is result, 2..nvtx+1 are the vertices
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(out) ::
     &     occ_vtx(ngastp,2,contr%nvtx+1)

      integer ::
     &     nvtx, idx, idxop, iblkop
      type(cntr_vtx), pointer ::
     &     vertex(:)
      integer, pointer ::
     &     op_occ(:,:,:)

      nvtx = contr%nvtx
      vertex => contr%vertex
      do idx = 1, nvtx+1
        if (idx.eq.1) then
          idxop = contr%idx_res 
          iblkop = contr%iblk_res 
        else
          idxop = vertex(idx-1)%idx_op
          iblkop = vertex(idx-1)%iblk_op
        end if
        if (idxop.eq.0) then
          occ_vtx(1:ngastp,1:2,idx) = 0
        else
          op_occ => op_info%op_arr(idxop)%op%ihpvca_occ
          occ_vtx(1:ngastp,1:2,idx) = op_occ(1:ngastp,1:2,iblkop)
        end if
      end do

      return
      end
