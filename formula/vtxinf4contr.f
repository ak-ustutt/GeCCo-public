*----------------------------------------------------------------------*
      subroutine vtxinf4contr(irestr_vtx,info_vtx,contr,op_info,ngas)
*----------------------------------------------------------------------*
*     set restriction for result and vertices of contr
*     and set info_vtx(2,nvtx+1): (total_spin,total_sym)
*     first enty is result, 2..nvtx+1 are the vertices
*     for occupations use occvtx4contr()
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     ngas
      integer, intent(out) ::
     &     irestr_vtx(2,ngas,2,2,contr%nvtx+1),
     &     info_vtx(2,contr%nvtx+1)

      integer ::
     &     nvtx, idx, idxop, iblkop
      type(cntr_vtx), pointer ::
     &     vertex(:)
      integer, pointer ::
     &     irestr_occ(:,:,:,:,:)
      type(operator), pointer ::
     &     op
      
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
          irestr_vtx(1:2,1:ngas,1:2,1:2,idx) = 0
          info_vtx(1,idx) = 0
          info_vtx(2,idx) = 1
        else
          op => op_info%op_arr(idxop)%op
          irestr_occ => op%igasca_restr
          if (.not.op%dagger) then
            irestr_vtx(1:2,1:ngas,1:2,1:2,idx) =
     &           irestr_occ(1:2,1:ngas,1:2,1:2,iblkop)
            info_vtx(1,idx) = op%mst
            info_vtx(2,idx) = op%gamt
          else
            irestr_vtx(1:2,1:ngas,1,1:2,idx) =
     &           irestr_occ(1:2,1:ngas,2,1:2,iblkop)
            irestr_vtx(1:2,1:ngas,2,1:2,idx) =
     &           irestr_occ(1:2,1:ngas,1,1:2,iblkop)
            info_vtx(1,idx) = -op%mst
            info_vtx(2,idx) = op%gamt
          end if
        end if
      end do

      return
      end
