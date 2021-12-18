*----------------------------------------------------------------------*
      subroutine occvtx4contr(mode,occ_vtx,contr,op_info)
*----------------------------------------------------------------------*
*     set occupations of operator blocks corresponding to result
*     and vertices of contraction contr
*     mode = 0:
*     first njoined occupations are result, njoined+1..njoined+nvtx 
*     are the vertices (njoined > 1, if result is a super-vertex)
*     mode = 1:
*     skip info for result vertices
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
     &     mode
      integer, intent(out) ::
     &     occ_vtx(ngastp,2,contr%nvtx+1) !CAUTION: correct last dimension
                                          !should be contr%nvtx+njoined !
      logical ::
     &     dagger
      integer ::
     &     nvtx, idx, idxop, iblkop, njoined
      type(cntr_vtx), pointer ::
     &     vertex(:)
      integer, pointer ::
     &     op_occ(:,:,:)

      if (mode.ne.0.and.mode.ne.1)
     &     call quit(1,'occvtx4contr','unknown mode')

      if (mode.eq.0) then
        idxop = contr%idx_res
        njoined = op_info%op_arr(idxop)%op%njoined
      else
        njoined = 0
      end if
      nvtx = contr%nvtx
      vertex => contr%vertex
      do idx = 1, nvtx+njoined
        if (idx.le.njoined) then
          idxop = contr%idx_res 
          dagger = contr%dagger
        else
          idxop = vertex(idx-njoined)%idx_op
          dagger = vertex(idx-njoined)%dagger
        end if
        if (idxop.eq.0) then
          occ_vtx(1:ngastp,1:2,idx) = 0
        else
          op_occ => op_info%op_arr(idxop)%op%ihpvca_occ
        end if
        if (idx.le.njoined) then
          if (.not.dagger) then
            iblkop = (contr%iblk_res-1)*njoined+idx         
          else
            iblkop = (contr%iblk_res)*njoined-idx+1         
          end if
        else
          iblkop = vertex(idx-njoined)%iblk_op 
        end if
        if (idxop.eq.0) then
          occ_vtx(1:ngastp,1:2,idx) = 0
        else
          if (.not.dagger) then
            occ_vtx(1:ngastp,1:2,idx) = op_occ(1:ngastp,1:2,iblkop)
          else
            occ_vtx(1:ngastp,1,idx) = op_occ(1:ngastp,2,iblkop)
            occ_vtx(1:ngastp,2,idx) = op_occ(1:ngastp,1,iblkop)
          end if
        end if
      end do

      return
      end
