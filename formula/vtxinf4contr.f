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
      include 'stdunit.h'
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
     &     nvtx, idx, idxmel, idxop, iblkop, njoined
      type(cntr_vtx), pointer ::
     &     vertex(:)
      integer, pointer ::
     &     irestr_occ(:,:,:,:,:,:)
      type(operator), pointer ::
     &     op
      type(me_list), pointer ::
     &     mel
      
      idxop = contr%idx_res
      njoined = op_info%op_arr(idxop)%op%njoined

      nvtx = contr%nvtx
      vertex => contr%vertex
      do idx = 1, nvtx+njoined
        if (idx.le.njoined) then
          idxop = contr%idx_res 
          iblkop = (contr%iblk_res-1)*njoined+idx
        else
          idxop = vertex(idx-njoined)%idx_op
          ! is already set to compound index in case of super-vertices:
          iblkop = vertex(idx-njoined)%iblk_op
        end if

        if (idxop.gt.0) then
          idxmel = op_info%op2list(idxop)
          if (idxmel.le.0) then
            call prt_contr2(luout,contr,op_info)
            write(luout,*) op_info%op2list(1:op_info%nops)
            call quit(1,'vtxinf4contr',
     &         'No list associated with operator '//
     &         trim(op_info%op_arr(idxop)%op%name))
          end if
        end if

        if (idxop.eq.0) then
          irestr_vtx(1:2,1:ngas,1:2,1:2,idx) = 0
          info_vtx(1,idx) = 0
          info_vtx(2,idx) = 1
        else
          op => op_info%op_arr(idxop)%op
          mel => op_info%mel_arr(idxmel)%mel
          irestr_occ => op%igasca_restr
          if (.not.op%dagger) then
            irestr_vtx(1:2,1:ngas,1:2,1:2,idx) =
     &           irestr_occ(1:2,1:ngas,1:2,1:2,1,iblkop)
            info_vtx(1,idx) = mel%mst
            info_vtx(2,idx) = mel%gamt
          else
            irestr_vtx(1:2,1:ngas,1,1:2,idx) =
     &           irestr_occ(1:2,1:ngas,2,1:2,1,iblkop)
            irestr_vtx(1:2,1:ngas,2,1:2,idx) =
     &           irestr_occ(1:2,1:ngas,1,1:2,1,iblkop)
            info_vtx(1,idx) = -mel%mst
            info_vtx(2,idx) = mel%gamt
          end if
        end if
      end do

      return
      end
