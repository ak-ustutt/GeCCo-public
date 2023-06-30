      subroutine contr_op_replace(all_found,contr,
     &     idx_opin,dagin,idx_opout,dagout,
     &     op_info)
*-----------------------------------------------------------------------
*     isolated sub from op_replace: for a given contraction, replace
*              operators with index idx_opin
*-----------------------------------------------------------------------
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     idx_opin, idx_opout
      logical, intent(inout) ::
     &     all_found
      logical, intent(in) ::
     &     dagin, dagout
      type(contraction), target, intent(inout) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      type(operator), pointer ::
     &     opin_pnt, opout_pnt
      integer ::
     &     idxin, idxout, idx, idx_form_op, idx_form_blk, idx_blk_out,
     &     nvtx, narc, nxarc, nfac, njoined, idx_join, ij, jdx,
     &     isuper
      integer, allocatable ::
     &     occ_temp(:,:,:), vtx_chng_idx(:)
      integer, pointer ::
     &     ivtx_reo(:), occ_vtx(:,:,:), vtx_where(:)
      logical ::
     &     reo, dag_form_op, dagi_xor_dago
      logical, pointer ::
     &     fix_vtx(:)

      integer, external ::
     &     iblk_occ
      logical, external ::
     &     opblk_restr_cmp

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'contr_op_replace')
        write(lulog,*) 'idx_opin, idx_opout: ',idx_opin, idx_opout
      end if

      idxin = idx_opin
      idxout = idx_opout
        
      opin_pnt => op_info%op_arr(idxin)%op
      opout_pnt => op_info%op_arr(idxout)%op
      ! Allocate a temporary array.
      njoined = opout_pnt%njoined
      allocate(occ_temp(ngastp,2,njoined),vtx_chng_idx(njoined))

      dagi_xor_dago = dagin
      if (dagout) dagi_xor_dago=.not.dagi_xor_dago

      ! make sure that the operators match:
      if (opin_pnt%njoined.ne.njoined)
     &     call quit(1,'contr_op_replace',
     &     'the shape of the operators does not match: '//
     &     trim(opin_pnt%name)//' '//trim(opout_pnt%name))


      all_found = .true.

      ! Loop over the contraction's vertices.
      nvtx = contr%nvtx
      do idx = 1,nvtx
          ! Check the operator index of the vertex.
          idx_form_op = contr%vertex(idx)%idx_op
          dag_form_op = contr%vertex(idx)%dagger

          ! If the index of the operator vertex equals that of the
          ! intermediate operator...
          idx_blk_out = -1
          if (idx_form_op.eq.idxin.and.(dag_form_op.eqv.dagin)) then
              idx_form_blk = contr%vertex(idx)%iblk_op

              occ_temp(1:ngastp,1:2,1)=
     &             opin_pnt%ihpvca_occ(1:ngastp,1:2,idx_form_blk)

              ! if necessary, look for further vertices of super-vertex
              if (njoined.gt.1) then
                 isuper = contr%svertex(idx)
                 if (contr%joined(0,isuper).ne.njoined) then
                    write(lulog,*) 'njoined,joined(0):',
     &                   njoined,contr%joined(0,isuper)
                    call quit(1,'op_replace','inconsistency')
                 end if
                 vtx_where => contr%joined(1:njoined,isuper)
                 if (.not.dagin) then
                    do ij = 2, njoined
                      occ_temp(1:ngastp,1:2,ij) =
     &                     opin_pnt%ihpvca_occ(1:ngastp,1:2,
     &                                              idx_form_blk-1+ij)
                    end do
                    idx_form_blk = (idx_form_blk-1)/njoined + 1
                 else
                    do ij = 1, njoined
                      occ_temp(1:ngastp,1:2,ij) =
     &                     opin_pnt%ihpvca_occ(1:ngastp,1:2,
     &                                        idx_form_blk-njoined+ij)
                   end do
                   idx_form_blk = (idx_form_blk-njoined)/njoined + 1
                end if
             end if

             ! apply the change directly:
             ! Locate the formal block's counterpart in the actual 
             ! operator. 
             idx_blk_out =
     &               iblk_occ(occ_temp,dagi_xor_dago,opout_pnt,
     &                        opin_pnt%blk_version(idx_form_blk))

             ! new: also check the restrictions
             if (idx_blk_out.gt.0) then
                if(.not.opblk_restr_cmp(opin_pnt,idx_form_blk,dagin,
     &                                    opout_pnt,idx_blk_out,dagout))
     &                 idx_blk_out = -1
             end if

             if (idx_blk_out.gt.0) then
                if (njoined.eq.1) then
                   contr%vertex(idx)%idx_op = idxout
                   contr%vertex(idx)%iblk_op =
     &                   idx_blk_out
                   contr%vertex(idx)%dagger =
     &                   dagout
                else if (.not.dagout) then
                   do ij = 1, njoined
                      jdx = vtx_where(ij)
                      contr%vertex(jdx)%idx_op = idxout
                      contr%vertex(jdx)%iblk_op =
     &                     (idx_blk_out-1)*njoined+ij
                      contr%vertex(jdx)%dagger =
     &                     dagout
                   end do
                else
                   do ij = 1, njoined
                      jdx = vtx_where(ij)
                      contr%vertex(jdx)%idx_op = idxout
                      contr%vertex(jdx)%iblk_op =
     &                     (idx_blk_out)*njoined-ij+1
                      contr%vertex(jdx)%dagger =
     &                     dagout
                   end do
                end if

             else
                all_found = .false.
             end if
          end if

          if (idx_blk_out.gt.0) then
             ! Ensure everything is properly set up.
             narc = contr%narc
             nfac = contr%nfac
             nxarc = contr%nxarc
             call resize_contr(contr,nvtx,narc,nxarc,nfac)

             call update_svtx4contr(contr)
             allocate(ivtx_reo(nvtx)) ! dummy
             call canon_contr(contr,.false.,ivtx_reo)
             deallocate(ivtx_reo)
          
             if(ntest.ge.100)then
               write(lulog,*) 'Operator-replaced contraction'
               call prt_contr2(lulog,contr,op_info)
             endif

          end if

      end do

      return
      end
