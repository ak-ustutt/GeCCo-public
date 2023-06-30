      subroutine op_replace(idx_opin,dagin,idx_opout,dagout,
     &     strict,leave,form,op_info)
*-----------------------------------------------------------------------
*     Routine which loops over a formula, form, replacing the operator, 
*     opin, with opout. This is useful for replacing formal 
*     intermediate-type operators with their actual counterparts.
*     adapted from GWR November 2007
*     now with idx_opin, idx_opout indices of operators and supporting
*     transposition
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
      logical, intent(in) ::
     &     strict, leave, dagin, dagout
      type(formula_item), target, intent(inout) ::
     &     form
      type(operator_info), intent(in) ::
     &     op_info

      type(operator), pointer ::
     &     opin_pnt, opout_pnt
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next
      integer ::
     &     idxin, idxout, idx, idx_form_op, idx_form_blk, idx_blk_out,
     &     ieqvfac, nvtx, narc, nxarc, nfac, njoined, idx_join, ij, jdx,
     &     isuper
      integer, allocatable ::
     &     occ_temp(:,:,:), vtx_chng_idx(:)
      integer, pointer ::
     &     ivtx_reo(:), occ_vtx(:,:,:), vtx_where(:)
      logical ::
     &     reo, change, remove, dag_form_op, dagi_xor_dago, all_found
      logical, pointer ::
     &     fix_vtx(:)

      integer, external ::
     &     idx_oplist2, iblk_occ
      logical, external ::
     &     opblk_restr_cmp

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'form_op_replace')
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
     &     call quit(1,'op_replace',
     &     'the shape of the operators does not match: '//
     &     trim(opin_pnt%name)//' '//trim(opout_pnt%name))

      form_pnt => form

      do
        ! save pointer to next node
        form_pnt_next => form_pnt%next

        ! Navigate to the correct parts of the formula.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.100) write(lulog,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.100)
     &         write(lulog,*) '[INIT TARGET]',form_pnt%target
        case(command_add_contribution)
          ! The necessary contractions are here.
c          write(lulog,*) '[ADD]'

          change = .false.
          remove = .false.
          ! Loop over the contraction's vertices.
          nvtx = form_pnt%contr%nvtx

          if (.true.) then
            ! simplified version for single vertex operators
            ! allows for several instances of the operator

            ! input to the following: form_pnt%contr, idxin, idxout,
            ! dagin, op_info
            ! output: form_pnt%contr (changed), found flage
            if (.true.) then
            call contr_op_replace(all_found,form_pnt%contr,
     &                            idxin,dagin,idxout,dagout,op_info)

C            else
C            ! to be deleted code:
C            do idx = 1,nvtx
C              ! Check the operator index of the vertex.
C              idx_form_op = form_pnt%contr%vertex(idx)%idx_op
C              dag_form_op = form_pnt%contr%vertex(idx)%dagger
C
C              ! If the index of the operator vertex equals that of the
C              ! intermediate operator...
C              if (idx_form_op.eq.idxin.and.(dag_form_op.eqv.dagin)) then
C                idx_form_blk = form_pnt%contr%vertex(idx)%iblk_op
C
C                occ_temp(1:ngastp,1:2,1)=
C     &             opin_pnt%ihpvca_occ(1:ngastp,1:2,idx_form_blk)
C
C                ! if necessary, look for further vertices of super-vertex
C                if (njoined.gt.1) then
C                  isuper = form_pnt%contr%svertex(idx)
C                  if (form_pnt%contr%joined(0,isuper).ne.njoined) then
C                    write(lulog,*) 'njoined,joined(0):',
C     &                   njoined,form_pnt%contr%joined(0,isuper)
C                    call quit(1,'op_replace','inconsistency')
C                  end if
C                  vtx_where => form_pnt%contr%joined(1:njoined,isuper)
C                  if (.not.dagin) then
C                    do ij = 2, njoined
C                      occ_temp(1:ngastp,1:2,ij) =
C     &                     opin_pnt%ihpvca_occ(1:ngastp,1:2,
C     &                                              idx_form_blk-1+ij)
C                    end do
C                    idx_form_blk = (idx_form_blk-1)/njoined + 1
C                  else
C                    do ij = 1, njoined
C                      occ_temp(1:ngastp,1:2,ij) =
C     &                     opin_pnt%ihpvca_occ(1:ngastp,1:2,
C     &                                        idx_form_blk-njoined+ij)
C                    end do
C                    idx_form_blk = (idx_form_blk-njoined)/njoined + 1
C                  end if
C                end if
C
C                ! apply the change directly:
C                ! Locate the formal block's counterpart in the actual 
C                ! operator. 
C                idx_blk_out =
C     &               iblk_occ(occ_temp,dagi_xor_dago,opout_pnt,
C     &                        opin_pnt%blk_version(idx_form_blk))
C
C                ! new: also check the restrictions
C                if (idx_blk_out.gt.0) then
C                  if(.not.opblk_restr_cmp(opin_pnt,idx_form_blk,dagin,
C     &                                    opout_pnt,idx_blk_out,dagout))
C     &                 idx_blk_out = -1
C                end if
C
C                if (idx_blk_out.le.0.and.strict) then
C                  write(lulog,*) trim(opin_pnt%name),
C     &                 ' block no. ', idx_form_blk
C                  call wrt_occ(lulog,occ_temp)
C                  call quit(1,'op_replace',
C     &                 'There is no block of '//trim(opout_pnt%name)//
C     &                 ' that corresponds to the present block of '//
C     &                 trim(opin_pnt%name)//'!')
C                else if (idx_blk_out.le.0) then
C                  ! not strict: remove that term if requested, 
C                  !             else leave it as is
C                  remove = .not.leave
C                end if
C
C                if (idx_blk_out.gt.0) then
C                  if (njoined.eq.1) then
C                    form_pnt%contr%vertex(idx)%idx_op = idxout
C                    form_pnt%contr%vertex(idx)%iblk_op =
C     &                   idx_blk_out
C                    form_pnt%contr%vertex(idx)%dagger =
C     &                   dagout
C                  else if (.not.dagout) then
C                    do ij = 1, njoined
C                      jdx = vtx_where(ij)
C                      form_pnt%contr%vertex(jdx)%idx_op = idxout
C                      form_pnt%contr%vertex(jdx)%iblk_op =
C     &                     (idx_blk_out-1)*njoined+ij
C                      form_pnt%contr%vertex(jdx)%dagger =
C     &                     dagout
C                    end do
C                  else
C                    do ij = 1, njoined
C                      jdx = vtx_where(ij)
C                      form_pnt%contr%vertex(jdx)%idx_op = idxout
C                      form_pnt%contr%vertex(jdx)%iblk_op =
C     &                     (idx_blk_out)*njoined-ij+1
C                      form_pnt%contr%vertex(jdx)%dagger =
C     &                     dagout
C                    end do
C                  end if
C
C                end if
C              end if
C
C            end do
C            ! up to here to be deleted
            end if

            ! handle the case if nothing was found:
            if (.not.all_found) then
               if (strict) then  
                  write(lulog,*) trim(opin_pnt%name)
                  call quit(1,'op_replace',
     &                 'There is no block of '//trim(opout_pnt%name)//
     &                 ' that corresponds to the present block of '//
     &                 trim(opin_pnt%name)//'!')
               else if (idx_blk_out.le.0) then
                  ! not strict: remove that term if requested, 
                  !             else leave it as is
                  remove = .not.leave
               end if
            end if

          else
            call quit(1,'op_replace','old route, does it work?')
            ! more complicated operator?
            ! will not neccessarly work if more than one instance
            ! is present (which happens rarely)
          end if

          if (remove) then
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          else if (change) then
            ! Locate the formal block's counterpart in the actual 
            ! operator. 
            idx_blk_out =
     &           iblk_occ(occ_temp,.false.,opout_pnt,
     &                    opin_pnt%blk_version((idx_form_blk-1)/
     &                    njoined+1))

            ! Replace the old indices with the new.
            do idx = 1,njoined
              idx_join = vtx_chng_idx(idx)
              form_pnt%contr%vertex(idx_join)%idx_op = idxout
              form_pnt%contr%vertex(idx_join)%iblk_op =
     &             (idx_blk_out-1)*njoined+idx
            enddo
          endif

          if (.not.remove) then
            ! Ensure everything is properly set up.
            narc = form_pnt%contr%narc
            nfac = form_pnt%contr%nfac
            nxarc = form_pnt%contr%nxarc
            call resize_contr(form_pnt%contr,nvtx,narc,nxarc,nfac)

            call update_svtx4contr(form_pnt%contr)
          
cmh            allocate(ivtx_reo(nvtx),fix_vtx(nvtx),
cmh     &           occ_vtx(ngastp,2,nvtx))
cmh            fix_vtx = .true.
cmh            call occvtx4contr(1,occ_vtx,form_pnt%contr,op_info)
cmh          
cmh            call topo_contr(ieqvfac,reo,ivtx_reo,form_pnt%contr,
cmh     &           occ_vtx,fix_vtx)
cmh          
cmh            call canon_contr(form_pnt%contr,reo,ivtx_reo)
cmh            deallocate(ivtx_reo,fix_vtx,occ_vtx)
            allocate(ivtx_reo(nvtx)) ! dummy
            call canon_contr(form_pnt%contr,.false.,ivtx_reo)
            deallocate(ivtx_reo)
          
            if(ntest.ge.100.and.change)then
              write(lulog,*) 'Operator-replaced contraction'
              call prt_contr2(lulog,form_pnt%contr,op_info)
            endif

          end if ! .not.remove

        case default
          ! just ignore unknown commands
c          write(lulog,*) 'command = ',form_pnt%command
c          call quit(1,'form_op_replace','command undefined here')
        end select

        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
