      subroutine form_op_replace(opin,opout,strict,form,op_info)
*-----------------------------------------------------------------------
*     Routine which loops over a formula, form, replacing the operator, 
*     opin, with opout. This is useful for replacing formal 
*     intermediate-type operators with their actual counterparts.
*     GWR November 2007
*     CAVEAT: opin and opout are strings (operator names)!
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

      character*(*), intent(in) ::
     &     opin, opout
      logical, intent(in) ::
     &     strict
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
     &     ieqvfac, nvtx, narc, nxarc, nfac, njoined, idx_join
      integer, allocatable ::
     &     occ_temp(:,:,:), vtx_chng_idx(:)
      integer, pointer ::
     &     ivtx_reo(:), occ_vtx(:,:,:)
      logical ::
     &     reo, change, remove
      logical, pointer ::
     &     fix_vtx(:)

      integer, external ::
     &     idx_oplist2, iblk_occ

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'form_op_replace')
        write(luout,*) 'opin:  "',trim(opin),'"'
        write(luout,*) 'opout: "',trim(opout),'"'
      end if


      ! Get the required operators.
      idxin = idx_oplist2(trim(opin),op_info)
      idxout = idx_oplist2(trim(opout),op_info)
      if (idxin.le.0.or.idxout.le.0) then
        write(luout,*) 'opin:  "',trim(opin) ,'" -> ',idxin
        write(luout,*) 'opout: "',trim(opout),'" -> ',idxout
        call quit(1,'form_op_replace','error')
      end if
        
      opin_pnt => op_info%op_arr(idxin)%op
      opout_pnt => op_info%op_arr(idxout)%op
      ! Allocate a temporary array.
      njoined = opout_pnt%njoined
      allocate(occ_temp(ngastp,2,njoined),vtx_chng_idx(njoined))

      ! make sure that the operators match:
      if (opin_pnt%njoined.ne.njoined)
     &     call quit(1,'form_op_replace',
     &     'the shape of the operators does not match: '//
     &     trim(opin)//' '//trim(opout))

      form_pnt => form

      do
        ! save pointer to next node
        form_pnt_next => form_pnt%next

        ! Navigate to the correct parts of the formula.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.100) write(luout,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.100)
     &         write(luout,*) '[INIT TARGET]',form_pnt%target
        case(command_add_contribution)
          ! The necessary contractions are here.
c          write(luout,*) '[ADD]'

          change = .false.
          remove = .false.
          ! Loop over the contraction's vertices.
          nvtx = form_pnt%contr%nvtx

          if (njoined.eq.1) then
            ! simplified version for single vertex operators
            ! allows for several instances of the operator
            do idx = 1,nvtx
              ! Check the operator index of the vertex.
              idx_form_op = form_pnt%contr%vertex(idx)%idx_op

              ! If the index of the operator vertex equals that of the
              ! intermediate operator...
              if (idx_form_op.eq.idxin) then
                idx_form_blk = form_pnt%contr%vertex(idx)%iblk_op

                occ_temp(1:ngastp,1:2,1)=
     &             opin_pnt%ihpvca_occ(1:ngastp,1:2,idx_form_blk)

                ! apply the change directly:
                ! Locate the formal block's counterpart in the actual 
                ! operator. 
                idx_blk_out =
     &               iblk_occ(occ_temp,.false.,opout_pnt)

                if (idx_blk_out.le.0.and.strict) then
                  write(luout,*) trim(opin),' block no. ', idx_form_blk
                  call wrt_occ(luout,occ_temp)
                  call quit(1,'form_op_replace',
     &                 'There is no block of '//trim(opout)//
     &                 ' that corresponds to the present block of '//
     &                 trim(opin)//'!')
                else if (idx_blk_out.le.0) then
                  ! not strict: remove that term
                  remove = .true.
                end if

                form_pnt%contr%vertex(idx)%idx_op = idxout
                form_pnt%contr%vertex(idx)%iblk_op =
     &               idx_blk_out
              end if

            end do

          else
            ! more complicated operator?
            ! will not neccessarly work if more than one instance
            ! is present (which happens rarely)
          idx_join = 0
          do idx = 1,nvtx
            ! Check the operator index of the vertex.
            idx_form_op = form_pnt%contr%vertex(idx)%idx_op

            ! If the index of the operator vertex equals that of the
            ! intermediate operator...
            if(idx_form_op.eq.idxin)then
              idx_form_blk = form_pnt%contr%vertex(idx)%iblk_op

              ! Keep a tally of the number of vertices of opin that have
              ! been found, and store their indices within the 
              ! contraction.
              idx_join = idx_join+1
              if(idx_join.gt.njoined)
     &             call quit(1,'form_op_replace','idx_join gt njoined')
              vtx_chng_idx(idx_join) = idx

              ! Copy the flagged vertex to a temporary array.
              occ_temp(1:ngastp,1:2,idx_join)=
     &             opin_pnt%ihpvca_occ(1:ngastp,1:2,idx_form_blk)

              ! Note that a change has to be made to the contraction.
              change = .true.
            endif
          enddo
          ! end of old code
          end if

          if (remove) then
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          else if (change) then
            ! Locate the formal block's counterpart in the actual 
            ! operator. 
            idx_blk_out =
     &           iblk_occ(occ_temp,.false.,opout_pnt)

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
          
            allocate(ivtx_reo(nvtx),fix_vtx(nvtx),
     &           occ_vtx(ngastp,2,nvtx))
            fix_vtx = .true.
            call occvtx4contr(1,occ_vtx,form_pnt%contr,op_info)
          
            call topo_contr(ieqvfac,reo,ivtx_reo,form_pnt%contr,
     &           occ_vtx,fix_vtx)
          
            call canon_contr(form_pnt%contr,reo,ivtx_reo)
            deallocate(ivtx_reo,fix_vtx,occ_vtx)
          
            if(ntest.ge.100.and.change)then
              write(luout,*) 'Operator-replaced contraction'
              call prt_contr2(luout,form_pnt%contr,op_info)
            endif

          end if ! .not.remove

        case default
          write(luout,*) 'command = ',form_pnt%command
          call quit(1,'form_op_replace','command undefined here')
        end select

        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
