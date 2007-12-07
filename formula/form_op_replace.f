      subroutine form_op_replace(opin,opout,form,op_info)
*-----------------------------------------------------------------------
*     Routine which loops over a formula, form, replacing the operator, 
*     opin, with opout. This is useful for replacing formal 
*     intermediate-type operators with their actual counterparts.
*     GWR November 2007
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

      character*256, intent(in) ::
     &     opin, opout
      type(formula_item), target, intent(inout) ::
     &     form
      type(operator_info), intent(in) ::
     &     op_info

      type(operator), pointer ::
     &     opin_pnt, opout_pnt
      type(formula_item), pointer ::
     &     form_pnt
      integer ::
     &     idxin, idxout, idx, idx_form_op, idx_form_blk, idx_blk_out,
     &     ieqvfac, nvtx, narc, nfac, njoined, idx_join
      integer, allocatable ::
     &     occ_temp(:,:,:), vtx_chng_idx(:)
      integer, pointer ::
     &     ivtx_reo(:), occ_vtx(:,:,:)
      logical ::
     &     reo, change
      logical, pointer ::
     &     fix_vtx(:)

      integer, external ::
     &     idx_oplist2, iblk_occ


      ! Get the required operators.
      idxin = idx_oplist2(trim(opin),op_info)
      opin_pnt => op_info%op_arr(idxin)%op
      idxout = idx_oplist2(trim(opout),op_info)
      opout_pnt => op_info%op_arr(idxout)%op
      ! Allocate a temporary array.
      njoined = opout_pnt%njoined
      allocate(occ_temp(ngastp,2,njoined),vtx_chng_idx(njoined))

      form_pnt => form

      do
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
          ! Loop over the contraction's vertices.
          nvtx = form_pnt%contr%nvtx
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

          if(change)then
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
          
          ! Ensure everything is properly set up.
          narc = form_pnt%contr%narc
          nfac = form_pnt%contr%nfac
          call resize_contr(form_pnt%contr,nvtx,narc,nfac)

          call update_svtx4contr(form_pnt%contr)
          
          allocate(ivtx_reo(nvtx),fix_vtx(nvtx),
     &         occ_vtx(ngastp,2,nvtx))
          fix_vtx = .true.
          call occvtx4contr(1,occ_vtx,form_pnt%contr,op_info)
          
          call topo_contr(ieqvfac,reo,ivtx_reo,form_pnt%contr,
     &         occ_vtx,fix_vtx)
          
          call canon_contr(form_pnt%contr,reo,ivtx_reo)
          deallocate(ivtx_reo,fix_vtx,occ_vtx)
          
          if(ntest.ge.100.and.change)then
            write(luout,*) 'Operator-replaced contraction'
            call prt_contr2(luout,form_pnt%contr,op_info)
          endif

        case default
          write(luout,*) 'command = ',form_pnt%command
          call quit(1,'form_op_replace','command undefined here')
        end select

        if(.not.associated(form_pnt%next))exit
        form_pnt => form_pnt%next

      enddo

      return
      end
