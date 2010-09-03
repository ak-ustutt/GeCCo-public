      subroutine r12_opt_truncation(flist,idxop1,idxop2,op_info)
*----------------------------------------------------------------------*
*     delete terms which contain operator blocks of wrong version
*     op1: only blocks with version = 1 allowed
*     op2: only blocks with version = 2 allowed
*
*     matthias, april 2009
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'def_del_list.h'
      include 'par_opnames_gen.h'

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     idxop1, idxop2

      logical ::
     &     delete
      integer ::
     &     nvtx, ivtx, njoined,
     &     idx_op, iblk_op, version

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      call quit(1,'r12_opt_truncation','call to obsolete routine')

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'r12_opt_truncation')
      endif

      form_pnt => flist
      do 
        form_pnt_next => form_pnt%next
        ! Locate actual formula items.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.1000) write(luout,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(luout,*) '[INIT_TARGET]'
        case(command_add_contribution)

          nvtx = form_pnt%contr%nvtx
          vertex => form_pnt%contr%vertex

          delete = .false.
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op
            njoined = op_info%op_arr(idx_op)%op%njoined
            version = op_info%op_arr(idx_op)%op%
     &                        blk_version((iblk_op-1)/njoined+1)
            delete = delete.or.idx_op.eq.idxop1.and.version.ne.1.or.
     &                         idx_op.eq.idxop2.and.version.ne.2
          end do

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(luout,*) 'Deleted formula item:'
              call prt_contr2(luout,form_pnt%contr,op_info)
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(luout,*)'command = ',form_pnt%command
          call quit(1,'r12_opt_truncation','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo


      return
      end
      
      
