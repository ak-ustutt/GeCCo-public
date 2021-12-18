      subroutine select_formal_blk(flist,mode,op_info)
*----------------------------------------------------------------------*
*     either extract or delete terms contributing to formal blocks
*
*     matthias, oct 2011
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      character(len=*), intent(in) ::
     &     mode

      logical ::
     &     delete, formal
      integer ::
     &     idx_op, iblk

      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'select_formal_blk')
        write(lulog,*) 'mode = ',trim(mode)
      endif

      select case(mode)
      case('delete','DELETE')
        delete = .true.
      case('extract','EXTRACT','keep','KEEP')
        delete = .false.
      case default
        call quit(1,'select_formal_blk','unknown mode: '//trim(mode))
      end select

      form_pnt => flist
      do 
        form_pnt_next => form_pnt%next
        ! Locate actual formula items.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.1000) write(lulog,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(lulog,*) '[INIT_TARGET]'
        case(command_add_contribution)

          idx_op = form_pnt%contr%idx_res
          iblk = form_pnt%contr%iblk_res

          ! is the result block a formal block?
          formal = op_info%op_arr(idx_op)%op%formal_blk(iblk)

          ! either delete formal blocks or non-formal blocks
          if (delete.eqv.formal) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*) 'Deleted formula item:'
              call prt_contr2(lulog,form_pnt%contr,op_info)
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'select_formal_blk','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
      
      
