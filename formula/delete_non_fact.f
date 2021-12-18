      subroutine delete_non_fact(ops,nops,form,op_info)
*-----------------------------------------------------------------------
*     Routine which scans a formula (form) searching for occurrences of 
*     the operators contained in ops(nops). If an operator is present,
*     the formula item is deleted: if none of the operators
*     is present, the item is retained.
*     GWR November 2007.
*-----------------------------------------------------------------------

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      
      integer, intent(in) ::
     &     nops
      character(len=64), intent(in) ::
     &     ops(nops)
      type(formula_item), target, intent(inout) ::
     &     form
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     form_pnt
      integer ::
     &     idx, jdx, nvtx
      logical ::
     &     redundant
      integer ::
     &     idx_ops(nops)
      integer, external ::
     &     idx_oplist2

      if(ntest.ge.100)then
        write(lulog,*)'============================='
        write(lulog,*)' Deleting non-factored terms '
        write(lulog,*)'============================='
      endif

      ! Associate the input operator names to their integer indices.
      do idx=1,nops
        idx_ops(idx) = idx_oplist2(trim(ops(idx)),op_info)
        if(idx_ops(idx).lt.0)
     &       call quit(1,'delete_non_fact','delete op not defined')
      enddo

      form_pnt => form
      do
        ! Locate actual formula items.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.1000) write(lulog,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(lulog,*) '[INIT_TARGET]'
        case(command_add_contribution)
          ! The necessary contractions are in these parts of the formula.
          if(ntest.ge.1000)then
            write(lulog,*) '[ADD]'
            call prt_contr2(lulog,form_pnt%contr,op_info)
          endif

          nvtx = form_pnt%contr%nvtx
          ! See if the contraction contains the prohibited operators.
          redundant = .false.
          del_loop: do idx=1,nops
            do jdx=1,nvtx
              if(form_pnt%contr%vertex(jdx)%idx_op.eq.
     &             idx_ops(idx))then
                redundant = .true.
                exit del_loop
              endif
            enddo
          enddo del_loop

          if(redundant)then
            call delete_fl_node(form_pnt)
            if(associated(form_pnt%contr))then
              call dealloc_contr(form_pnt%contr)
              deallocate(form_pnt%contr)
            endif
          endif

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'delete_non_fact','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt%next))exit
        form_pnt => form_pnt%next

      enddo

      if(ntest.ge.100)then
        call write_title(lulog,wst_title,'Non-redundant R12 Terms')
        call print_form_list(lulog,form,op_info)
      endif

      return
      end
