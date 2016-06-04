      subroutine select_terms(flist,
     &     idxop_res,
     &     ninclude,idxop_incl,iblk_incl,
     &     ninclude_or,idxop_incl_or,iblk_incl_or,
     &     nexclude,idxop_excl,iblk_excl,
     &     op_info)
*----------------------------------------------------------------------*
*     preliminary start-up for truncated CC expansions
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

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     ninclude, ninclude_or, nexclude,
     &     idxop_res,
     &     idxop_incl(ninclude), iblk_incl(ninclude),
     &     idxop_incl_or(ninclude), iblk_incl_or(ninclude),
     &     idxop_excl(nexclude), iblk_excl(nexclude)

      logical ::
     &     delete
      integer ::
     &     nvtx, ivtx,
     &     idx_op, iblk_op, idx
      character*64 ::
     &     op_name
      logical ::
     &     found

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'select_terms')
        write(lulog,*) 'ninclude: ',ninclude
        write(lulog,*) 'op:  ',idxop_incl(1:ninclude)
        write(lulog,*) 'blk: ',iblk_incl(1:ninclude)
        write(lulog,*) 'ninclude_or: ',ninclude_or
        write(lulog,*) 'op:  ',idxop_incl_or(1:ninclude_or)
        write(lulog,*) 'blk: ',iblk_incl_or(1:ninclude_or)
        write(lulog,*) 'nexclude: ',nexclude
        write(lulog,*) 'op:  ',idxop_excl(1:nexclude)
        write(lulog,*) 'blk: ',iblk_excl(1:nexclude)
      endif

      form_pnt => flist
      do 
        form_pnt_next => form_pnt%next
        ! Locate actual formula items.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.1000) write(lulog,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(lulog,*) '[INIT_TARGET]'
          form_pnt%target = idxop_res
        case(command_add_contribution)

          if (ntest.ge.1000) then
            write(lulog,*) 'current item:'
            call prt_contr2(lulog,form_pnt%contr,op_info)
          end if

          nvtx = form_pnt%contr%nvtx
          vertex => form_pnt%contr%vertex

          delete = .false.
          do ivtx = 1, nvtx
            if (delete) exit
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op

            do idx = 1, nexclude
              delete = delete.or.(idx_op.eq.idxop_excl(idx).and.
     &                            ( 0     .eq.iblk_excl(idx).or.
     &                             iblk_op.eq.iblk_excl(idx)))
            end do
          end do

          if (ntest.ge.1000)
     &         write(lulog,*) 'after exclude: delete=',delete

          do idx = 1, ninclude
            if (delete) exit
            idx_op = idxop_incl(idx)
            iblk_op = iblk_incl(idx)

            found = .false.
            do ivtx = 1, nvtx
              found = found.or.(idx_op.eq.vertex(ivtx)%idx_op.and.
     &                          (iblk_op.eq.0.or.
     &                           iblk_op.eq.vertex(ivtx)%iblk_op))
            end do
            delete = delete.or..not.found
          end do

          if (ntest.ge.1000)
     &         write(lulog,*) 'after include: delete=',delete

          found = ninclude_or.eq.0
          do ivtx = 1, nvtx
            if (delete) exit
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op

            do idx = 1, ninclude_or
              found = found.or.(idx_op.eq.idxop_incl_or(idx).and.
     &                          ( 0     .eq.iblk_incl_or(idx).or.
     &                           iblk_op.eq.iblk_incl_or(idx)))
            end do
            if (found) exit
          end do
          delete = delete.or..not.found

          if (ntest.ge.1000)
     &         write(lulog,*) 'after include_or: delete=',delete

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*)'Deleted formula item'
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          else
            form_pnt%target = idxop_res
            form_pnt%contr%idx_res =
     &           idxop_res
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'select_terms','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
      
      
