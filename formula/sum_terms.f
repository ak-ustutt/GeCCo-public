*----------------------------------------------------------------------*
      subroutine sum_terms(fl_tgt,op_info)
*----------------------------------------------------------------------*
*     look for equal terms and sum them
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula_item_array.h'
      include 'def_formula_item_list.h'
      
      type(formula_item), target, intent(inout) ::
     &     fl_tgt
      ! only for debug output:
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idxop_tgt, iblk_tgt
      type(formula_item), pointer ::
     &     fl_tgt_pnt, fl_tgt_pnt_next, fl_tgt_current

      logical, external ::
     &     cmp_contr

      if (ntest.ge.100) then
        write(luout,*) '====================='
        write(luout,*) ' info from sum_terms'
        write(luout,*) '====================='
      end if

      if (fl_tgt%command.ne.command_set_target_init)
     &       call quit(1,'sum_terms',
     &       'target formula definition must start with [INIT]')

      fl_tgt_current => fl_tgt
      ! loop over target items
      tgt_loop: do

        ! new operator target ?
        if (fl_tgt_current%command.eq.command_set_target_init) then
          idxop_tgt = fl_tgt_current%target
          if (.not.associated(fl_tgt_current%next))
     &         call quit(1,'sum_terms',
     &         'unexpected end of list (target)')
          if (ntest.ge.100) then
            write(luout,'(70("="))')
            write(luout,*) 'New operator target: ',idxop_tgt
            write(luout,'(70("="))')
          end if
          fl_tgt_current => fl_tgt_current%next
        end if

        iblk_tgt = fl_tgt_current%contr%iblk_res
        if (ntest.ge.100) then
          write(luout,*) 'current term:'
          call prt_contr2(luout,fl_tgt_current%contr,op_info)
        end if

        if (.not.associated(fl_tgt_current%next))
     &         call quit(1,'sum_terms',
     &         'unexpected end of list (target)')
        fl_tgt_pnt => fl_tgt_current%next
        search_loop: do
          ! we search within a single domain only
          if (fl_tgt_pnt%command.ne.command_add_contribution)
     &       exit search_loop

          ! as the node might get deleted, better save the next pointer
          if (.not.associated(fl_tgt_pnt%next))
     &         call quit(1,'sum_terms',
     &         'unexpected end of list (target, inner loop)')
          fl_tgt_pnt_next => fl_tgt_pnt%next
          
          if (fl_tgt_pnt%contr%iblk_res.eq.iblk_tgt) then
            if (fl_tgt_pnt%contr%idx_res.ne.idxop_tgt)
     &           call quit(1,'sum_terms',
     &           'suspicious change of operator target')

            if (cmp_contr(fl_tgt_pnt%contr,
     &                    fl_tgt_current%contr,.true.)) then
              if (ntest.ge.100) then
                write(luout,*) 'found equal term:'
                call prt_contr2(luout,fl_tgt_current%contr,
     &               op_info)
                write(luout,*) 'now summing and deleting'
              end if
              fl_tgt_current%contr%fac =
     &             fl_tgt_current%contr%fac + fl_tgt_pnt%contr%fac
              call delete_fl_node(fl_tgt_pnt)
c dbg
c              print *,'mark top'
c              call print_form_list(luout,fl_tgt,op_info)
c              print *,'mark bot'
c dbg
              deallocate(fl_tgt_pnt)
            end if
          end if

          fl_tgt_pnt => fl_tgt_pnt_next
        end do search_loop

        ! advance to next item
        ! end of formula list?
        if (.not.associated(fl_tgt_current%next)) exit tgt_loop
        ! go to next item
        fl_tgt_current => fl_tgt_current%next
        if (fl_tgt_current%command.eq.command_end_of_formula)
     &                                             exit tgt_loop

      end do tgt_loop
        
      return
      end
