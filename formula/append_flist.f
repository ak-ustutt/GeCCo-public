*----------------------------------------------------------------------*
      subroutine append_flist(fl1,fl2,factor,
     &                               op_info)
*----------------------------------------------------------------------*
*     input: fl1, fl2 formula lists (fl1 can be empty)
*            factor for fl2
*     on output: a concatenated formula on fl1 
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
     &     fl1, fl2
      real(8), intent(in) ::
     &     factor
      ! only for debug output:
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     list_empty
      integer ::
     &     tgt_current
      type(formula_item), pointer ::
     &     fl1_current, fl1_end,
     &     fl2_current, fl2_start

      if (ntest.ge.100) then
        write(lulog,*) '=============='
        write(lulog,*) ' append_flist '
        write(lulog,*) '=============='
      end if

      ! quick check for empty 2nd list
      ! (a) direct END
      if (fl2%command.eq.command_end_of_formula)
     /       write(lulog,*) 'quick exit 1'
      if (fl2%command.eq.command_end_of_formula) return
      ! (b) INIT END
      if (fl2%command.eq.command_set_target_init) then
         if (.not.associated(fl2%next))
     &        call quit(1,'append_flist','no END for 2nd form?')
         if (fl2%next%command.eq.command_end_of_formula)
     /       write(lulog,*) 'quick exit 2'
         if (fl2%next%command.eq.command_end_of_formula) return
      end if

      fl1_current => fl1
      tgt_current = 0

      ! go to end of fl1
      loop1: do

        if (fl1_current%command.eq.command_set_target_init) then
          tgt_current = fl1_current%target
        end if
        if (fl1_current%command.eq.command_end_of_formula) then
          exit loop1
        end if

        if (.not.associated(fl1_current%next))
     &       call quit(1,'append_flist','Sudden end of list')

        fl1_current => fl1_current%next

      end do loop1
c dbg
      write(lulog,*) 'tgt_current = ',tgt_current
c dbg

      list_empty = (.not.associated(fl1_current%prev))

      if (.not.list_empty) then
        if (.not.associated(fl1_current%prev))
     &        call quit(1,'append_flist','Inconsistent lists?')
        fl1_end => fl1_current%prev  ! omit END
        call delete_fl_node(fl1_current)

        if (associated(fl1_current)) deallocate(fl1_current)

      end if

      fl2_start => fl2%next%prev
      fl2_current => fl2%next%prev
      ! look whether fl2 starts with INIT
      if (fl2_current%command.eq.command_set_target_init) then

c dbg
        write(lulog,*) 'compare to: ',fl2_current%target
c dbg
        if (tgt_current == fl2_current%target) then
          if (list_empty)
     &        call quit(1,'append_flist','Should not land here')
          if (.not.associated(fl2_current%next))
     &        call quit(1,'append_flist','Sudden end of list(2)')
c dbg
        write(lulog,*) 'removing init. next cmd: ',
     &         fl2_current%next%command
c dbg

          fl2_start => fl2_current%next

          call delete_fl_node(fl2_current)

          ! this here fails, not clear why ... accepting minor mem leak for now
          !if (associated(fl2_current)) deallocate(fl2_current)
        end if
      end if

      ! join the two lists
      if (list_empty) then
        call delete_fl_node(fl1_current)
        fl1%command = fl2_start%command
        fl1%target = fl2_start%target
        if (fl2_start%command.ne.command_set_target_init)
     &      call quit(1,'append_flist','not expected')
        fl1%prev => null()
        fl1%next => fl2_start%next
      else
        fl1_end%next => fl2_start
        fl2_start%prev => fl1_end
c dbg
        write(lulog,*) 'at the link:'
        write(lulog,*) ' a: ',fl1_end%command,fl1_end%next%command
        write(lulog,*) ' b: ',fl2_start%prev%command,fl2_start%command
c dbg
      end if

      if (factor==1d0) return

      ! loop over second list an apply factor
      fl2_current => fl2_start
      loop2: do

        if (fl2_current%command.eq.command_add_contribution) then
          if (.not.associated(fl2_current%contr))
     &       call quit(1,'append_flist','contr not associated?')
          fl2_current%contr%fac =
     &            fl2_current%contr%fac * factor
        end if
        if (fl2_current%command.eq.command_end_of_formula) then
          exit loop2
        end if

        if (.not.associated(fl2_current%next))
     &       call quit(1,'append_flist','Sudden end of list(3)')

        fl2_current => fl2_current%next

      end do loop2

      return
      end
