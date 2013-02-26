*----------------------------------------------------------------------*
      subroutine del_zero_terms(fl_tgt,mode,op_info,thrsh)
*----------------------------------------------------------------------*
*     look for terms with factor below threshold and delete them
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 000

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
      real(8) ::
     &     thrsh
      character(len=*), intent(in) ::
     &     mode

      integer ::
     &     idxop_tgt, iterm
      type(formula_item), pointer ::
     &     fl_tgt_pnt, fl_tgt_pnt_next, fl_tgt_current

      ! if first entry is an [END]: do nothing
      if (fl_tgt%command.eq.command_end_of_formula) return

      if (ntest.ge.100) then
        write(luout,*) '========================='
        write(luout,*) ' info from del_zero_terms'
        write(luout,*) '========================='
        write(luout,*) ' mode = ',trim(mode)
      end if

      if (fl_tgt%command.ne.command_set_target_init) then
        if (fl_tgt%command.ne.command_add_contribution)
     &       call quit(1,'del_zero_terms',
     &         'must start with [INIT] or [ADD]')
        idxop_tgt = fl_tgt%target
      end if

      ! first sum identical terms?
      if (trim(mode).eq.'sum'.or.trim(mode).eq.'SUM')
     &   call sum_terms(fl_tgt,op_info)

      iterm = 0
      fl_tgt_current => fl_tgt
      ! loop over target items
      tgt_loop: do

        ! new operator target ?
        if (fl_tgt_current%command.eq.command_set_target_init) then
          idxop_tgt = fl_tgt_current%target
          if (.not.associated(fl_tgt_current%next))
     &         call quit(1,'del_zero_terms',
     &         'unexpected end of list (target)')
          if (ntest.ge.100) then
            write(luout,'(70("="))')
            write(luout,*) 'New operator target: ',idxop_tgt
            write(luout,'(70("="))')
          end if
          fl_tgt_current => fl_tgt_current%next
        end if

        iterm = iterm+1
        if (ntest.ge.100) then
          write(luout,*) 'current term: # ',iterm
          call prt_contr2(luout,fl_tgt_current%contr,op_info)
        end if

        fl_tgt_pnt => fl_tgt_current%next

        ! delete term with factor below threshold
        if (abs(fl_tgt_current%contr%fac).lt.thrsh) then
          call delete_fl_node(fl_tgt_current)
          deallocate(fl_tgt_current)
        end if

        ! advance to next item
        fl_tgt_current => fl_tgt_pnt
        if (.not.associated(fl_tgt_current))
     &         call quit(1,'del_zero_terms',
     &         'unexpected end of list (target)')
        if (fl_tgt_current%command.eq.command_end_of_formula)
     &                                             exit tgt_loop

      end do tgt_loop
c dbg
c      print *,' del_zero_terms yields',iterm,' terms.'
c dbgend
        
      return
      end
