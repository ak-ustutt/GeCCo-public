*----------------------------------------------------------------------*
      subroutine expand_subexpr(fl_tgt,fl_intm,mode,op_info)
*----------------------------------------------------------------------*
*     input: a definition of an intermediate on fl_intm
*            a target formula on fl_tgt
*     mode: if there are no terms for an intermediate block:
*           mode=0: quit with error message
*           mode=1: delete term
*           mode=2: ignore (keep term)
*     find all occurences of the intermediate vertex in fl_tgt and
*     expand them by the definition given in fl_intm
*     on output: an expanded formula on fl_tgt
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
      
      type(formula_item), target, intent(in) ::
     &     fl_intm
      type(formula_item), target, intent(inout) ::
     &     fl_tgt
      integer, intent(in)::
     &     mode
      ! only for debug output:
      type(operator_info), intent(in) ::
     &     op_info
      logical ::
     &     success, advance, adj_intm
      integer ::
     &     idxop_tgt, iblk_tgt, idxop_intm, iblk_intm, ivtx, nterms,
     &     njoined, iterm, ivtx_old
      type(formula_item), pointer ::
     &     fl_tgt_current, fl_tgt_current_next, fl_intm_pnt, fl_expand
      type(formula_item_list), pointer ::
     &     fpl_intm_c2blk
      type(operator), pointer ::
     &     op_intm

      integer, external ::
     &     vtx_in_contr

      if (ntest.ge.100) then
        write(luout,*) '==============================='
        write(luout,*) ' expand_subexpr messing around'
        write(luout,*) '==============================='
      end if

      if (fl_tgt%command.ne.command_set_target_init)
     &       call quit(1,'expand_subexpr',
     &       'target formula definition must start with [INIT]')

      ! first item should define new operator target
      if (fl_intm%command.ne.command_set_target_init)
     &     call quit(1,'expand_subexpr',
     &     'intermediate definition must start with [INIT]')

      idxop_intm = abs(fl_intm%target)
      adj_intm = fl_intm%target.lt.0
      op_intm => op_info%op_arr(idxop_intm)%op
      njoined = op_intm%njoined

      fl_tgt_current => fl_tgt
      ! loop over target items
      iterm = 0
      ivtx_old = 0
      tgt_loop: do

        ! new operator target ?
        if (fl_tgt_current%command.eq.command_set_target_init) then
          idxop_tgt = fl_tgt_current%target
          if (.not.associated(fl_tgt_current%next))
     &         call quit(1,'expand_subexpr',
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

        ! is intermediate vertex contained in terms?
        ! If we had the same term before: search for next appearance
        ivtx = vtx_in_contr(idxop_intm,adj_intm,ivtx_old+1,
     &                      fl_tgt_current%contr)
        ivtx_old = 0
        advance = .true.

        if (ntest.ge.100) then
          write(luout,*) 'idxop_intm = ',idxop_intm
          write(luout,*) 'ivtx       = ',ivtx
        end if

        if (ivtx.gt.0) then
          iblk_intm = fl_tgt_current%contr%vertex(ivtx)%iblk_op
          if (njoined.gt.1)
     &         iblk_intm = (iblk_intm-1)/njoined+1

          ! if yes, look for first block in intermediate ...
          fl_intm_pnt => fl_intm
          success = .false.
          do while (associated(fl_intm_pnt%next))
            fl_intm_pnt => fl_intm_pnt%next
            if (fl_intm_pnt%command.eq.command_end_of_formula)
     &           exit
            if (fl_intm_pnt%contr%iblk_res.eq.iblk_intm) then
              success = .true.
              exit
            end if
          end do
          
          if (success) then
          ! ... and collect all terms
          allocate(fpl_intm_c2blk)
          call init_formula_plist(fpl_intm_c2blk)
c dbg
c          print *,'bef. calling collect_contr2block:'
c          call print_form_list(luout,fl_intm_pnt,op_info)
c dbg
          call collect_contr2block(fpl_intm_c2blk,nterms,
     &                                       fl_intm_pnt,op_info)          

          ! generate new terms
          allocate(fl_expand)
          call init_formula(fl_expand)
          call expand_term(fl_expand,nterms,ivtx,
     &         njoined,fl_tgt_current,fpl_intm_c2blk,.false.,op_info)
          end if
          
          if (nterms.gt.0.and.success) then
c dbg
            iterm = iterm + 1
c            if (mod(iterm,10).eq.0) print *,'insertion # ',iterm
c dbgend
            if (ntest.ge.100) then
              write(luout,*) 'inserting ',nterms,' new terms'
c dbg
c              print *,'new terms:'
c              call print_form_list(luout,fl_expand,op_info)
c dbg
            end if
c dbg
c            ! already remove forbidden terms
c            call select_mrcc_lag(fl_expand,(/'H','T','CUM'/),3,
c     &                           '---',op_info)
c            if (nterms.le.1) print *,' got',nterms,'term(s)'
c dbgend

            ! sum terms in expanded formula (saves time)
            if (nterms.gt.1) call sum_terms(fl_expand,op_info)

            ! and replace current term
            call replace_fl_node(fl_tgt_current,fl_expand)
            if (associated(fl_tgt_current%contr)) then
              call dealloc_contr(fl_tgt_current%contr)
              deallocate(fl_tgt_current%contr)
            end if
            deallocate(fl_tgt_current)

            fl_tgt_current => fl_expand
c dbg
c            print *,'current target formula'
c            call print_form_list(luout,fl_tgt,op_info)
c dbg
            ! we re-visit the generated terms (for multiple expansions)
            advance = .false.
          else
            advance = .true.
            select case(mode)
            case(0)
              call quit(1,'expand_subexpr',
     &                  'no terms! block not defined?')
            case(1)
              fl_tgt_current => fl_tgt_current%prev
              fl_tgt_current_next => fl_tgt_current%next
              call delete_fl_node(fl_tgt_current_next)
              deallocate(fl_tgt_current_next)
            case default
              ! check same term again
              ivtx_old = ivtx
              advance = .false.
            end select
          end if

          if (success) then
          call dealloc_formula_plist(fpl_intm_c2blk)
          deallocate(fpl_intm_c2blk)
          end if

        end if

        if (advance) then
          ! advance to next item
          ! end of formula list?
          if (.not.associated(fl_tgt_current%next)) exit tgt_loop
          ! go to next item
          fl_tgt_current => fl_tgt_current%next
          if (fl_tgt_current%command.eq.command_end_of_formula)
     &                                             exit tgt_loop
        end if

      end do tgt_loop
        
      return
      end
