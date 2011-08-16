*----------------------------------------------------------------------*
      subroutine factor_out_subexpr2(fl_tgt,fl_intm,nrpl,op_info)
*----------------------------------------------------------------------*
*     input: a definition of an intermediate on fl_intm
*            a target formula on fl_tgt
*     find all occurences of the intermediate in fl_tgt and
*     modify fl_tgt accordingly
*     on output: a reduced formula on fl_tgt, 
*                the number of replacements on nrpl (info)
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
      integer, intent(out) ::
     &     nrpl
      ! only for debug output:
      type(operator_info), intent(in) ::
     &     op_info

      integer, parameter ::
     &     nmod_max = 10
      integer ::
     &     nterms, nposs, idxop_tgt, iblk_tgt, iterm, iblk_int,
     &     nmod, idx, imod(nmod_max)
      real(8) ::
     &     xmod(nmod_max)
      type(contraction) ::
     &     contr_rpl
      type(formula_item), pointer ::
     &     fl_tgt_current, fl_intm_pnt
      type(formula_item_list), pointer ::
     &     fpl_intm_c2blk,
     &     fpl_intm_start,
     &     fpl_intm_current,
     &     fpl_intm_in_tgt,
     &     fpl_intm_in_tgt_pnt
      logical ::
     &     success
      integer, pointer ::
     &     iposs_blk(:)
      integer, external ::
     &     idxlist

      if (ntest.ge.100) then
        write(luout,*) '==================================='
        write(luout,*) ' factor_out_subexpr messing around'
        write(luout,*) '==================================='
      end if

      nrpl = 0   

      call init_contr(contr_rpl)

      if (fl_tgt%command.ne.command_set_target_init)
     &       call quit(1,'factor_out_subexpr',
     &       'target formula definition must start with [INIT]')

      ! first item should define new operator target
      if (fl_intm%command.ne.command_set_target_init)
     &     call quit(1,'factor_out_subexpr',
     &     'intermediate definition must start with [INIT]')

      fl_tgt_current => fl_tgt
      ! loop over target items
      tgt_loop: do
c dbg
c        success = .false.
c dbgend

        ! new operator target ?
        if (fl_tgt_current%command.eq.command_set_target_init) then
          idxop_tgt = fl_tgt_current%target
          if (.not.associated(fl_tgt_current%next))
     &         call quit(1,'factor_out_subexpr',
     &         'unexpected end of list (target)')
          if (ntest.ge.100) then
            write(luout,'(70("="))')
            write(luout,*) 'New operator target: ',idxop_tgt
            write(luout,'(70("="))')
          end if
          ! for a save exit (although this normally should not happen):
          if (.not.associated(fl_tgt_current%next)) exit tgt_loop

          fl_tgt_current => fl_tgt_current%next
          cycle tgt_loop
        else if (fl_tgt_current%command.eq.command_end_of_formula) then
          ! may happen for empty formulae
          exit tgt_loop
        end if

        if (.not.associated(fl_tgt_current%contr))
     &       call quit(1,'factor_out_subexpr','I''m confused ...')

        iblk_tgt = fl_tgt_current%contr%iblk_res
        if (ntest.ge.100) then
          write(luout,*) 'current term:'
          call prt_contr2(luout,fl_tgt_current%contr,op_info)
        end if

        ! ------------------------------------------------------------
        ! make a list of intermediates which could be contained
        ! in the current term, ordered according to largest
        ! "overlap" with current term (i.e. number of common factors)
        ! and sequence of definition in intermediate formula
        ! (first has highest priority)
        ! ------------------------------------------------------------
        allocate(fpl_intm_start)
        call find_possible_subexpr(nposs,fpl_intm_start,
     &       fl_tgt_current,fl_intm,op_info)

        if (ntest.ge.100) then
          write(luout,*) '# of possible replacements: ',nposs
        end if

        ! anything found?
        if (nposs.gt.0) then

          ! -----------------------------------------------------
          ! loop over possible replacements
          ! and find out whether the remaining terms are 
          ! present, as well
          ! -----------------------------------------------------
          fpl_intm_current => fpl_intm_start
          allocate(iposs_blk(op_info%op_arr(fpl_intm_current%item%
     &                        contr%idx_res)%op%n_occ_cls))
          iposs_blk = 0
          poss_loop: do
            fl_intm_pnt => fpl_intm_current%item
            iblk_int = fl_intm_pnt%contr%iblk_res
            iposs_blk(iblk_int) = iposs_blk(iblk_int) + 1

            if (ntest.ge.100) then
              write(luout,'(x,2(a,i4))') 'poss. # ',iposs_blk(iblk_int),
     &           ' of block ',iblk_int,' (starts with:)'
              call prt_contr2(luout,fl_intm_pnt%contr,op_info)
            end if

            ! collect all contributions with same result block index
            ! in contraction list cl_intm_c2blk
            allocate(fpl_intm_c2blk)
            call init_formula_plist(fpl_intm_c2blk)
            call collect_contr2block(fpl_intm_c2blk,nterms,
     &                                       fl_intm_pnt,op_info)
            allocate(fpl_intm_in_tgt)
            call init_formula_plist(fpl_intm_in_tgt)

            ! find all blocks which originate from a contraction with
            ! current intermediate block
            ! collect blocks on array fa_intm_in_tgt
            ! the new contraction with the intermediate is on
            ! contr_rpl
            call find_contr_w_intm2(success,fpl_intm_in_tgt,contr_rpl,
     &         fl_tgt_current,fpl_intm_c2blk,iposs_blk(iblk_int),
     &         nmod_max,nmod,imod,xmod,
     &         op_info)

            if (ntest.ge.100.and..not.success) then
              write(luout,*) 'NO SUCCESS'
            end if

            if (success) then
              if (ntest.ge.100) then
                write(luout,*) 'SUCCESS'
                write(luout,*) 'now modifying formula list'
              end if
              nrpl = nrpl + 1  ! increment counter
              ! replace first node with new intermediate and delete
              ! all other nodes
              call copy_contr(contr_rpl,fpl_intm_in_tgt%item%contr)
              fpl_intm_in_tgt_pnt => fpl_intm_in_tgt
              iterm = 1
              do while(associated(fpl_intm_in_tgt_pnt%next))
                fpl_intm_in_tgt_pnt => fpl_intm_in_tgt_pnt%next
                iterm = iterm + 1
                idx = idxlist(iterm,imod,nmod,1)
                if (idx.gt.0) then ! just change factor
                  fpl_intm_in_tgt_pnt%item%contr%fac = xmod(idx)
                  cycle
                end if
                call delete_fl_node(fpl_intm_in_tgt_pnt%item)
                if (associated(fpl_intm_in_tgt_pnt%item%contr)) then
                  call dealloc_contr(fpl_intm_in_tgt_pnt%item%contr)
                  deallocate(fpl_intm_in_tgt_pnt%item%contr)
                end if
                deallocate(fpl_intm_in_tgt_pnt%item)
              end do
              
            end if

            call dealloc_formula_plist(fpl_intm_in_tgt)
            call dealloc_formula_plist(fpl_intm_c2blk)
            deallocate(fpl_intm_in_tgt)
            deallocate(fpl_intm_c2blk)
            if (success) exit poss_loop
            if (.not.associated(fpl_intm_current%next))
     &           exit poss_loop
            fpl_intm_current => fpl_intm_current%next

          end do poss_loop
          deallocate(iposs_blk)

        end if ! any possible replacement found for current term

        call dealloc_formula_plist(fpl_intm_start)
        deallocate(fpl_intm_start)

c dbg
c        ! if success: try again with same term!
c        if (.not.success) then
c dbgend
        ! advance to next item
        ! end of formula list?
        if (.not.associated(fl_tgt_current%next)) exit tgt_loop
        ! go to next item
        fl_tgt_current => fl_tgt_current%next
        if (fl_tgt_current%command.eq.command_end_of_formula)
     &                                             exit tgt_loop
c dbg
c        end if
c dbgend

      end do tgt_loop

      call dealloc_contr(contr_rpl)

      return
      end
