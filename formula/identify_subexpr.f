*----------------------------------------------------------------------*
      subroutine identify_subexpr2(fl_intm,fl_tgt,op_info)
*----------------------------------------------------------------------*
*     input: a definition of an intermediate on fl_intm
*            a target formula on fl_tgt
*     find all occurences of the intermediat in fl_tgt and
*     modify fl_tgt accordingly
*     on output: an reduced formula on fl_tgt
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

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
      ! only for debug output:
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     nterms, nposs, idxop_tgt, iblk_tgt, iterm, iposs
      type(contraction) ::
     &     contr_rpl
      type(formula_item), pointer ::
     &     fl_tgt_current, fl_intm_pnt
      type(formula_item_list), pointer ::
     &     fpl_intm_c2blk,
     &     fpl_intm_start,
     &     fpl_intm_current
      type(formula_item_array), pointer ::
     &     fpa_intm_in_tgt(:)
      logical ::
     &     success

      if (ntest.ge.100) then
        write(luout,*) '================================='
        write(luout,*) ' identify_subexpr messing around'
        write(luout,*) '================================='
      end if

      call init_contr(contr_rpl)

      if (fl_tgt%command.ne.command_set_target_init)
     &       call quit(1,'identify_subexpr',
     &       'target formula definition must start with [INIT]')

      ! first item should define new operator target
      if (fl_intm%command.ne.command_set_target_init)
     &     call quit(1,'identify_subexpr',
     &     'intermediate definition must start with [INIT]')

      fl_tgt_current => fl_tgt
      ! loop over target items
      tgt_loop: do

        ! new operator target ?
        if (fl_tgt_current%command.eq.command_set_target_init) then
          idxop_tgt = fl_tgt_current%target
          if (.not.associated(fl_tgt_current%next))
     &         call quit(1,'identify_subexpr',
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
          call prt_contr2(luout,fl_tgt_current%contr,op_info%op_arr)
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
          iposs = 0
          poss_loop: do
            iposs = iposs+1
            fl_intm_pnt => fpl_intm_current%item

            if (ntest.ge.100) then
              write(luout,*) 'poss. # ',iposs,' (starts with:)'
              call prt_contr2(luout,fl_intm_pnt%contr,op_info%op_arr)
            end if

            ! collect all contributions with same result block index
            ! in contraction list cl_intm_c2blk
            allocate(fpl_intm_c2blk)
            call collect_contr2block(fpl_intm_c2blk,nterms,
     &                                       fl_intm_pnt,op_info)
            allocate(fpa_intm_in_tgt(nterms))

            ! find all blocks which originate from a contraction with
            ! current intermediate block
            ! collect blocks on array fa_intm_in_tgt
            ! the new contraction with the intermediate is on
            ! contr_rpl
            call find_contr_w_intm(success,fpa_intm_in_tgt,contr_rpl,
     &         fl_tgt_current,fpl_intm_c2blk,nterms,
     &         op_info)

            if (ntest.ge.100.and..not.success) then
              write(luout,*) 'NO SUCCESS'
            end if

            if (success) then
              if (ntest.ge.100) then
                write(luout,*) 'SUCCESS'
                write(luout,*) 'now modifying formula list'
              end if
              ! replace first node with new intermediate and delete
              ! all other nodes
              call copy_contr(contr_rpl,fpa_intm_in_tgt(1)%item%contr)
              do iterm = 2, nterms
                call delete_fl_node(fpa_intm_in_tgt(iterm)%item)
                deallocate(fpa_intm_in_tgt(iterm)%item)
              end do
              
            end if

            deallocate(fpa_intm_in_tgt)
            call dealloc_formula_plist(fpl_intm_c2blk)
            deallocate(fpl_intm_c2blk)
            if (success) exit poss_loop
            if (.not.associated(fpl_intm_current%next))
     &           exit poss_loop
            fpl_intm_current => fpl_intm_current%next

          end do poss_loop

        end if ! any possible replacement found for current term

        call dealloc_formula_plist(fpl_intm_start)
        deallocate(fpl_intm_start)

        ! advance to next item
        ! end of formula list?
        if (.not.associated(fl_tgt_current%next)) exit tgt_loop
        ! go to next item
        fl_tgt_current => fl_tgt_current%next
        if (fl_tgt_current%command.eq.command_end_of_formula)
     &                                             exit tgt_loop

      end do tgt_loop
        
      call dealloc_contr(contr_rpl)
c dbg
c        call quit(1,'test','exit 2')
c dbg

      return
      end
