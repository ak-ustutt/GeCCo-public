*----------------------------------------------------------------------*
      subroutine find_possible_subexpr(nposs,fpl_intm_start,
     &     fl_tgt,fl_intm,op_info)
*----------------------------------------------------------------------*
*     given a term of a target expression (on fl_tgt)
*     and a list of definitions for intermediates (on fl_intm),
*     find all possible intermediates that may be used for factoring
*     the target formula, i.e. all intermediates that potentially
*     are contained in the present term
*     the possibilities are ranked wrt the length, e.g.
*
*     target ABCD
*     rank1   BCD
*     rank2  ABC   (because it was defined second in fl_intm)
*     rank3    CD
*     rank4  A
*
*     Andreas, June 2007
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'

      integer, parameter ::
     &     ntest = 00,
     &     maxposs = 128

      integer, intent(out) ::
     &     nposs
      type(formula_item_list), target ::
     &     fpl_intm_start
      type(formula_item), target ::
     &     fl_tgt, fl_intm
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     first
      integer ::
     &     iblk_intm, idx, jdx, len_i
      integer ::
     &     len_list(maxposs)
      type(formula_item), pointer ::
     &     fl_intm_pnt, fl_intm_stblk
      type(formula_item_list), pointer ::
     &     fpl_intm_start_pnt,
     &     fpl_intm_new

      logical, external ::
     &     contr_in_contr

      if (ntest.ge.100) then
        write(luout,*) '================================='
        write(luout,*) 'output from find_possible_subexpr'
        write(luout,*) '================================='
      end if

      nposs = 0
      nullify(fpl_intm_start%prev)
      nullify(fpl_intm_start%next)
      first = .true.
      fl_intm_pnt => fl_intm
      ! loop over all blocks of the intermediate
      do

        if (fl_intm_pnt%command.eq.command_set_target_init) then
          if (.not.associated(fl_intm_pnt%next))
     &         call quit(1,'find_possible_subexpr',
     &         'unexpected end of list (intermediate)')
          fl_intm_pnt => fl_intm_pnt%next
          if (fl_intm_pnt%command.ne.command_add_contribution)
     &         call quit(1,'find_possible_subexpr',
     &         'expected [ADD] (intermediate)')
          iblk_intm = fl_intm_pnt%contr%iblk_res
          ! new block -> remember the start item
          fl_intm_stblk => fl_intm_pnt
        else if (fl_intm_pnt%command.eq.command_add_contribution) then
          if (iblk_intm.gt.fl_intm_pnt%contr%iblk_res)
     &         call quit(1,'find_possible_subexpr',
     &         'expected intermediate blocks in increasing sequence')
          ! new block? remember the start item
          if (iblk_intm.ne.fl_intm_pnt%contr%iblk_res)
     &         fl_intm_stblk => fl_intm_pnt
          iblk_intm = fl_intm_pnt%contr%iblk_res          
        else if (fl_intm_pnt%command.eq.command_end_of_formula) then
          exit
        else
          call quit(1,'find_possible_subexpr',
     &         'unexpected command')
        end if

        ! is intermediate contained in current term?
        if (contr_in_contr(fl_intm_pnt%contr,fl_tgt%contr)) then
          nposs = nposs+1
          if (nposs.gt.maxposs)
     &         call quit(1,'find_possible_subexpr',
     &         'More than <maxposs> possibilities: Are you sure?')
          len_i = fl_intm_pnt%contr%nvtx
          if (first) then
            first = .false.
            fpl_intm_start%item => fl_intm_stblk
            nullify(fpl_intm_start%prev)
            nullify(fpl_intm_start%next)
            len_list(1) = len_i
          else
            ! compare with possibilities on list
            idx = 0
            fpl_intm_start_pnt => fpl_intm_start
            do
              idx = idx+1
              if (len_i.gt.len_list(idx)) then
                do jdx = nposs, idx+1, -1
                  len_list(jdx)=len_list(jdx-1)
                end do
                len_list(idx) = len_i
                exit
              else if (idx.eq.nposs) then
                len_list(idx) = len_i
                exit
              end if
              ! point to node, after which we might append
              if (idx.gt.1) then
                if (.not.associated(fpl_intm_start_pnt%next)) 
     &             call quit(1,'find_possible_subexpr','inconsistency')
                fpl_intm_start_pnt => fpl_intm_start_pnt%next 
              end if
            end do
            ! allocate a new instance
            allocate(fpl_intm_new)
            if (idx.eq.1) then
              ! first node special:
              ! copy entries of first node to second node
              fpl_intm_new%item => fpl_intm_start%item
              fpl_intm_new%next => fpl_intm_start%next
              ! and insert it as second node
              fpl_intm_new%prev => fpl_intm_start
              ! put new entry to first node
              fpl_intm_start%item => fl_intm_stblk
              fpl_intm_start%next => fpl_intm_new
            else
              ! new entry on new node
              fpl_intm_new%item => fl_intm_stblk
              ! and insert into list
              fpl_intm_new%prev => fpl_intm_start_pnt
              fpl_intm_new%next => fpl_intm_start_pnt%next
              if (associated(fpl_intm_start_pnt%next)) then
                fpl_intm_new%next%prev => fpl_intm_new
              end if
              fpl_intm_start_pnt%next => fpl_intm_new
            end if
          end if  
        end if


        ! advance intermediate pointer
        if (.not.associated(fl_intm_pnt%next)) exit
        fl_intm_pnt => fl_intm_pnt%next

      end do

      if (nposs.gt.0.and.ntest.ge.100) then
        fpl_intm_start_pnt => fpl_intm_start
        idx = 0
        do
          idx = idx+1
          write(luout,*) 'possibility # ',idx
          call prt_contr2(luout,fpl_intm_start_pnt%item%contr,
     &         op_info)
          if (.not.associated(fpl_intm_start_pnt%next)) exit
          fpl_intm_start_pnt => fpl_intm_start_pnt%next
        end do
      end if

      return
      end
