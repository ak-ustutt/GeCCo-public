*----------------------------------------------------------------------*
      subroutine collect_contr2block(fpl_intm_c2blk,nterms,
     &                               fl_intm,op_info)
*----------------------------------------------------------------------*
*     starting at first item of fl_intm, collect all terms that 
*     contribute to the same block on pointer list fpl_intm_c2blk
*     return number of terms on nterms
*     on entry, the head of fpl_intm should be allocated
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), target, intent(in) ::
     &     fl_intm
      type(formula_item_list), target, intent(inout) ::
     &     fpl_intm_c2blk
      integer, intent(out) ::
     &     nterms
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     fl_pnt
      type(formula_item_list), pointer ::
     &     fpl_c2blk_pnt

      integer ::
     &     iterm, idxop_intm, iblk_intm

      if (ntest.ge.100) then
        write(luout,*) '==============================='
        write(luout,*) ' info from collect_contr2block'
        write(luout,*) '==============================='
      end if

      if (fl_intm%command.ne.command_add_contribution)
     &       call quit(1,'collect_contr2block',
     &     '[ADD] expected on first item')

      nterms = 1 ! at least the current item on fl_intm
c      nullify(fpl_intm_c2blk%prev)
      nullify(fpl_intm_c2blk%next)

      fpl_intm_c2blk%item => fl_intm

      idxop_intm = fl_intm%target
      iblk_intm  = fl_intm%contr%iblk_res

      if (ntest.ge.100) then
        write(luout,*) 'idxop_intm, iblk_intm: ',idxop_intm, iblk_intm
      end if

      fl_pnt => fl_intm
      fpl_c2blk_pnt => fpl_intm_c2blk
      do while (associated(fl_pnt%next))
        fl_pnt => fl_pnt%next
        if (fl_pnt%command.eq.command_end_of_formula) exit
        if (fl_pnt%command.eq.command_set_target_init) exit
        if (fl_pnt%command.ne.command_add_contribution)
     &       call quit(1,'collect_contr2block',
     &       'only [ADD],[INIT],[END] expected')
        if (fl_pnt%contr%idx_res.ne.idxop_intm) then
          write(luout,*) 'scanning for operator ',idxop_intm
          call prt_contr2(luout,fl_pnt%contr,op_info)
          call quit(1,'collect_contr2block',
     &       'suspicious change of operator index')
        end if
        if (fl_pnt%contr%iblk_res.eq.iblk_intm) then
          ! store pointer to this item on list
          nterms = nterms+1
          allocate(fpl_c2blk_pnt%next)
          fpl_c2blk_pnt%next%prev => fpl_c2blk_pnt
          fpl_c2blk_pnt => fpl_c2blk_pnt%next
          fpl_c2blk_pnt%item => fl_pnt
          nullify(fpl_c2blk_pnt%next)
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'nterms = ',nterms
        fpl_c2blk_pnt => fpl_intm_c2blk
        iterm = 1
        do
          write(luout,*) 'term #',iterm
          call prt_contr2(luout,fpl_c2blk_pnt%item%contr,op_info)
          if (.not.associated(fpl_c2blk_pnt%next)) exit
          fpl_c2blk_pnt => fpl_c2blk_pnt%next
          iterm = iterm+1
        end do
      end if

      return
      end
