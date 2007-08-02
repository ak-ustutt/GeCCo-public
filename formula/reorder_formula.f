*----------------------------------------------------------------------*
      subroutine reorder_formula(form,op_info)
*----------------------------------------------------------------------*
*     reorder formula
*     a) group all contributions to one result block
*     b) ... further reodering .... ?
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

      type(formula_item), intent(inout), target ::
     &     form
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idx_tgt, iblk_tgt, nblk, nfound
      type(formula_item), pointer ::
     &     form_pnt
      type(formula_item_list), target ::
     &     fpl_reo
      type(formula_item_list), pointer ::
     &     fpl_reo_pnt

      if (ntest.ge.100) then
        write(luout,*) '========================='
        write(luout,*) ' reorder_formula at work'
        write(luout,*) '========================='
      end if

      if (form%command.eq.command_end_of_formula) return
      idx_tgt = form%target

      if (idx_tgt.gt.0) then
        nblk = op_info%op_arr(idx_tgt)%op%n_occ_cls
      else
        nblk = 1
      end if

      if (nblk.gt.1) then

        call init_formula_plist(fpl_reo)
        fpl_reo_pnt => fpl_reo
        ! keep [INIT]
        if (form%command.ne.command_add_contribution) then
          fpl_reo_pnt%item => form
        end if

        tgt_loop: do iblk_tgt = 1, nblk
          ! look for first term which contributes to this block
          form_pnt => form
          if (form_pnt%command.ne.command_add_contribution) then
            form_pnt => form_pnt%next
          end if
          search: do
            if (form_pnt%command.ne.command_add_contribution)
     &         cycle tgt_loop
            if (form_pnt%target.ne.idx_tgt .or.
     &           form_pnt%contr%idx_res.ne.idx_tgt) then
              write(luout,*) 'idx_tgt, target, idx_res:',
     &             idx_tgt, form_pnt%target, form_pnt%contr%idx_res
              call quit(1,'reorder_formula',
     &             'suspicious change of target')
            end if
            if (form_pnt%contr%iblk_res.eq.iblk_tgt) exit search
            form_pnt => form_pnt%next
          end do search
          ! get a new enty on pointer list
          if (associated(fpl_reo_pnt%item)) then
            call new_formula_plist_entry(fpl_reo_pnt)
            fpl_reo_pnt => fpl_reo_pnt%next
          end if
          ! collect pointers to all contributions with same block
          call collect_contr2block(fpl_reo_pnt,nfound,form_pnt,op_info)
          ! advance pointer
          do while(associated(fpl_reo_pnt%next))
            fpl_reo_pnt => fpl_reo_pnt%next
          end do
          if (ntest.ge.100) then
            write(luout,*) 'iblk, # terms: ',iblk_tgt,nfound
          end if
        end do tgt_loop
        
        call relink_formula_list(form,fpl_reo)

        call dealloc_formula_plist(fpl_reo)        

        if (ntest.ge.100) then
          call print_form_list(luout,form,op_info)
        end if

      end if

      return
      end

