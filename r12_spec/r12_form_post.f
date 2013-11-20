*----------------------------------------------------------------------*
      subroutine r12_form_post(form,nterms,
     &     idxtbar,idxcbar,idxham,idxtop,idxc12, iprint,
     &     op_info)
*----------------------------------------------------------------------*
*     resort terms according to tbar blocks and count terms
*     print some statistics
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'

      type(formula_item), intent(inout), target ::
     &     form
      integer, intent(in) ::
     &     idxtbar, idxcbar, idxham, idxtop, idxc12, iprint
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(out) ::
     &     nterms
      
      integer ::
     &     n_commu(1:15)
      type(formula_item_list), target ::
     &     fpl_reo
      type(formula_item_list), pointer ::
     &     fpl_pnt
      type(formula_item), pointer ::
     &     form_pnt
      type(operator), pointer ::
     &     op_tbar, op_cbar
      integer ::
     &     nblk_tbar, nblk_cbar, iblk

      nterms = 0
      ! initialize reodering list
      call init_formula_plist(fpl_reo)
      fpl_pnt => fpl_reo
      ! point to "initialize" item
      if (form%command.ne.command_set_target_init)
     &     call quit(1,'cc_form_post','wrong assumptions? (1)')
      fpl_pnt%item => form
      ! and advance to first "add"
      form_pnt => form%next
      if (form_pnt%command.ne.command_add_contribution)
     &     call quit(1,'cc_form_post','wrong assumptions? (2)')

      op_tbar => op_info%op_arr(abs(idxtbar))%op
      nblk_tbar = op_tbar%n_occ_cls
      op_cbar => op_info%op_arr(abs(idxcbar))%op
      nblk_cbar = op_cbar%n_occ_cls

      if (iprint.gt.0) then
        write(lulog,'(x,78("-"))')
        write(lulog,'(x,a)')
     &  ' L    number  n-fold commutators'
        write(lulog,'(x,a)')
     &  'class  of      0   1       2           3                4'
        write(lulog,'(x,a)')
     &  '      terms        T   C  TT  TC  CC TTT TTC TCC CCC TTTT '//
     &       'TTTC TTCC TCCC CCCC'
        write(lulog,'(x,78("-"))')
      end if

      ! loop over Tbar blocks (0 is the unit-projection)
      do iblk = 0, nblk_tbar
        if (iblk.eq.0) then
          call collect_terms_w_op(fpl_pnt,form_pnt,
     &         2,(/idxtbar,idxcbar/),-1,0)          
          if (associated(fpl_pnt%next)) fpl_pnt => fpl_pnt%next
        else
          call collect_terms_w_op(fpl_pnt,form_pnt,1,idxtbar,iblk,1)
          if (associated(fpl_pnt%next)) fpl_pnt => fpl_pnt%next
        end if

        call r12_count_terms(fpl_pnt,idxtop,idxc12,n_commu,op_info)

        do while(associated(fpl_pnt).and.associated(fpl_pnt%next))
          fpl_pnt => fpl_pnt%next
        end do

        if (iprint.ge.10) then
c          write(lulog,'("??",x,"L",i2,3x,i6,x,15(i4))')
          write(lulog,'("??",x,"L",i1,x,i6,x,10(i4),5(i5))')
     &       iblk,sum(n_commu(1:15)),n_commu(1:15)
        else if (iprint.gt.0) then
          write(lulog,'(3x,"L",i1,x,i6,x,10(i4),5(i5))')
     &       iblk,sum(n_commu(1:15)),n_commu(1:15)
        end if
        nterms = nterms + sum(n_commu(1:15))

      end do
      ! loop over Cbar blocks
      do iblk = 1, nblk_cbar
        call collect_terms_w_op(fpl_pnt,form_pnt,1,idxcbar,iblk,1)
        if (associated(fpl_pnt%next)) then
          fpl_pnt => fpl_pnt%next
          call r12_count_terms(fpl_pnt,idxtop,idxc12,n_commu,op_info)
        else
          n_commu = 0
        end if

        do while(associated(fpl_pnt%next))
          fpl_pnt => fpl_pnt%next
        end do

        if (iprint.ge.10) then
          write(lulog,'("??",x,"C",i1,x,i6,x,10(i4),5(i5))')
c          write(lulog,'("??",x,"C",i2,3x,i6,x,15(i4))')
     &       iblk+1,sum(n_commu(1:15)),n_commu(1:15)
        else if (iprint.gt.0) then
          write(lulog,'(3x,"C",i1,x,i6,x,10(i4),5(i5))')
     &       iblk+1,sum(n_commu(1:15)),n_commu(1:15)
        end if
        nterms = nterms + sum(n_commu(1:15))

      end do

      if (iprint.gt.0)
     &     write(lulog,'(x,78("-"))')

      call relink_formula_list(form,fpl_reo)

      call dealloc_formula_plist(fpl_reo)

      return
      end
