*----------------------------------------------------------------------*
      subroutine cc_form_post(form,nterms,idxtbar,idxham,idxtop,op_info)
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
      include 'explicit.h'

      type(formula_item), intent(inout), target ::
     &     form
      integer, intent(in) ::
     &     idxtbar, idxham, idxtop
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(out) ::
     &     nterms
      
      integer ::
     &     n_commu(0:4)
      type(formula_item_list), target ::
     &     fpl_reo
      type(formula_item_list), pointer ::
     &     fpl_pnt
      type(formula_item), pointer ::
     &     form_pnt
      type(operator), pointer ::
     &     op_tbar
      integer ::
     &     nblk_tbar, iblk

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

      op_tbar => op_info%op_arr(idxtbar)%op
      nblk_tbar = op_tbar%n_occ_cls

      if (iprlvl.gt.0) then
        write(luout,'(2x,42("-"))')
        write(luout,'(3x,a)') '  L    number of   n-fold commutators'
        write(luout,'(3x,a)') 'class    terms    0    1    2    3    4'
        write(luout,'(2x,42("-"))')
      end if

      ! loop over Tbar blocks (0 is the unit-projection)
      do iblk = 0, nblk_tbar
        if (iblk.eq.0) then
          call collect_terms_w_op(fpl_pnt,form_pnt,idxtbar,-1,0)          
          fpl_pnt => fpl_pnt%next
        else
          call collect_terms_w_op(fpl_pnt,form_pnt,idxtbar,iblk,1)
          fpl_pnt => fpl_pnt%next
        end if
c dbg
c        print *,'iblk = ',iblk
c        print *,'first collected term:'
c        print *,' command = ',fpl_pnt%item%command
c        if (fpl_pnt%item%command.eq.4) then
c          call prt_contr2(luout,fpl_pnt%item%contr,op_info%op_arr)
c        end if
c dbg
        if(explicit)then
          call r12_count_terms(fpl_pnt,idxtop,n_commu,op_info%op_arr)
        else  
          call cc_count_terms(fpl_pnt,idxtop,n_commu)
        endif  

        do while(associated(fpl_pnt%next))
          fpl_pnt => fpl_pnt%next
        end do

        if (iprlvl.ge.10) then
          write(luout,'(2x,"??",i4,3x,i6,x,5(x,i4))')
     &       iblk,sum(n_commu(0:4)),n_commu(0:4)
        else if (iprlvl.gt.0) then
          write(luout,'(4x,i4,3x,i6,x,5(x,i4))')
     &       iblk,sum(n_commu(0:4)),n_commu(0:4)
        end if
        nterms = nterms + sum(n_commu(0:4))

      end do

      call relink_formula_list(form,fpl_reo)

      call dealloc_formula_plist(fpl_reo)

      return
      end
