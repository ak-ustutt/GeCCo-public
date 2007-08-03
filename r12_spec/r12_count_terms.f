*----------------------------------------------------------------------*
      subroutine r12_count_terms(fpl_terms,idxtop,n_commu,ops)
*----------------------------------------------------------------------*
*     count number of commutators for terms pointed to by fpl_terms
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'
      include 'def_operator_array.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     idxtop
      type(formula_item_list), target, intent(inout) ::
     &     fpl_terms
      integer, intent(out) ::
     &     n_commu(0:4)

      type(operator_array), intent(in) ::
     &     ops(*)
      type(formula_item_list), pointer ::
     &     fpl_pnt
      integer ::
     &     nfound_top, nvtx, ivtx
      type(contraction), pointer ::
     &     cur_contr

      n_commu(0:4) = 0

      fpl_pnt => fpl_terms

      do
        if (fpl_pnt%item%command.ne.command_add_contribution) exit

        cur_contr => fpl_pnt%item%contr
        nvtx = cur_contr%nvtx
        nfound_top = 0
        if(ntest.gt.0)then
          call prt_contr2(luout,cur_contr,ops)
        endif  
        
        do ivtx = 1, nvtx
          if (cur_contr%vertex(ivtx)%idx_op.eq.idxtop)
     &         nfound_top = nfound_top+1
        end do
        if (nfound_top.gt.4)
     &       call quit(1,'cc_count_terms','> 4-fold commutator?')
        n_commu(nfound_top) = n_commu(nfound_top)+1 
        
        if (.not.associated(fpl_pnt%next)) exit
        fpl_pnt => fpl_pnt%next
      end do

      return
      end
