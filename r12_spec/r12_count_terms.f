*----------------------------------------------------------------------*
      subroutine r12_count_terms(fpl_terms,idxtop,idxc12,
     &                           n_per_class,op_info)
*----------------------------------------------------------------------*
*     count number of commutators for terms pointed to by fpl_terms
*
*     1  T  TT  TTT  TTTT
*        C  TC  TTC  TTTC
*           CC  TCC  TTCC
*               CCC  TCCC
*                    CCCC   -> 15 classes
*
*     adr = n(T+C)*(n(T+C)+1)/2 + n(C)+1
*
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

      integer, intent(in) ::
     &     idxtop, idxc12
      type(formula_item_list), target, intent(inout) ::
     &     fpl_terms
      integer, intent(out) ::
     &     n_per_class(15)

      type(operator_info), intent(in) ::
     &     op_info
      type(formula_item_list), pointer ::
     &     fpl_pnt
      integer ::
     &     nfound_top, nfound_c12, icls, nvtx, ivtx
      type(contraction), pointer ::
     &     cur_contr

      n_per_class(1:15) = 0

      if (.not.associated(fpl_terms%item)) print *,'WARN'
      if (.not.associated(fpl_terms%item)) return

      fpl_pnt => fpl_terms

      do
        if (.not.associated(fpl_pnt).or.
     &      fpl_pnt%item%command.ne.command_add_contribution) exit

        cur_contr => fpl_pnt%item%contr
        nvtx = cur_contr%nvtx
        nfound_top = 0
        nfound_c12 = 0
        if(ntest.gt.0)then
          call prt_contr2(lulog,cur_contr,op_info)
        endif  
        
        do ivtx = 1, nvtx
          if (cur_contr%vertex(ivtx)%idx_op.eq.idxtop)
     &         nfound_top = nfound_top+1
        end do
        do ivtx = 1, nvtx
          if (cur_contr%vertex(ivtx)%idx_op.eq.idxc12)
     &         nfound_c12 = nfound_c12+1
        end do
        if (nfound_top.gt.4.or.nfound_c12.gt.4)
     &       call quit(1,'cc_count_terms','> 4-fold commutator?')
        icls = (nfound_top+nfound_c12)*(nfound_top+nfound_c12+1)/2
     &         + nfound_c12 + 1
        n_per_class(icls) = n_per_class(icls)+1 
        
        if (.not.associated(fpl_pnt%next)) exit
        fpl_pnt => fpl_pnt%next
      end do

      return
      end
