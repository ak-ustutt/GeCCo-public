*----------------------------------------------------------------------*
      subroutine set_RC_contrib(flist,ansatz,approx,
     &     irdag,icint,
     &     idx_intm,idx_op,nop,op_info)
*----------------------------------------------------------------------*
*     set contributions to B matrix arising from
*       R12 C 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(inout), target ::
     &     flist
      character*(*) ::
     &     approx
      integer, intent(in) ::
     &     ansatz,nop,irdag,icint,
     &     idx_intm,idx_op(nop)
      type(operator_info), intent(inout) ::
     &     op_info

      type(formula_item), pointer ::
     &     flist_pnt

      integer ::
     &     idx_1, idx_2, idx_prj, njoined_c, njoined_intm

      if (irdag.gt.nop.or.
     &    icint.gt.nop) then
        write(lulog,*) 'idx: ',irdag,icint
        write(lulog,*) 'nop: ',nop
        call quit(1,'set_RC_contrib',
     &         'not enough operators on input list')
      end if

      if (idx_op(irdag).le.0.or.
     &    idx_op(icint).le.0) then
        write(lulog,*) 'idx: ',idx_op(irdag),idx_op(icint)
        call quit(1,'set_RC_contrib',
     &         'operator(s) not on input list')
      end if

      njoined_c = op_info%op_arr(idx_op(icint))%op%njoined
      if (njoined_c.ne.1)
     &     call quit(1,'set_RC_contrib','not yet prepared for this')
      njoined_intm = op_info%op_arr(idx_intm)%op%njoined

      !----------------------------------!
      ! - R12+ C                    !
      !----------------------------------!
      idx_1 = idx_op(irdag)
      idx_2 = idx_op(icint)

      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do

      ! should automatically give the correct connection (only 1)
      if (njoined_intm.eq.1) then
        call expand_op_product2(flist_pnt,idx_intm,
     &       -1d0,4,3,
     &       (/idx_intm,-idx_1,idx_2,idx_intm/),
     &       (/1       ,2     ,3    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       0,0,
     &       (/2,3,2,7/),1,  ! force P1P2 contraction
     &       .false.,op_info)
      else if (njoined_intm.eq.2) then
        call expand_op_product2(flist_pnt,idx_intm,
     &       -1d0,6,3,
     &       (/idx_intm,-idx_1,idx_intm,idx_intm,idx_2,idx_intm/),
     &       (/1       ,2     ,1       ,1       ,3    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       (/2,6, 1,5/),2,    ! avoid cross contrib. to external lines
     &       0,0,
     &       .false.,op_info)
      else
        call quit(1,'set_RC_contrib','unexpected: njoined_intm>2')
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'R^+.C contribution'
        call print_form_list(lulog,flist_pnt,op_info)
      end if

      return
      end

      
