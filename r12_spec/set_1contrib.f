*----------------------------------------------------------------------*
      subroutine set_1contrib(flist,fac,idx,
     &     idx_intm,idx_op,nop,op_info)
*----------------------------------------------------------------------*
*     set terms arising from 1 in Q = 1 - P
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
      real(8), intent(in) ::
     &     fac
      integer, intent(in) ::
     &     nop,idx,idx_intm,idx_op(nop)
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     flist_pnt

      integer ::
     &     njoined_intm, njoined_op

      if (idx.gt.nop) then
        write(luout,*) 'idx, nop: ',idx,nop
        call quit(1,'set_1contrib','not enough operators on input list')
      end if

      njoined_intm = op_info%op_arr(idx_intm)%op%njoined
      njoined_op   = op_info%op_arr(idx_op(idx))%op%njoined

      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      if (njoined_intm.eq.njoined_op) then
        ! add blocks of AB
        call set_primitive_formula(flist_pnt,idx_op(idx),
     &       fac,idx_intm,.false.,op_info) 
      else if (njoined_intm.lt.njoined_op) then
        ! generate the appropriate self-contraction
        call expand_op_product2(flist_pnt,idx_intm,
     &     1d0,4,2,
     &     (/idx_intm,idx_op(idx),idx_op(idx),idx_intm/),
     &     (/1       ,2   , 2, 1       /),       
     &     -1, -1,
     &     0,0,
     &     0,0,
     &     0,0,
     &     op_info)
      else
        call quit(1,'set_1contrib','not prepared for this case!')        
      end if

      if (ntest.ge.100) then
        write (luout,*) 'generated in set_1contrib:'
        call print_form_list(luout,flist_pnt,op_info)
      end if

      return
      end

      
