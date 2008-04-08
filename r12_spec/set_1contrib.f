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

      if (idx.gt.nop) then
        write(luout,*) 'idx, nop: ',idx,nop
        call quit(1,'set_1contrib','not enough operators on input list')
      end if

      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! add blocks of AB
      call set_primitive_formula(flist_pnt,idx_op(idx),
     &       fac,idx_intm,.false.,op_info) 

      return
      end

      
