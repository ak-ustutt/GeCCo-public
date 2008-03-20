*----------------------------------------------------------------------*
      subroutine set_Pcontrib(flist,ansatz,
     &     iopa,iopb,
     &     idx_intm,idx_op,nop,op_info)
*----------------------------------------------------------------------*
*     set terms arising from P in Q = 1 - P
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
      integer, intent(in) ::
     &     ansatz,nop,idx_intm,idx_op(nop),iopa,iopb
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     flist_pnt

      integer ::
     &     idx_1, idx_2, idx_prj

      if (iopa.gt.nop.or.iopb.gt.nop) then
        write(luout,*) 'idx, nop: ',iopa,iopb,nop
        call quit(1,'set_Pcontrib','not enough operators on input list')
      end if

      idx_1 = idx_op(iopa)
      idx_2 = idx_op(iopb)

      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      idx_prj = 2
      if (ansatz.gt.1) idx_prj = 4
      call expand_op_product2(flist_pnt,idx_intm,
     &     -1d0,6,3,
     &     (/idx_intm,-idx_1,idx_intm,idx_intm,idx_2,idx_intm/),
     &     (/1       ,2   ,1       ,1       ,3    ,1       /),       
     &     -1, -1,
     &     0,0,
     &     (/2,6, 1,5/),2,      ! avoid cross contrib. to external lines
     &     (/2,5,2,idx_prj/),1, ! def. of projector
     &     op_info)

      if (ntest.ge.100) then
        write(luout,*) 'result after expand_op_product'
        call print_form_list(luout,flist,op_info)
      end if

      return
      end

      
