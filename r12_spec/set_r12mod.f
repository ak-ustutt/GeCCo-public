*----------------------------------------------------------------------*
      subroutine set_r12mod(flist,type,
     &     idx_intm,idx_op,nop,op_info)
*----------------------------------------------------------------------*
*     set formulae for modified R12-integral lists
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
      character*(*), intent(in) ::
     &     type
      integer, intent(in) ::
     &     nop,idx_intm,idx_op(nop)
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     flist_pnt

      integer ::
     &     idx_1, idx_2

      if (trim(type).ne.'RB'.or.
     &    trim(type).ne.'RT'.or.
     &    trim(type).ne.'RV'.or.
     &     nop.ne.3) then
        write(luout,*) 'type, nop: ',trim(type),nop
        call quit(1,'set_r12mod','misguided call?')
      end if

      if (trim(type).eq.'RB') then
        idx_1 = idx_op(1)
        idx_2 = idx_op(3)
      else
        idx_1 = idx_op(3)
        idx_2 = idx_op(1)
      end if

      ! ------------------------------------
      ! contributions from extended R12-list
      ! ------------------------------------
      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! set the product
      call expand_op_product2(flist_pnt,idx_intm,
     &       1d0,4,3,
     &       (/idx_intm,idx_1,idx_2,idx_intm/),
     &       (/1       ,2    ,3    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       0,0,
     &       (/2,3,1,0/),1,! one-particle RI
     &       op_info)

      ! ------------------------------------
      ! contributions from extended R12-list
      ! ------------------------------------
      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      if (trim(type).eq.'RB') then
        idx_1 = idx_op(2)
        idx_2 = idx_op(3)
      else
        idx_1 = idx_op(3)
        idx_2 = idx_op(2)
      end if
      ! set the product
      call expand_op_product2(flist_pnt,idx_intm,
     &       1d0,4,3,
     &       (/idx_intm,idx_1,idx_2,idx_intm/),
     &       (/1       ,2    ,3    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       0,0,
     &       (/2,3,1,0/),1,! one-particle RI
     &       op_info)
      
      if (ntest.ge.100) then
        write(luout,*) 'type: ',trim(type)
        write(luout,*) 'result after expand_op_product'
        call print_form_list(luout,flist,op_info)
      end if

      return
      end

      
