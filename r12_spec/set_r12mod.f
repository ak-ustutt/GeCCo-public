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

      logical ::
     &     new
      integer ::
     &     idx_1, idx_2, nproj, iproj(12)

      if ((trim(type).ne.'RB'.and.
     &     trim(type).ne.'RD'.and.
     &     trim(type).ne.'RT'.and.
     &     trim(type).ne.'RV').or.
     &     nop.lt.3) then
        write(luout,*) 'type, nop: ',trim(type),nop
        call quit(1,'set_r12mod','misguided call?')
      end if

      new = .false.

      if (trim(type).eq.'RB') then
        idx_1 = idx_op(1) ! R12
        idx_2 = idx_op(3) ! one-int
      else
        idx_1 = idx_op(3) ! one-int
        idx_2 = idx_op(1) ! R12
      end if

      if (trim(type).eq.'RD') then 
        iproj(1:4) = (/2,3,1,0/) ! full one-particle RI
        nproj = 1
      end if
      if (trim(type).eq.'RV') then
        new = .true.  ! use new descriptor based setup
        iproj(5:8) = (/2,3,1,IPART/) ! only P contraction
                                    ! (and X contraction, see below)
        iproj(1:4) = (/2,4,1,-IHOLE/) ! only H contribution to open line
        iproj(9:12)= (/1,3,1,IEXTR/) ! only X contribution to open line
        nproj = 3
      end if

      if (new) goto 2308
      ! ------------------------------------
      ! contributions from usual R12-list
      ! ------------------------------------
      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! set the product
c dbg
      print *,'relevant call'
c dbg
      call expand_op_product2(flist_pnt,idx_intm,
     &       1d0,4,3,
     &       (/idx_intm,idx_1,idx_2,idx_intm/),
     &       (/1       ,2    ,3    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       0,0,
     &       iproj,nproj,
     &       .false.,op_info)

      if (trim(type).eq.'RV') then

        ! add X contraction contributions
        iproj(1:4) = (/2,4,1,-IHOLE/) ! only H contribution to open line
        iproj(5:8) = (/2,3,1,IEXTR/)
        nproj = 2
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
     &       iproj,nproj,
     &       .false.,op_info)
      end if

      if (trim(type).eq.'RB') then
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
     &       (/2,3,1,0/),1,     ! one-particle RI
     &       .false.,op_info)
      
      end if

      if (ntest.ge.100) then
        write(luout,*) 'type: ',trim(type)
        write(luout,*) 'result after expand_op_product'
        call print_form_list(luout,flist,op_info)
      end if

      return

 2308 continue ! new part

      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! contraction over P, restrict R12 to PX part
      call expand_op_product3(flist_pnt,idx_intm,
     &       1d0,4,3,
     &       (/idx_intm,idx_1,idx_2,idx_intm/),
     &       (/1       ,2    ,3    ,1       /),
     &       -1, -1,
     &       0,0,
     &       0,0,
     &       (/'2,3,,P            ',
     &         '3,,PX,[HPVX][HPVX]'/),2,
     &       op_info)

      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! contraction over X, restrict R12 to PX and XX parts
      call expand_op_product3(flist_pnt,idx_intm,
     &       1d0,4,3,
     &       (/idx_intm,idx_1,idx_2,idx_intm/),
     &       (/1       ,2    ,3    ,1       /),
     &       -1, -1,
     &       0,0,
     &       0,0,
     &       (/'2,3,,X               ',
     &         '3,,[PX]X,[HPVX][HPVX]'/),2,
     &       op_info)


      if (ntest.ge.100) then
        write(luout,*) 'type: ',trim(type)
        write(luout,*) 'result after expand_op_product'
        call print_form_list(luout,flist,op_info)
      end if

      end

      
