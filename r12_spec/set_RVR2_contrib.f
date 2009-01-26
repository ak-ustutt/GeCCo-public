*----------------------------------------------------------------------*
      subroutine set_RVR2_contrib(flist,ansatz,approx,
     &     irint,ivint,
     &     idx_intm,idx_op,nop,op_info)
*----------------------------------------------------------------------*
*     set contributions to B matrix arising from
*       R12^+ V12 (v1+v2) Q12 R12 (analoguous to RC contribution,
*       if no explicit C intermediate is defined
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
     &     ansatz,nop,irint,ivint,
     &     idx_intm,idx_op(nop)
      type(operator_info), intent(inout) ::
     &     op_info

      type(formula_item), pointer ::
     &     flist_pnt

      integer ::
     &     idx_1, idx_2, idx_prj, njoined_c, njoined_intm

      if (irint.gt.nop.or.
     &    ivint.gt.nop) then
        write(luout,*) 'idx: ',irint,ivint
        write(luout,*) 'nop: ',nop
        call quit(1,'set_RVR2_contrib',
     &         'not enough operators on input list')
      end if

      if (idx_op(irint).le.0.or.
     &    idx_op(ivint).le.0) then
        write(luout,*) 'idx: ',idx_op(irint),idx_op(ivint)
        call quit(1,'set_RVR2_contrib',
     &         'operator(s) not on input list')
      end if

      njoined_intm = op_info%op_arr(idx_intm)%op%njoined

      !----------------------------------!
      ! - R12+ C                    !
      !----------------------------------!
      idx_1 = idx_op(irint)
      idx_2 = idx_op(ivint)

      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do

      ! should automatically give the correct connection (only 1)
      if (njoined_intm.eq.1) then
        call expand_op_product2(flist_pnt,idx_intm,
     &       -1d0,5,4,
     &       (/idx_intm,-idx_1,idx_2,idx_1,idx_intm/),
     &       (/1       ,2     ,3    ,4    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       0,0,
     &       (/2,3,1,IPART,      ! force P contraction
     &         2,4,1,IPART,      ! force P contraction
     &         3,4,1,IEXTR/),3,  ! force X contraction
     &       op_info)
      else if (njoined_intm.eq.2) then
        call quit(1,'set_RVR2_contrib','unchecked route')
        call expand_op_product2(flist_pnt,idx_intm,
     &       -1d0,7,4,
     &       (/idx_intm,-idx_1,idx_intm,idx_intm,idx_2,idx_1,idx_intm/),
     &       (/1       ,2     ,1       ,1       ,3     ,4   ,1       /),       
     &       -1, -1,
     &       0,0,
     &       (/2,7, 1,6/),2,    ! avoid cross contrib. to external lines
     &       0,0,
     &       op_info)
      else
        call quit(1,'set_RVR2_contrib','unexpected: njoined_intm>2')
      end if

      if (ntest.ge.100) then
        write(luout,*) 'R^+.V.R (2) contribution'
        call print_form_list(luout,flist_pnt,op_info)
      end if

      return
      end

      
