*----------------------------------------------------------------------*
      subroutine set_Ycontrib(flist,ansatz,approx,
     &     irdag,irtilde,
     &     idx_intm,idx_op,nop,op_info)
*----------------------------------------------------------------------*
*     set Y type contributions to B matrix
*       R12 Q12 Rtilde
*     skip for HY1
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
     &     ansatz,nop,irdag,irtilde,
     &     idx_intm,idx_op(nop)
      type(operator_info), intent(inout) ::
     &     op_info

      type(formula_item), pointer ::
     &     flist_pnt, flist_pnt0

      integer ::
     &     idx_1, idx_2, idx_prj, njoined_intm


      if (approx(1:1).eq.'A') return
      if ( approx(8:10).eq.'HY1' .or.
     &    (approx(8:10).eq.'HY2'.and.ansatz.eq.1)) return
      if (approx(8:10).ne.'HY2'.and.approx(8:10).ne.'   ')
     &     call quit(1,'set_Ycontrib',
     &     'Hybrid approximations: only HY1/HY2 is implemented')

      if (irdag.gt.nop.or.
     &    irtilde.gt.nop) then
        write(lulog,*) 'idx: ',irdag,irtilde
        write(lulog,*) 'nop: ',nop
        call quit(1,'set_Ycontrib',
     &         'not enough operators on input list')
      end if

      if (idx_op(irdag).le.0.or.
     &    idx_op(irtilde).le.0) then
        write(lulog,*) 'idx: ',idx_op(irdag),idx_op(irtilde)
        call quit(1,'set_Ycontrib',
     &         'operator(s) not on input list')
      end if

      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      flist_pnt0 => flist

      njoined_intm = op_info%op_arr(idx_intm)%op%njoined

      !----------------------------------!
      ! - R12+ Rtilde                    !
      !----------------------------------!
      idx_1 = idx_op(irdag)
      idx_2 = idx_op(irtilde)

      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! set Q projector
      idx_prj = 1
      if (ansatz.gt.1) idx_prj = 3
      if (approx(8:10).eq.'HY2') idx_prj = 11
      if (njoined_intm.eq.1) then
c        call expand_op_product2(flist_pnt,idx_intm,
c     &       -1d0,4,3,
c     &       (/idx_intm,-idx_1,idx_2,idx_intm/),
c     &       (/1       ,2     ,3    ,1       /),       
c     &       -1, -1,
c     &       (/2,3/),1,   ! force additional "normal" contractions
c     &       0,0,
c     &       (/2,3,2,idx_prj/),1, ! def. of projector
c     &       .false.,op_info)
c        flist_pnt => flist
c        do while(associated(flist_pnt%next))
c          flist_pnt => flist_pnt%next
c        end do
        ! this gives the terms with ONLY the projector as contraction
        call expand_op_product2(flist_pnt,idx_intm,
     &       -1d0,4,3,
     &       (/idx_intm,-idx_1,idx_2,idx_intm/),
     &       (/1       ,2     ,3    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       0,0,
     &       (/2,3,2,idx_prj/),1, ! def. of projector
     &       .false.,op_info)
      else if (njoined_intm.eq.2) then
        call expand_op_product2(flist_pnt,idx_intm,
     &       -1d0,6,3,
     &       (/idx_intm,-idx_1,idx_intm,idx_intm,idx_2,idx_intm/),
     &       (/1       ,2     ,1       ,1       ,3    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       (/2,6, 1,5/),2,    ! avoid cross contrib. to external lines
     &       (/2,5,2,idx_prj/),1, ! def. of projector
     &       .false.,op_info)
      else
          call quit(1,'set_Ycontrib','unexpected: njoined_intm>2')
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'Ytilde contribution'
        call print_form_list(lulog,flist_pnt0,op_info)
      end if

      return
      end

      
