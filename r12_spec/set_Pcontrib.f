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
     &     flist_pnt, flist_pnt0

      integer ::
     &     idx_1, idx_2, idx_prj, njoined_intm

      if (iopa.gt.nop.or.iopb.gt.nop) then
        write(lulog,*) 'idx, nop: ',iopa,iopb,nop
        call quit(1,'set_Pcontrib','not enough operators on input list')
      end if

      idx_1 = idx_op(iopa)
      idx_2 = idx_op(iopb)

      njoined_intm = op_info%op_arr(idx_intm)%op%njoined

      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      flist_pnt0 => flist_pnt
      idx_prj = 2
      if (ansatz.gt.1) idx_prj = 4
      if (njoined_intm.eq.1) then
c        ! terms which need additional contraction on top of 
c        ! that caused by the projector
c        call expand_op_product2(flist_pnt,idx_intm,
c     &     -1d0,4,3,
c     &     (/idx_intm,-idx_1,idx_2,idx_intm/),
c     &     (/1       ,2   ,3    ,1       /),       
c     &     -1, -1,
c     &     (/2,3/),1,
c     &     0,0,
c     &     (/2,3,2,idx_prj/),1, ! def. of projector
c     &     .false.,op_info)
c        flist_pnt => flist
c        do while(associated(flist_pnt%next))
c          flist_pnt => flist_pnt%next
c        end do
        ! try terms that are contracted by the projector only
        call expand_op_product2(flist_pnt,idx_intm,
     &     -1d0,4,3,
     &     (/idx_intm,-idx_1,idx_2,idx_intm/),
     &     (/1       ,2   ,3    ,1       /),       
     &     -1, -1,
     &     0,0,
     &     0,0,
     &     (/2,3,2,idx_prj/),1, ! def. of projector
     &     .false.,op_info)
      else if (njoined_intm.eq.2) then
        call expand_op_product2(flist_pnt,idx_intm,
     &     -1d0,6,3,
     &     (/idx_intm,-idx_1,idx_intm,idx_intm,idx_2,idx_intm/),
     &     (/1       ,2   ,1       ,1       ,3    ,1       /),       
     &     -1, -1,
     &     0,0,
     &     (/2,6, 1,5/),2,      ! avoid cross contrib. to external lines
     &     (/2,5,2,idx_prj/),1, ! def. of projector
     &     .false.,op_info)
      else
        call quit(1,'set_Pcontrib','unexpected: njoined_intm>2')
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'result for P contribution'
        call print_form_list(lulog,flist_pnt0,op_info)
      end if

      return
      end

      
