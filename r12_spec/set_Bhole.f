*----------------------------------------------------------------------*
      subroutine set_Bhole(flist,ansatz,
     &     iop0,iopa,iopb,ioph,
     &     idx_intm,idx_op,nop,op_info,orb_info)
*----------------------------------------------------------------------*
*     set hole B-matrix arising in "fixed-amplitude" calculations
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00
      character(6), parameter ::
     &     op_scr_f  = '_SCR_F'

      type(formula_item), intent(inout), target ::
     &     flist
      integer, intent(in) ::
     &     ansatz,nop,idx_intm,idx_op(nop),iop0,iopa,iopb,ioph
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(formula_item), pointer ::
     &     flist_pnt, flist_pnt0

      integer ::
     &     idx_0, idx_1, idx_2, idx_f, idx_h, idx_prj,
     &     njoined_intm, ndef
      type(operator), pointer ::
     &     opf_pnt, op_pnt
      integer ::
     &     occ_def(ngastp,2)
      integer, external ::
     &     idx_oplist2

      if (iopa.gt.nop.or.iopb.gt.nop.or.
     &    iop0.gt.nop.or.ioph.gt.nop) then
        write(luout,*) 'idx, nop: ',iop0,iopa,iopb,ioph,nop
        call quit(1,'set_Bhole','not enough operators on input list')
      end if

      ! dummy operator: 1 particle part of H
      call add_operator(op_scr_f,op_info)
      idx_f = idx_oplist2(op_scr_f,op_info)
      opf_pnt => op_info%op_arr(idx_f)%op
      ndef = 1
      occ_def = 0
      occ_def(IHOLE,1) = 1
      occ_def(IHOLE,2) = 1
      call set_uop(opf_pnt,op_scr_f,.false.,
     &       occ_def,ndef,orb_info)

      idx_0 = idx_op(iop0)
      idx_1 = idx_op(iopa)
      idx_2 = idx_op(iopb)
      idx_h = idx_op(ioph)

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
        ! terms which need additional contraction on top of 
        ! that caused by the projector and the contraction to F
        call expand_op_product2(flist_pnt,idx_intm,
     &     1d0,5,3,
     &     (/idx_intm,idx_0,idx_f,idx_0,idx_intm/),
     &     (/1       ,2    ,3    ,2    ,1       /),       
     &     -1, -1,
     &     (/2,3, 2,4, 3,4/),3,
     &     0,0,
     &     0,0, 
     &     .false.,op_info)
        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        end do
        call expand_op_product2(flist_pnt,idx_intm,
     &     -1d0,5,4,
     &     (/idx_intm,-idx_1,idx_f,idx_2,idx_intm/),
     &     (/1       ,2     ,3    ,4    ,1       /),       
     &     -1, -1,
     &     (/2,3, 2,4, 3,4/),3,
     &     0,0,
     &     (/2,4,2,idx_prj/),1, ! def. of projector
     &     .false.,op_info)
        ! try terms that are contracted by the projector/F only
        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        end do
        call expand_op_product2(flist_pnt,idx_intm,
     &     1d0,5,3,
     &     (/idx_intm,idx_0,idx_f,idx_0,idx_intm/),
     &     (/1       ,2    ,3    ,2    ,1       /),       
     &     -1, -1,
     &     (/2,3, 3,4/),2,
     &     (/2,4/),1,
     &     0,0,
     &     .false.,op_info)
        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        end do
        call expand_op_product2(flist_pnt,idx_intm,
     &     -1d0,5,4,
     &     (/idx_intm,-idx_1,idx_f,idx_2,idx_intm/),
     &     (/1       ,2     ,4    ,3    ,1       /),       
     &     -1, -1,
     &     (/2,3, 3,4/),2,
     &     0,0,
     &     (/2,4,2,idx_prj/),1,            ! def. of projector
     &     .false.,op_info)
      else
        call quit(1,'set_Bhole','njoined>1 needed? not available yet!')
      end if

      ! F -> H replacement:
      op_pnt => op_info%op_arr(idx_h)%op
      call form_op_replace(opf_pnt%name,op_pnt%name,.false.,
     &                     flist_pnt0,op_info)
      if (ntest.ge.100) then
        write(luout,*) 'result for Bhole'
        call print_form_list(luout,flist_pnt0,op_info)
      end if

      ! get rid of scratch operator
      call del_operator(op_scr_f,op_info)

      return
      end

      
