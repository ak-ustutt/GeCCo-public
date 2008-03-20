*----------------------------------------------------------------------*
      subroutine set_C_fr(flist,approx,
     &     iham,ir12,irbar,irtilde,
     &     idx_intm,idx_op,nop,op_info,orb_info)
*----------------------------------------------------------------------*
*     set contributions to C matrix which arise form R12 and Fock matrix
*
*     approx A/B: -f^c_a r^{ij}_{cb} + rbar^{ij}_{ab} - rtilde^{ij}_{ab} 
*
*     approx C:   f^{all}_{a} r^{ij}_{a,all}
*
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
      character*(*) ::
     &     approx
      integer, intent(in) ::
     &     nop,ir12,irbar,irtilde,iham,
     &     idx_intm,idx_op(nop)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(formula_item), pointer ::
     &     flist_pnt

      real(8) ::
     &     fac
      integer ::
     &     idx_1, idx_2, idx_prj, idx_f, ndef
      type(operator), pointer ::
     &     opf_pnt, op_pnt
      integer ::
     &     occ_def(ngastp,2)

      integer, external ::
     &     idx_oplist2


      if (iham.gt.nop.or.
     &    ir12.gt.nop) then
          write(luout,*) 'idx: ',iham,ir12
          write(luout,*) 'nop: ',nop
          call quit(1,'set_C_fr',
     &         'not enough operators on input list')
      end if

      if (idx_op(iham).le.0.or.
     &    idx_op(ir12).le.0) then
        write(luout,*) 'idx: ',idx_op(iham),idx_op(ir12)
        call quit(1,'set_C_fr',
     &         'operator(s) not on input list')
      end if

      ! dummy operator: 1 particle part of H
      call add_operator(op_scr_f,op_info)
      idx_f = idx_oplist2(op_scr_f,op_info)
      opf_pnt => op_info%op_arr(idx_f)%op
      ndef = 1
      occ_def = 0
      occ_def(1,1:2) = 1
      call set_uop(opf_pnt,op_scr_f,.false.,
     &     occ_def,ndef,orb_info)

      idx_1 = idx_op(ir12)
      idx_2 = idx_op(iham)

      if (approx(1:1).eq.'C') then
        fac = 1d0
        idx_prj = 0
      else
        fac = -1d0
        idx_prj = IPART
      end if

      ! go to end of list
      flist_pnt => flist
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      call expand_op_product2(flist_pnt,idx_intm,
     &       fac,4,3,
     &       (/idx_intm,idx_f,idx_1,idx_intm/),
     &       (/1       ,2    ,3    ,1       /),       
     &       -1, -1,
     &       0,0,
     &       0,0,  
     &       (/2,3,1,idx_prj/),1,
     &       op_info)

      ! F -> H replacement:
      op_pnt => op_info%op_arr(idx_2)%op
      call form_op_replace(opf_pnt%name,op_pnt%name,flist_pnt,op_info)

      if (ntest.ge.100) then
        write(luout,*) 'approx = ',trim(approx)
        write(luout,*) 'C: F.R12 contribution:'
        call print_form_list(luout,flist_pnt,op_info)
      end if

      call del_operator(op_scr_f,op_info)

      if (approx(1:1).eq.'C') return

      if (irbar.gt.nop.or.
     &    irtilde.gt.nop) then
          write(luout,*) 'idx: ',irbar,irtilde
          write(luout,*) 'nop: ',nop
          call quit(1,'set_C_fr',
     &         'not enough operators on input list')
      end if

      if (idx_op(irbar).le.0.or.
     &    idx_op(irtilde).le.0) then
        write(luout,*) 'idx: ',idx_op(irbar),idx_op(irtilde)
        call quit(1,'set_C_fr',
     &         'operator(s) not on input list')
      end if

      ! add Rbar and Rtilde in case of the other approximations

      idx_1 = idx_op(irbar)
      idx_2 = idx_op(irtilde)

      ! advance to end of list
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! add blocks of Rbar
      call set_primitive_formula(flist_pnt,idx_1,
     &       1d0,idx_intm,.false.,op_info) 
      
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! add blocks of Rtilde
      call set_primitive_formula(flist_pnt,idx_2,
     &       1d0,idx_intm,.false.,op_info) 
      

      return
      end

      
