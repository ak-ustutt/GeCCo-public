      subroutine set_pint_contract(flist,ansatz,
     &     idx_opsin,nopsin,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     Routine used to form the various contractions necessary to form
*     the CABS approximation to the P-intermediate.
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_orbinf.h'

      type(formula_item), target, intent(inout) ::
     &     flist
      integer, intent(in) ::
     &     ansatz, nopsin, idx_opsin(nopsin)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info


      integer ::
     &     idx, idx_shape, idx_temp, ndef, idx_prj
      character ::
     &     name*(form_maxlen_label*2)

      type(operator), pointer ::
     &     op
      type(formula_item), pointer ::
     &     form_pnt

      integer, external ::
     &     idx_oplist2


      if(ntest.ge.100)then
        write(luout,*)'P-Intermediate Contraction'
        write(luout,*)'Constituent operators: '
        do idx = 1, nopsin
          write(luout,*)trim(op_info%op_arr(idx_opsin(idx))%op%name)
        enddo
      endif

      ! Get indices of input operators.
      idx_shape = idx_opsin(1)

      ! Point to the formula and move to the end of the list.
      form_pnt => flist
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F2G part:
      call set_1contrib(form_pnt,1d0,6,
     &     idx_shape,idx_opsin,nopsin,op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the (FG) contracted with F part.
      idx_prj = 8
      call expand_op_product2(form_pnt,idx_shape,
     &     -2d0,7,3,
     &     (/idx_shape,idx_opsin(4),idx_opsin(4),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,           2,           2,
     &               1,        1,           3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/3,6,2,idx_prj/),1,
     &     op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{pq}.G_{pq}^{rs}.R_{rs}^{ij}.
      idx_prj = 8
      call expand_op_product2(form_pnt,idx_shape,
     &     1d0,7,4,
     &     (/idx_shape,-idx_opsin(2),idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,3,2,idx_prj,3,6,2,idx_prj/),2,
     &     op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{kl}^{p"m}.V_{p"m}^{ij}.
      idx_prj = 9
      if(ansatz.gt.1) idx_prj = 10
      call expand_op_product2(form_pnt,idx_shape,
     &     -2d0,7,3,
     &     (/idx_shape,-idx_opsin(2),idx_opsin(5),
     &       idx_shape,idx_shape,idx_opsin(5),idx_shape/),
     &     (/        1,            2,           3,
     &               1,        1,           3,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,3,2,idx_prj/),1,
     &     op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"q}.G_{p"q}^{r"s}.R_{r"s}^{ij}.
      idx_prj = 9
      if(ansatz.gt.1) idx_prj = 10
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,7,4,
     &     (/idx_shape,-idx_opsin(2),idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,3,2,idx_prj,3,6,2,idx_prj/),2,
     &     op_info)

      if(ntest.ge.100)then
        write(luout,*)'Final formula: P-Int.'
        call print_form_list(luout,flist,op_info)
      endif

      return
      end
