      subroutine set_zint_contract(flist,ansatz,
     &     idx_opsin,nopsin,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     Routine used to form the various contractions necessary to form
*     the CABS approximation to the Z-intermediate.
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 100

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
     &     idx, idx_shape, nops ,nvtx, idx_temp, ndef, idx_prj,
     &     idx_prj2, idx_prj3
      character ::
     &     name*(form_maxlen_label*2)

      integer, pointer ::
     &     idxarr(:), svtxarr(:), occ_def(:,:,:)
      character(4), parameter ::
     &     op_temp = 'TEMP'

      type(operator), pointer ::
     &     op
      type(formula_item), pointer ::
     &     form_pnt

      integer, external ::
     &     idx_oplist2


      if(ntest.ge.100)then
        write(luout,*)'Z-Intermediate Contraction'
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

      ! Add the G^{p'q}_{km}.FF_{p'l}^{ij} terms.
      idx_prj = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     1d0,9,3,
     &     (/idx_shape,idx_opsin(3),idx_shape,idx_shape,idx_opsin(4),
     &       idx_shape,idx_shape,idx_opsin(4),idx_shape/),
     &     (/        1,           2,        1,        1,           3,
     &               1,        1,           3,        1/),
     &     -1,-1,
     &     (/2,6/),1,
     &     0,0,
     &     (/2,5,1,idx_prj,3,5,1,idx_prj,2,4,1,idx_prj/),3,
     &     op_info)

      form_pnt => flist
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      idx_prj = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     1d0,9,3,
     &     (/idx_shape,idx_opsin(3),idx_shape,idx_shape,idx_opsin(4),
     &       idx_shape,idx_shape,idx_opsin(4),idx_shape/),
     &     (/        1,           2,        1,        1,           3,
     &               1,        1,           3,        1/),
     &     -1,-1,
     &     (/2,6/),1,
     &     0,0,
     &     (/2,5,1,idx_prj,3,5,1,idx_prj/),2,
     &     op_info)

      form_pnt => flist
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      do idx = 1, 2
        idx_prj  = 1
        idx_prj2 = 2*idx
        call expand_op_product2(form_pnt,idx_shape,
     &       1d0,9,3,
     &       (/idx_shape,idx_opsin(3),idx_shape,idx_shape,idx_opsin(4),
     &       idx_shape,idx_shape,idx_opsin(4),idx_shape/),
     &       (/        1,           2,        1,        1,           3,
     &       1,        1,           3,        1/),
     &       -1,-1,
     &       (/2,6/),1,
     &       0,0,
     &       (/2,4,1,idx_prj,2,5,1,idx_prj2,3,5,1,idx_prj/),3,
     &       op_info)

        form_pnt => flist
        do while(associated(form_pnt%next))
          form_pnt => form_pnt%next
        enddo

        call expand_op_product2(form_pnt,idx_shape,
     &       1d0,9,3,
     &       (/idx_shape,idx_opsin(3),idx_shape,idx_shape,idx_opsin(4),
     &       idx_shape,idx_shape,idx_opsin(4),idx_shape/),
     &       (/        1,           2,        1,        1,           3,
     &       1,        1,           3,        1/),
     &       -1,-1,
     &       (/2,6/),1,
     &       0,0,
     &       (/2,5,1,idx_prj2,3,5,1,idx_prj/),2,
     &       op_info)

        form_pnt => flist
        do while(associated(form_pnt%next))
          form_pnt => form_pnt%next
        enddo

      enddo

      ! Add the F_{kl}^{p"q"}.G_{p"m}^{r"s}.R_{r"q"}^{ij}.
 100  idx_prj = 4
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj,5,8,1,idx_prj/),3,
     &     op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"q"}.G_{p"q}^{r"s}.R_{r"q"}^{ij}.
      idx_prj = 4
      idx_prj2= 2
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj2,5,8,1,idx_prj/),3,
     &     op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{aq"}.G_{am}^{pb}.R_{bq"}^{ij}.
      idx_prj = 4
      idx_prj2= 2
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj2,2,8,1,idx_prj,5,8,1,idx_prj2/),3,
     &     op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"q"}.G_{p"m}^{pa}.R_{aq"}^{ij}.
      idx_prj = 4
      idx_prj2= 2
      call expand_op_product2(form_pnt,idx_shape,
     &     -2d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj,5,8,1,idx_prj2/),3,
     &     op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"q}.G_{p"m}^{pr}.R_{rq}^{ij}.
      idx_prj  = 4
      idx_prj2 = 2
      idx_prj3 = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     -2d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj2,5,8,1,idx_prj2/),3,
     &     op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      call expand_op_product2(form_pnt,idx_shape,
     &     -2d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj2,5,8,1,idx_prj3/),3,
     &     op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      call expand_op_product2(form_pnt,idx_shape,
     &     -2d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj3,5,8,1,idx_prj2/),3,
     &     op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      call expand_op_product2(form_pnt,idx_shape,
     &     -2d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj3,5,8,1,idx_prj3/),3,
     &     op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"q"}.G_{p"m}^{pn}.R_{nq"}^{ij}.
      idx_prj = 4
      idx_prj2= 1
      call expand_op_product2(form_pnt,idx_shape,
     &     -2d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj,5,8,1,idx_prj2/),3,
     &     op_info)

      ! Point to the formula and move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Add the F_{ij}^{p"q"}.G_{p"m}^{pn}.R_{nq"}^{ij}.
      idx_prj  = 2
      idx_prj2 = 1
      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj,5,8,1,idx_prj/),3,
     &     op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj,5,8,1,idx_prj2/),3,
     &     op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj2,5,8,1,idx_prj/),3,
     &     op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj2,2,8,1,idx_prj,5,8,1,idx_prj/),3,
     &     op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj,2,8,1,idx_prj2,5,8,1,idx_prj2/),3,
     &     op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj2,2,8,1,idx_prj,5,8,1,idx_prj2/),3,
     &     op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj2,2,8,1,idx_prj2,5,8,1,idx_prj/),3,
     &     op_info)

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      call expand_op_product2(form_pnt,idx_shape,
     &     -1d0,9,4,
     &     (/idx_shape,-idx_opsin(2),idx_shape,idx_shape,idx_opsin(3),
     &       idx_shape,idx_shape,idx_opsin(2),idx_shape/),
     &     (/        1,            2,        1,        1,           3,
     &               1,        1,           4,        1/),
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,5,1,idx_prj2,2,8,1,idx_prj2,5,8,1,idx_prj2/),3,
     &     op_info)

 200  if(ntest.ge.100)then
        write(luout,*)'Final formula: Z-Int.'
        call print_form_list(luout,flist,op_info)
      endif

c dbg
c      stop
c dbg

      return
      end
