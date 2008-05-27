      subroutine set_p3f_contract(flist,
     &     idx_opsin,nopsin,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     Routine used to contract the V-intermediate with the F12 integral
*     in order to form P3G.
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
     &     nopsin,idx_opsin(nopsin)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info


      integer ::
     &     idx, idx_shape, nops ,nvtx, idx_temp, ndef
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
        write(luout,*)'P3F-intermediate contraction'
        write(luout,*)'Constituent operators: '
        do idx = 1, nopsin
          write(luout,*)trim(op_info%op_arr(idx_opsin(idx))%op%name)
        enddo
      endif

      ! Get indices of input operators.
      idx_shape = idx_opsin(1)

      ! Point to the formula.
      form_pnt => flist
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Temporary operator.
      call add_operator(op_temp,op_info)
      idx_temp = idx_oplist2(op_temp,op_info)
      op => op_info%op_arr(idx_temp)%op
      ndef = 2
      allocate(occ_def(ngastp,2,ndef))
      occ_def = 0
      occ_def(IHOLE,1,1) = 1
      occ_def(IPART,2,1) = 1
      occ_def(IHOLE,1,2) = 1
      occ_def(IEXTR,2,2) = 1

      call set_uop2(op,op_temp,
     &     occ_def,ndef,1,orb_info)
      deallocate(occ_def)

      ! Arrays for input to expand_op_product2.
      allocate(idxarr(9),svtxarr(9))
      idxarr(1:9) = idx_shape
      idxarr(2) = -idx_opsin(2)
      idxarr(5) = idx_temp
      idxarr(8) = idx_opsin(2)

      svtxarr(1:9) = 1
      svtxarr(2) = 2
      svtxarr(5) = 3
      svtxarr(8) = 4
      
      call expand_op_product2(form_pnt,idx_shape,
     &     1d0,9,4,
     &     idxarr,
     &     svtxarr,
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     (/2,8,1,4/),1,
     &     op_info)

      ! Replace temporary operator with actual F12 integral.
      call form_op_replace2(op_temp,
     &     trim(op_info%op_arr(idx_opsin(3))%op%name),.false.,
     &     flist,op_info)

      if(ntest.ge.100)then
        write(luout,*)'Final formula (P3F intermediate)'
        call print_form_list(luout,flist,op_info)
      endif

      call del_operator(op_temp,op_info)
      deallocate(idxarr,svtxarr)

      return
      end
