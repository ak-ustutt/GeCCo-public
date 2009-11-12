      subroutine set_p3g_contract(form,title,
     &     shape,opsin,nopsin,
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

      type(formula), intent(inout) ::
     &     form
      integer, intent(in) ::
     &     nopsin
      character(*), intent(in) ::
     &     title,shape,opsin(nopsin)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info


      integer ::
     &     idx, idx_shape, nops ,nvtx, idx_temp, ndef
      character ::
     &     name*(form_maxlen_label*2)

      integer ::
     &     idx_opsin(nopsin)
      integer, pointer ::
     &     idxarr(:), svtxarr(:), occ_def(:,:,:)
      character(4), parameter ::
     &     op_temp = 'TEMP'

      type(formula_item), target ::
     &     flist

      type(operator), pointer ::
     &     op
      type(formula_item), pointer ::
     &     form_pnt

      integer, external ::
     &     idx_oplist2


      call quit(1,'set_p3g_contract','obsolete?')

      if(ntest.ge.100)then
        write(luout,*)'Test Contraction'
        write(luout,*)'Shape: ',trim(shape)
      endif

      ! Get indices of input operators.
      idx_shape = idx_oplist2(shape,op_info)

      do idx = 1, nopsin
        idx_opsin(idx) = idx_oplist2(opsin(idx),op_info)
      enddo


      ! Initialise the formula.
      call init_formula(flist)
      call new_formula_item(flist,command_set_target_init,
     &     idx_shape)
      form_pnt => flist

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Temporary operator.
      call add_operator(op_temp,op_info)
      idx_temp = idx_oplist2(op_temp,op_info)
      op => op_info%op_arr(idx_temp)%op
      allocate(occ_def(ngastp,2,1))
      ndef = 1
      occ_def = 0
      occ_def(IHOLE,1,1) = 2
      occ_def(IPART,2,1) = 1
      occ_def(IEXTR,2,1) = 1

      call set_uop2(op,op_temp,
     &     occ_def,ndef,1,orb_info)
      deallocate(occ_def)

      ! Arrays for input to expand_op_product2.
      allocate(idxarr(9),svtxarr(9))
      idxarr(1:9) = idx_shape
      idxarr(2) = idx_temp
c      idxarr(2) = -idx_opsin(1)
      idxarr(5) = idx_opsin(2)
      idxarr(8) = idx_opsin(2)

      svtxarr(1:9) = 1
      svtxarr(2) = 2
      svtxarr(5) = 3
      svtxarr(8) = 3
      
      call expand_op_product2(form_pnt,idx_shape,
     &     1d0,9,3,
     &     idxarr,
     &     svtxarr,
     &     -1,-1,
     &     (/2,5/),1,
     &     (/2,8/),1,
     &     0,0,
     &     op_info)

      ! Replace temporary operator with actual F12 integral.
      call form_op_replace2(op_temp,opsin(1),.true.,
     &     flist,op_info)

      ! Write to file.
      form%comment = trim(title)
      write(name,'(a,".fml")') trim(form%label)
      call file_init(form%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form%fhand,flist,form%comment)

      if(ntest.ge.100)then
        write(luout,*)'Final formula (test contraction)'
        call print_form_list(luout,flist,op_info)
      endif

      call del_operator(op_temp,op_info)
      deallocate(idxarr,svtxarr)
      call dealloc_formula_list(flist)

      return
      end
