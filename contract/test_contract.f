      subroutine test_contract(form,title,
     &     shape,opsin,nopsin,
     &     op_info)

      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      type(formula), intent(inout) ::
     &     form
      integer, intent(in) ::
     &     nopsin
      character(*), intent(in) ::
     &     title,shape,opsin(nopsin)
      type(operator_info), intent(inout) ::
     &     op_info


      integer ::
     &     idx, idx_shape, nops ,nvtx
      character ::
     &     name*(form_maxlen_label*2)

      integer ::
     &     idx_opsin(nopsin)
      integer, pointer ::
     &     idxarr(:), svtxarr(:)

      type(formula_item), target ::
     &     flist

      type(formula_item), pointer ::
     &     form_pnt

      integer, external ::
     &     idx_oplist2


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

      ! Arrays for input to expand_op_product2.
      nops = nopsin + 1
      nvtx = nopsin + 2
      allocate(idxarr(nvtx),svtxarr(nvtx))

      idxarr(1) = idx_shape
      idxarr(nvtx) = idx_shape

      svtxarr(1) = 1
      svtxarr(nvtx) = 1

      do idx = 1, nopsin
        idxarr(1+idx) = idx_opsin(idx)
        svtxarr(1+idx) = 1+idx
      enddo 

      ! Perform thec contraction.
      call expand_op_product2(form_pnt,idx_shape,
     &     1d0,nvtx,nops,
     &     idxarr,
     &     svtxarr,
     &     -1,-1,
     &     0,0,
     &     0,0,
     &     0,0,
     &     op_info)
      

      ! Write to file.
      form%comment = trim(title)
      write(name,'(a,".fml")') trim(form%label)
      call file_init(form%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form%fhand,flist,form%comment)

      if(ntest.ge.100)then
        write(luout,*)'Final formula (test contraction)'
        call print_form_list(luout,flist,op_info)
      endif

      deallocate(idxarr,svtxarr)
      call dealloc_formula_list(flist)

      return
      end
