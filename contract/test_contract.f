      subroutine test_contract(form,title,
     &     shape,opsin,nopsin,
     &     op_info)

      implicit none

      integer, parameter ::
     &     ntest = 00

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
        write(lulog,*)'Test Contraction'
        write(lulog,*)'Shape: ',trim(shape)
      endif

      ! Get indices of input operators.
      idx_shape = idx_oplist2(shape,op_info)

      do idx = 1, nopsin
        idx_opsin(idx) = idx_oplist2(opsin(idx),op_info)
      enddo


      ! Initialise the formula.
c      call init_formula(flist)
c      call new_formula_item(flist,command_set_target_init,
c     &     idx_shape)
      form_pnt => flist

      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

      ! Arrays for input to expand_op_product2.
c      nops = nopsin + 1
c      nvtx = nopsin + 2
c      allocate(idxarr(nvtx),svtxarr(nvtx))

c      idxarr(1) = idx_shape
c      idxarr(nvtx) = idx_shape

c      svtxarr(1) = 1
c      svtxarr(nvtx) = 1

c      do idx = 1, nopsin
c        idxarr(1+idx) = idx_opsin(idx)
c        svtxarr(1+idx) = 1+idx
c      enddo 

c      ! Perform thec contraction.
c      call expand_op_product2(form_pnt,idx_shape,
c     &     1d0,nvtx,nops,
c     &     idxarr,
c     &     svtxarr,
c     &     -1,-1,
c     &     0,0,
c     &     0,0,
c     &     0,0,
c     &     .false.,op_info)
 
      allocate(idxarr(9),svtxarr(9))
      idxarr(1:9) = idx_shape
      idxarr(2) = -idx_opsin(1)
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
     &     .false.,op_info)

      ! Write to file.
c      form%comment = trim(title)
c      write(name,'(a,".fml")') trim(form%label)
c      call file_init(form%fhand,name,ftyp_sq_unf,0)
c      call write_form_list(form%fhand,flist,form%comment)

      if(ntest.ge.100)then
        write(lulog,*)'Final formula (test contraction)'
        call print_form_list(lulog,flist,op_info)
      endif

      deallocate(idxarr,svtxarr)
      call dealloc_formula_list(flist)

      return
      end
