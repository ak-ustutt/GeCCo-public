      subroutine set_pint_contract(flist,
     &     idx_opsin,nopsin,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     Routine used to form the various contractions necessary to form
*     the CABS approximation to the P-intermediate.
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

c      type(formula_item), target ::
c     &     flist

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
      call set_1contrib(form_pnt,1d0,7,
     &     idx_shape,idx_opsin,nopsin,op_info)

      ! Move to the end of the list.
      do while(associated(form_pnt%next))
        form_pnt => form_pnt%next
      enddo

c dbg
      goto 100
c dbg

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
c      call form_op_replace2(op_temp,opsin(1),.true.,
c     &     flist,op_info)

c      ! Write to file.
c      form%comment = trim(title)
c      write(name,'(a,".fml")') trim(form%label)
c      call file_init(form%fhand,name,ftyp_sq_unf,0)
c      call write_form_list(form%fhand,flist,form%comment)

 100  if(ntest.ge.100)then
        write(luout,*)'Final formula: P-Int.)'
        call print_form_list(luout,flist,op_info)
      endif

c dbg
      stop
c dbg

      call del_operator(op_temp,op_info)
      deallocate(idxarr,svtxarr)
c      call dealloc_formula_list(flist)

      return
      end
