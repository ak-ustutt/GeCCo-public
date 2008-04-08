      subroutine inter_full_contr(form_list,idx_res,
     &     fac,nops,idx_op,op_info)
*-----------------------------------------------------------------------
*     Routine used to form the full contraction between an intermediate 
*     type operator and the necessary coefficients and then write it to
*     the requested formula.
*     form_list: the formula to which the result is written.
*     idx_res: the index of the operator which describes the outer arcs 
*              of the result.
*     fac: numerical factor by which the contraction is pre-multiplied.
*     idx_op(nops): indices of all operators to be contracted.
*     iblk_min_in(nops): first block of each operator to be considered.
*     iblk_max_in(nops): last block of each operator to be considered.
*     op_info: information on all operators.
*
*     GWR Sept. 2007.
*-----------------------------------------------------------------------
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'ifc_operators.h'
      include 'ifc_baserout.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(inout), target::
     &     form_list
      real(8), intent(in)::
     &     fac
      integer, intent(in)::
     &     nops, idx_res, idx_op(nops)
      type(operator_info), intent(in)::
     &     op_info

      ! Local variables.
      type(formula_item), pointer ::
     &     form_pnt
      type(contraction) ::
     &     proto
      type(operator), pointer ::
     &     op_res, op_int, op_spc

      integer ::
     &     nvtx, narc, idx, ivtx, sidx, n_occ_cls, iocc, iblk_off, ispc,
     &     idx_spc, iarc
      integer ::
     &     occ_temp(ngastp,2)

      integer, external ::
     &     iblk_occ

      if(ntest.ge.100)then
        write(luout,*) '==========================='
        write(luout,*) 'Info from inter_full_contr'
        write(luout,*) '==========================='
        write(luout,*) 'idx_res = ',idx_res
        write(luout,*) 'fac, nops: ',fac, nops
        write(luout,*) 'idx_op: ',idx_op(1:nops)
      endif

      form_pnt => form_list

      op_res => op_info%op_arr(idx_res)%op
      op_int => op_info%op_arr(idx_op(1))%op 
      n_occ_cls=op_int%n_occ_cls

      ! Initialize the contraction.
      call init_contr(proto)
      nvtx=2*nops-1
      narc=2*(nops-1)
      call resize_contr(proto,nvtx,narc,0,0)
      proto%nvtx=nvtx
      proto%narc=narc
      proto%fac=fac
      proto%idx_res=idx_res

      ! Further information on supervertices.
      proto%nsupvtx = nops
      
      ! Intermediate.
      proto%joined(0,1)=nops
      proto%joined(1:nvtx,1)=0
      ivtx=1
      do idx=1,nops
        proto%joined(idx,1)=ivtx
        ivtx=ivtx+2
      enddo  
      proto%svertex(1:nvtx:2)=1
      
      ! Spacers.
      sidx=2
      do idx=2,nops
        proto%joined(0,idx)=1
        proto%joined(1,idx)=sidx
        proto%joined(2:nvtx,idx)=0

        proto%svertex(sidx)=idx
        sidx=sidx+2
      enddo

      ! Loop over the blocks of the intermediate.
      do iocc=1,n_occ_cls

        iblk_off=(iocc-1)*nops

        ! Set the first intermediate vertex.
c        do idx=1,nops
          proto%vertex(1)%idx_op=idx_op(1)
          proto%vertex(1)%iblk_op=1
c        enddo

        ! Check that the outer arcs of the intermediate correspond to 
        ! those of the resultant.
        occ_temp=iocc_xdn(1,op_info%op_arr(idx_op(1))%
     &       op%ihpvca_occ(1:ngastp,1:2,iblk_off+1))
     &       +iocc_xdn(2,op_info%op_arr(idx_op(1))%
     &       op%ihpvca_occ(1:ngastp,1:2,iblk_off+nops))

        if(iblk_occ(occ_temp,op_res%dagger,op_res).le.0)
     &       call quit(1,'inter_full_contr','intermed/result mismatch') 

        ! Loop over the vertices of the intermediate in order to determine
        ! the correct blocks of the spacer operators to use in the
        ! contractions.
        ivtx=1
        iarc=0
        do ispc=2,nops

          ! Form the shape of the necessary spacer operator.
          occ_temp=iocc_xdn(2,op_info%op_arr(idx_op(1))%
     &         op%ihpvca_occ(1:ngastp,1:2,iblk_off+ispc-1))
     &         + iocc_xdn(1,op_info%op_arr(idx_op(1))%
     &         op%ihpvca_occ(1:ngastp,1:2,iblk_off+ispc))

          op_spc => op_info%op_arr(idx_op(ispc))%op
          idx_spc=iblk_occ(occ_temp,op_spc%dagger,op_spc)
          if(idx_spc.le.0)
     &         call quit(1,'inter_full_contr','unrecognised spacer blk')

          ! Add the spacer vertex to the contraction.
          ivtx=ivtx+1
          proto%vertex(ivtx)%idx_op=idx_op(ispc)
          proto%vertex(ivtx)%iblk_op=idx_spc
          
          ! Add the next intermediate vertex.
          ivtx=ivtx+1
          proto%vertex(ivtx)%idx_op=idx_op(1)
          proto%vertex(ivtx)%iblk_op=ispc

          ! Set the actual contractions.
          iarc=iarc+1
          proto%arc(iarc)%link(1)=ivtx-2
          proto%arc(iarc)%link(2)=ivtx-1
          proto%arc(iarc)%occ_cnt=iocc_xdn(2,op_info%op_arr(idx_op(1))%
     &         op%ihpvca_occ(1:ngastp,1:2,iblk_off+ispc-1))

          iarc=iarc+1
          proto%arc(iarc)%link(1)=ivtx-1
          proto%arc(iarc)%link(2)=ivtx
          proto%arc(iarc)%occ_cnt=iocc_dagger(iocc_xdn(1,
     &         op_info%op_arr(idx_op(1))%
     &         op%ihpvca_occ(1:ngastp,1:2,iblk_off+ispc)))

        enddo

c        if(ntest.ge.100)then
c          call prt_contr2(luout,proto,op_info)
c        endif

        ! Save the contraction to the formula.
        form_pnt%command = command_add_contribution
        form_pnt%target = proto%idx_res
        allocate(form_pnt%contr,form_pnt%next)
        call init_contr(form_pnt%contr)
        call copy_contr(proto,form_pnt%contr)
        form_pnt%next%prev => form_pnt
        form_pnt => form_pnt%next
        form_pnt%next => null()
        form_pnt%contr => null()
        form_pnt%interm => null()
        form_pnt%command = command_end_of_formula

      enddo
      
      if(ntest.ge.100)then
        call print_form_list(luout,form_list,op_info)
      endif

      call dealloc_contr(proto)

      return
      end
