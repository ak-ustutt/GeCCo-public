      subroutine mult_op(nops_in,idxop_in,idxop_out,
     &     op_info,orb_info,str_info)
*----------------------------------------------------------------------*
*     Wrapper routine to direct the multiplication of the matrix 
*     representation of nop_in operators.
*     The array, idx_op_in, contains the indices of those operators. 
*     The multiplication proceeds along the array:
*     op1*op2*op3*.... = op_out
*     
*     GWR November 2007
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'
      include 'def_operator.h'
      include 'def_operator_list.h'
      include 'def_operator_array.h'
      include 'def_file_list.h'
      include 'def_file_array.h'
      include 'def_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      integer, intent(in) ::
     &     nops_in, idxop_in(nops_in), idxop_out
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info

      type(file_array) ::
     &     ffin(nops_in)
      type(operator_array) ::
     &     ops_in(nops_in)
      type(filinf), pointer ::
     &     ffout
      type(operator), pointer ::
     &     op_out

      integer ::
     &     iop, iocc_cls, nocc_cls, idx, njoined, join_off
      logical ::
     &     closeit(nops_in), closeout
      integer ::
     &     opin_temp(ngastp,2), opout_temp(ngastp,2)

      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      integer, external ::
     &     idx_oplist2
      logical, external ::
     &     iocc_equal

      call atim_csw(cpu0,sys0,wall0)


      if(ntest.ge.100)then
        write(luout,*)'================================'
        write(luout,*)' Operator matrix multiplication '
        write(luout,*)'================================'
      endif

      closeit(1:nops_in) = .false.
      do iop = 1, nops_in
        ffin(iop)%fhand => op_info%opfil_arr(idxop_in(iop))%fhand
        if(ffin(iop)%fhand%unit.le.0)then
          call file_open(ffin(iop)%fhand)
          closeit(iop) = .true.
        endif
        ops_in(iop)%op => op_info%op_arr(idxop_in(iop))%op
      enddo

      op_out => op_info%op_arr(idxop_out)%op
      if(.not.associated(op_info%opfil_arr(idxop_out)%fhand))then
        call assign_file_to_op(idxop_out,.true.,ffout,
     &                         1,1,1,
     &                         0,op_info)
      endif

      ffout => op_info%opfil_arr(idxop_out)%fhand
      closeout = .false.
      if(ffout%unit.le.0)then
        call file_open(ffout)
        closeout = .true.
      endif
        
      ! Check to see whether the input and output operators have the 
      ! same shapes in each block. Assumes that the equivalent blocks
      ! are ordered the same in all operators.
      opin_temp(1:ngastp,1:2)=0
      opout_temp(1:ngastp,1:2)=0
      nocc_cls = op_out%n_occ_cls
      njoined = op_out%njoined

      do iop = 1, nops_in
        if(nocc_cls.ne.ops_in(iop)%op%n_occ_cls)
     &       call quit(1,'mult_op','Unequal no. of occupation classes')
        if(njoined.ne.ops_in(iop)%op%njoined)
     &       call quit(1,'mult_op','Unequal no. of supervertices')

        do iocc_cls = 1, nocc_cls
          join_off = (iocc_cls-1)*njoined

          do idx=1,njoined

            opin_temp(1:ngastp,1:2) =
     &           ops_in(iop)%op%ihpvca_occ(1:ngastp,1:2,join_off+idx)
            opout_temp(1:ngastp,1:2) = 
     &           op_out%ihpvca_occ(1:ngastp,1:2,join_off+idx)
            
            if (.not.iocc_equal(opin_temp,ops_in(iop)%op%dagger,
     &           opout_temp,op_out%dagger)) then
              write(luout,*)'Operator in: ',trim(ops_in(iop)%op%name)
              write(luout,*)'Operator out: ',trim(op_out%name)
              call quit(1,'mult_op','in and out incompatible: occs.')
            endif  
        
          enddo
        enddo
      enddo

      ! Call the multiplication routine.
c      call op_mult(-1d0,ffin(1),ops_in(1)%op,1,
c     &     ffout,op_out,1,orb_info)
      call op_mult2(-1d0,ffin,ops_in,nops_in,ffout,op_out,
     &     op_info,orb_info)

      if(ntest.ge.1000)then
        write(luout,*) 'Product operator: ',trim(op_out%name)
        call wrt_op_file(luout,5,ffout,op_out,1,
     &       nocc_cls,str_info,orb_info)
      endif

c      stop

      do iop = 1, nops_in
        if(closeit(iop))
     &       call file_close_keep(ffin(iop)%fhand)
      enddo
      ! Close the output if required, bearing in mind that it could be 
      ! the same as an input file.
      if(closeout.and.ffout%unit.gt.0)
     &     call file_close_keep(ffout)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'time for multiplication',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
