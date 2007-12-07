      subroutine inv_op(idx_op_inp,idx_op_inv,
     &     op_info,orb_info,str_info)
*----------------------------------------------------------------------*
*     Wrapper subroutine used in the inversion of the matrix 
*     representation of an operator, op_inp.
*     The resultant is op_inv.
*     GWR November 2007.
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
     &     idx_op_inp, idx_op_inv
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info

      type(operator), pointer ::
     &     op_inp, op_inv
      type(filinf), pointer ::
     &     ffinp, ffinv

      integer ::
     &     njoined, join_off, idx, nocc_cls, iocc_cls
      integer ::
     &     opinp_temp(ngastp,2), opinv_temp(ngastp,2)
      logical ::
     &     closeit
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      logical, external ::
     &     iocc_equal


      call atim_csw(cpu0,sys0,wall0)

      ! Point to the relevant operators and their associated files.
      op_inp => op_info%op_arr(idx_op_inp)%op
      ffinp => op_info%opfil_arr(idx_op_inp)%fhand
      if(.not.associated(ffinp))
     &     call quit(1,'inv_op','no file handle for '//
     &     trim(op_inp%name))

      op_inv => op_info%op_arr(idx_op_inv)%op
      call assign_file_to_op(idx_op_inv,.true.,ffinv,
     &                       1,1,1,
     &                       0,op_info)
      ffinv => op_info%opfil_arr(idx_op_inv)%fhand

      if(ffinv%unit.le.0)then
        call file_open(ffinv)
      endif
      closeit = .false.
      if(ffinp%unit.le.0)then
        call file_open(ffinp)
        closeit = .true.
      endif

      if(ntest.ge.100)then
        write(luout,*) '===================='
        write(luout,*) ' Operator inversion '
        write(luout,*) '===================='
        write(luout,*) 'To be inverted: ',trim(op_inp%name)
        write(luout,*) 'The inverse: ',trim(op_inv%name)
      endif

      ! Check that the two operators have the same shape.
      njoined = op_inp%njoined
      if(njoined.ne.op_inv%njoined)
     &     call quit(1,'inv_op','in and out incompatible: njoined')
      nocc_cls = op_inp%n_occ_cls
      if(nocc_cls.ne.op_inv%n_occ_cls)
     &     call quit(1,'inv_op','in and out incompatible: nocc_cls')
      

      ! Compare the block structures of the operators. The assumption is 
      ! that the comparable blocks are in the same order. 
      opinp_temp(1:ngastp,1:2)=0
      opinv_temp(1:ngastp,1:2)=0

      do iocc_cls = 1, nocc_cls
        join_off = (iocc_cls-1)*njoined

        do idx=1,njoined
          opinp_temp(1:ngastp,1:2) =
     &         op_inp%ihpvca_occ(1:ngastp,1:2,join_off+idx)
          opinv_temp(1:ngastp,1:2) = 
     &         op_inv%ihpvca_occ(1:ngastp,1:2,join_off+idx)

          if (.not.iocc_equal(opinp_temp,op_inp%dagger,
     &         opinv_temp,op_inv%dagger)) then
            call quit(1,'inv_op','in and out incompatible: occs.')
          endif  

        enddo
      enddo

      ! Call the actual inversion routine.
      call invert(ffinp,op_inp,ffinv,op_inv,nocc_cls,
     &     op_info,orb_info)

      if(ntest.ge.1000)then
        call wrt_op_file(luout,5,ffinv,op_inv,1,
     &       op_inv%n_occ_cls,str_info,orb_info)
      endif

      call file_close_keep(ffinv)
      if(closeit) then
        call file_close_keep(ffinp)
      endif

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'time for inversion',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
