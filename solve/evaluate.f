*----------------------------------------------------------------------*
      subroutine evaluate(
     &     nop_out,idxop_out,
     &     nop_in,idxop_in,
     &     ffform_opt,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     just evaluate given formula
*
*     the formula file ffform_opt describes how to calculate 
*     the matrix trial-vector products and the r.h.s.
*
*     nop_out               total number of updated operators returned
*     idxop_out(1..nop_opt) indices of operators
*     
*     nop_in                number of operators that need be set on input
*     idxop_in(nop_in)      indices of operators
*
*     op_info:  operator definitions and files
*     str_info: string information (to be passed to subroutines)
*     orb_info: orbital space information (to be passed)
*
*     andreas, Aug. 2007
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nop_in, nop_out,
     &     idxop_in(nop_in),
     &     idxop_out(nop_out)
      type(filinf), intent(inout) ::
     &     ffform_opt
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf) ::
     &     orb_info

      integer ::
     &     ifree, iop
      real(8) ::
     &     energy, xresnrm, xret(10*nop_out)
      type(filinf) ::
     &     ffdum

      ifree = mem_setmark('evaluate')

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'entered evaluate')
      end if

      ! assign files to output operators
      do iop = 1, nop_out

        ! standard names (ffdum is a dummy)
        call assign_file_to_op(idxop_out(iop),.true.,ffdum,
     &                         1,1,1,
     &                         0,op_info)

      end do

      ! call the scheduler
      call frm_sched(xret,ffform_opt,
     &               op_info,str_info,strmap_info,orb_info)

      if (iprlvl.ge.5) then
        call write_title(luout,wst_title,'norms of output operators')
        do iop = 1, nop_out
          write(luout,'(4x,i4," - ",g12.6)') iop, xret(iop)
        end do
      end if

      if (ntest.ge.1000) then
        do iop = 1, nop_out
          write(luout,*) 'dump of result for ',
     &         trim(op_info%op_arr(idxop_out(iop))%op%name)
          call wrt_op_file(luout,5,
     &       op_info%opfil_arr(idxop_out(iop))%fhand,
     &       op_info%op_arr(idxop_out(iop))%op,
     &       1,op_info%op_arr(idxop_out(iop))%op%n_occ_cls,
     &       str_info,orb_info)
        end do
      end if

      ifree = mem_flushmark()

      return
      end

