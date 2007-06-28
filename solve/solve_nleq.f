*----------------------------------------------------------------------*
      subroutine solve_nleq(
     &     nop_opt,nop_out,idxop_out,idxfil_out,
     &     nop_in,idxop_in,idxfil_in,
     &     ffform_opt,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     solve non-linear equations
*
*     the formula file ffform_opt describes how to calculate 
*     energy and residual
*
*     nop_opt               number of operators to be optimized
*     nop_out               total number of modified operators returned
*     idxop_out(1..nop_opt) indices of those operators in ops(nops)
*     idxop_out(nop_opt+1,..)   residuals
*     idxop_out(2*nop_opt+1,...) other operators updated/modified
*                           energy, intermediates for reuse
*     idxfil_out(2,nop_out) indices and records
*     
*     nop_in                number of operators that need be set on input
*     idxop_in(nop_in)      indices of those operators in ops(nops)
*                           the first nop_out entries define the
*                           preconditioners
*     idxfil_in(2,nop_in)   indices and records
*
*     op_info:  operator definitions and files
*     str_info: string information (to be passed to subroutines)
*     orb_info: orbital space information (to be passed)
*
*----------------------------------------------------------------------*
      implicit none             ! what else ?

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'def_optimize_info.h'
      include 'def_optimize_status.h'
      include 'ifc_memman.h'

      integer, intent(in) ::
     &     nop_in, nop_out, nop_opt,
     &     idxop_in(nop_in), idxfil_in(2,nop_in),
     &     idxop_out(nop_out), idxfil_out(2,nop_out)
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

      logical ::
     &     conv
      integer ::
     &     imacit, imicit, imicit_tot, iprint, task, ifree, iop, nintm
      real(8) ::
     &     energy, xresnrm, xret(10)
      type(operator_array), pointer ::
     &     op_opt(:)
      type(file_array), pointer ::
     &     ffopt(:), ffgrd(:), ffdia(:),
     &     ff_trv(:), ff_h_trv(:)   ! not yet needed
      type(optimize_info) ::
     &     opti_info
      type(optimize_status) ::
     &     opti_stat

      ifree = mem_setmark('solve_nleq')

      ! number of permanent intermediates generated:
      nintm = nop_out - nop_opt*2

      allocate(ffopt(nop_opt),ffdia(nop_opt),
     &     ffgrd(nop_opt),op_opt(nop_opt))
      do iop = 1, nop_opt
        ffopt(iop)%fhand => op_info%opfil_arr(idxfil_out(1,iop))%fhand
        ffgrd(iop)%fhand =>
     &              op_info%opfil_arr(idxfil_out(1,iop+nop_opt))%fhand
        op_opt(iop)%op   => op_info%op_arr(idxop_out(iop))%op
        ffdia(iop)%fhand => op_info%opfil_arr(idxfil_in(1,iop))%fhand
      end do

      ! for savety:
      allocate(ff_trv(1),ff_h_trv(1))

      do iop = 1, nop_opt
        ! open result vector file(s)
        call file_open(ffopt(iop)%fhand)
        ! open corresponding residuals ...
        call file_open(ffgrd(iop)%fhand)
        ! ... and corresponding preconditioner(s)
        if (ffdia(iop)%fhand%unit.le.0)
     &       call file_open(ffdia(iop)%fhand)
      end do

      ! further input operators
      do iop = nop_opt+1,nop_in
        call file_open(op_info%opfil_arr(idxfil_in(1,iop))%fhand)
      end do
      
      ! get initial amplitudes
      do iop = 1, nop_opt
        call zeroop(ffopt(iop)%fhand,op_opt(iop)%op)
      end do

      call set_opti_info(opti_info,nop_opt,op_opt)

      ! start optimization loop
      imacit = 0
      imicit = 0
      imicit_tot = 0
      task = 0
      opt_loop: do while(task.lt.8)

        call optcont
     &       (imacit,imicit,imicit_tot,
     &       task,conv,
     &       energy,xresnrm,
     &       ffopt,ffgrd,ffdia,
     &       ff_trv,ff_h_trv,
     &       opti_info,opti_stat)
        ! 1 - get energy
        ! 2 - get residual
        if (iand(task,1).eq.1.or.iand(task,2).eq.2) then
          call frm_sched(xret,ffform_opt,
     &         op_info,str_info,strmap_info,orb_info)
          ! intermediates should be generated first, so
          energy =  xret(nintm+1)
          xresnrm = xret(nintm+2)
        end if

        write(luout,'(">>>",i3,f24.12,x,g10.4)')imacit,energy,xresnrm

      end do opt_loop

      ifree = mem_flushmark()
      deallocate(ffopt,ffdia,ffgrd,op_opt,ff_trv,ff_h_trv)
      return
      end

