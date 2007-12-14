*----------------------------------------------------------------------*
      subroutine solve_nleq(
     &     nop_opt,nop_out,idxop_out,
     &     nop_in,idxop_in,
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
*     
*     nop_in                number of operators that need be set on input
*     idxop_in(nop_in)      indices of those operators in ops(nops)
*                           the first nop_out entries define the
*                           preconditioners
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

      logical ::
     &     conv
      integer ::
     &     imacit, imicit, imicit_tot, iprint, task, ifree, iop, nintm
      real(8) ::
     &     energy, xresnrm(nop_opt), xret(10)
      type(operator_array), pointer ::
     &     op_opt(:), op_grd(:)
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

      do iop = 1, nintm
        call assign_file_to_op(idxop_out(nop_opt*2+iop),
     &                                             .true.,ffopt(iop),
     &                         1,1,1,
     &                         0,op_info)
      end do

      allocate(ffopt(nop_opt),ffdia(nop_opt),
     &     ffgrd(nop_opt),op_opt(nop_opt),op_grd(nop_opt))
      do iop = 1, nop_opt
        call assign_file_to_op(idxop_out(iop),.true.,ffopt(iop),!<-dummy 
     &                         1,1,1,
     &                         0,op_info)
        ffopt(iop)%fhand => op_info%opfil_arr(idxop_out(iop))%fhand
        call assign_file_to_op(idxop_out(iop+nop_opt),.true.,ffgrd(iop),!<-dummy
     &                         1,1,1,
     &                         0,op_info)
        ffgrd(iop)%fhand =>
     &              op_info%opfil_arr(idxop_out(iop+nop_opt))%fhand
        op_grd(iop)%op   => op_info%op_arr(idxop_out(iop+nop_opt))%op
        op_opt(iop)%op   => op_info%op_arr(idxop_out(iop))%op
        ffdia(iop)%fhand => op_info%opfil_arr(idxop_in(iop))%fhand
      end do

      ! for safety reasons, we allocate the two guys
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
      
      ! get initial amplitudes
      do iop = 1, nop_opt
c dbg
c        print *,'zeroing: ',trim(op_opt(iop)%op%name),
c     &       ' file: ',trim(ffopt(iop)%fhand%name)
c dbg
        call zeroop(ffopt(iop)%fhand,op_opt(iop)%op)
      end do

      call set_opti_info(opti_info,1,nop_opt,1,op_opt)

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
     &       opti_info,opti_stat,
     &       op_info,orb_info)
        ! 1 - get energy
        ! 2 - get residual
        if (iand(task,1).eq.1.or.iand(task,2).eq.2) then
          call frm_sched(xret,ffform_opt,
     &         op_info,str_info,strmap_info,orb_info)
          ! intermediates should be generated first, energy
          ! is expected to be the last "intermediate"
          energy =  xret(nintm)
          xresnrm(1:nop_opt) = xret(nintm+1:nintm+nop_opt)
c          if (nopt.eq.2) xres2nrm = xret(nintm+2)
        end if

        if (nop_opt.eq.1)
     &       write(luout,'(">>>",i3,f24.12,x,g10.4)')
     &       imacit,energy,xresnrm(1)
        if (nop_opt.eq.2)
     &       write(luout,'(">>>",i3,f24.12,2(x,g10.4))')
     &       imacit,energy,xresnrm(1:2)

      end do opt_loop

      ifree = mem_flushmark()
      deallocate(ffopt,ffdia,ffgrd,op_opt,op_grd,ff_trv,ff_h_trv)
      return
      end

