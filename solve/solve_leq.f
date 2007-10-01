*----------------------------------------------------------------------*
      subroutine solve_leq(
     &     nop_opt,nroots,nop_out,idxop_out,idxfil_out,
     &     nop_in,idxop_in,idxfil_in,
     &     ffform_opt,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     solve linear equations  Mx = -g
*
*     the formula file ffform_opt describes how to calculate 
*     the matrix trial-vector products and the r.h.s.
*
*     nop_opt               number of x operators to be solved for
*                           in case of coupled equations
*     nroots                number of roots per x operator
*     nop_out               total number of updated operators returned
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
      implicit none             ! for sure

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'def_optimize_info.h'
      include 'def_optimize_status.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     nop_in, nop_out, nop_opt, nroots,
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
     &     iter, iprint, task, ifree, iop, nintm, irequest,
     &     nrequest, nvectors, iroot
      real(8) ::
     &     energy, xresnrm, xret(10)
      type(operator_array), pointer ::
     &     op_opt(:)
      type(file_array), pointer ::
     &     ffdia(:), ff_rhs(:), ff_trv(:),
     &     ffopt(:), ff_mvp(:)
      type(optimize_info) ::
     &     opti_info
      type(optimize_status) ::
     &     opti_stat
      integer, pointer ::
     &     irecmvp(:), irectrv(:)

      character ::
     &     fname*256

      ifree = mem_setmark('solve_nleq')

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'entered solve_leq')
        write(luout,*) 'nop_opt = ',nop_opt
        write(luout,*) 'nroots  = ',nroots
      end if

      if (nop_opt.gt.1)
     &     call quit(1,'solve_leq','did not yet consider coupled LEQs')

      allocate(op_opt(nop_opt))
      do iop = 1, nop_opt
        ! pointer array for operators:
        op_opt(iop)%op   => op_info%op_arr(idxop_out(iop))%op
      end do

      call set_opti_info(opti_info,2,nop_opt,nroots,op_opt)

      ! number of permanent intermediates generated:
      nintm = nop_out - nop_opt*3

      allocate(ff_trv(nop_opt),ff_mvp(nop_opt),ffdia(nop_opt),
     &     ff_rhs(nop_opt))

      nvectors = opti_info%maxsbsp

      do iop = 1, nop_opt
        ! assing trial-vectors to result operator
        allocate(ff_trv(iop)%fhand)
        write(fname,'("trv_",i3.3,".da")') iop
        call file_init(ff_trv(iop)%fhand,fname,ftyp_da_unf,lblk_da)
        call assign_file_to_op(idxop_out(iop),.false.,ff_trv(iop)%fhand,
     &                         1,1,nvectors,
     &                         0,op_info)
        ! assign matrix-vector products to transformation
        allocate(ff_mvp(iop)%fhand)
        write(fname,'("mvp_",i3.3,".da")') iop
        call file_init(ff_mvp(iop)%fhand,fname,ftyp_da_unf,lblk_da)
        call assign_file_to_op(idxop_out(nop_opt+iop),
     &                                        .false.,ff_mvp(iop)%fhand,
     &                         1,1,nvectors,
     &                         0,op_info)

        ! standard name for r.h.s.
        call assign_file_to_op(idxop_out(iop+2*nop_opt),.true.,
     &                                                  ff_rhs(iop),
     &                         1,1,nroots,
     &                         0,op_info)
        ff_rhs(iop)%fhand =>
     &              op_info%opfil_arr(idxop_out(iop+2*nop_opt))%fhand
        ! diagonal exists already
        ffdia(iop)%fhand => op_info%opfil_arr(idxop_in(iop))%fhand

      end do

      ! records with trial vectors and Mv-products, requested by leq_control:
      ifree = mem_alloc_int(irectrv,nroots,'rectrv')
      ifree = mem_alloc_int(irecmvp,nroots,'recmvp')

      do iop = 1, nop_opt
        ! open result vector file(s)
        call file_open(ff_trv(iop)%fhand)
        ! open corresponding matrix vector products ...
        call file_open(ff_mvp(iop)%fhand)
        ! right hand sides ...
        call file_open(ff_rhs(iop)%fhand)
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
        do iroot = 1, nroots
          
          call switch_opfile_record(idxop_out(iop),iroot,op_info)
          call zeroop(ff_trv(iop)%fhand,op_opt(iop)%op)

        end do
      end do

      ! start optimization loop
      iter = 0
      task = 0
      opt_loop: do while(task.lt.8)

        call leq_control
     &       (iter,
     &       task,conv,xresnrm,
     &       nrequest,irectrv,irecmvp,
     &       ffopt,ff_trv,ff_mvp,ff_rhs,ffdia,
     &       opti_info,opti_stat)

        if (iter.gt.1)
     &       write(luout,'(">>>",i3,24x,x,g10.4)')iter-1,xresnrm

        ! 4 - get residual
        if (iand(task,4).eq.4) then
          ! preliminary solution: 
          !   outside loop over requested Mv-products
          do irequest = 1, nrequest
            do iop = 1, nop_opt
              call switch_opfile_record(idxop_out(iop        ),
     &             irectrv(irequest),op_info)
              call switch_opfile_record(idxop_out(iop+nop_opt),
     &             irecmvp(irequest),op_info)
            end do
c dbg
c            call write_title(luout,wst_dbg_subr,'INPUT VECTOR:')
c            call wrt_op_file(luout,4,ff_trv(1)%fhand,
c     &             op_info%op_arr(idxop_out(1))%op,
c     &          1,op_info%op_arr(idxop_out(1))%op%n_occ_cls,
c     &             str_info,orb_info)
c
c dbg
            call frm_sched(xret,ffform_opt,
     &           op_info,str_info,strmap_info,orb_info)
c dbg
c            if (iter.eq.1) then
c              call write_title(luout,wst_dbg_subr,'RHS:')
c              call wrt_op_file(luout,4,ff_rhs(1)%fhand,
c     &             op_info%op_arr(idxop_out(1+2*nop_opt))%op,
c     &          1,op_info%op_arr(idxop_out(1+2*nop_opt))%op%n_occ_cls,
c     &             str_info,orb_info)
c            else
c              call write_title(luout,wst_dbg_subr,'MVP:')
c              call wrt_op_file(luout,4,ff_mvp(1)%fhand,
c     &             op_info%op_arr(idxop_out(1+nop_opt))%op,
c     &          1,op_info%op_arr(idxop_out(1+nop_opt))%op%n_occ_cls,
c     &             str_info,orb_info)
c            end if
c dbg
          end do
        end if

      end do opt_loop

      do iop = 1, nop_opt
        call file_close_delete(ff_trv(iop))
        call file_close_delete(ff_mvp(iop))
        deallocate(ff_trv(iop)%fhand)
        deallocate(ff_mvp(iop)%fhand)
      end do
      deallocate(ff_trv,ff_rhs,ffdia,ff_mvp,op_opt)

      ifree = mem_flushmark()

      return
      end

