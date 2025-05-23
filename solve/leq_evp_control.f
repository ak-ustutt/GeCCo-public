*----------------------------------------------------------------------*
      subroutine leq_evp_control(modestr,iter,
     &                   task,conv,xrsnrm,xeig,
     &                   use_s,
     &                   nrequest,irectrv,irecmvp,irecmet,
     &                   me_opt,me_scr,me_trv,me_mvp,me_met,me_rhs,
     &                   me_dia,me_ext,
     &                   me_special,nspecial,nextra,idxspc,
c     &                   ffopt,fftrv,ffmvp,ffmet,ffrhs,ffdia,
     &                   flist,depend,
     &                   fspc,nspcfrm,
     &                   opti_info,opti_stat,
     &                   orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*
* Control routine for linear equations and eigenvalue problems.
*
* The optimization strategy is set in opti_info, see def_optinf.h 
* for further comments.
*
* It keeps track of the iterations and expects to be
* called with iter set to zero at the first time.
*
* leq_control() requests information via the task parameter (bit flags):
*     
*     1: unused
*     2: unused
*     4: calculate matrix-trialvector product
*     8: exit
*
* Files: 
* a) passed to slave routines
*   fftrv   -- new trialvectors, on convergence:
*              solution vectors to linear equations
*
* b) passed from slave routines
*   ffmvp  -- matrix-trialvector product
*   ffrhs  -- right-hand side (for leq)
*
* c) own book keeping files
*   ff_h_trv_old  -- previous H/A matrix-vector product
*   ff_newstp -- current new step (scratch)
*   ffst_sbsp -- subspace of previous steps
*   ffgv_sbsp -- subspace of previous gradient/vector function diff.s
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'
      include 'def_file_array.h'
      include 'def_optimize_info.h'
      include 'def_optimize_status.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_dependency_info.h'
      include 'molpro_out.h'
      
* parameters
      integer, parameter ::
     &     ntest = 000
      character(len=10), parameter ::
     &     name_alg_leq(0:1) =
     &     (/"  CONJGRAD","  SUBSPACE"/),
     &     name_alg_evp(0:1) =
     &     (/"  xxxxxxxx","  DAVIDSON"/)

      integer, intent(out) ::
     &     task
      logical, intent(out) ::
     &     conv
      character(*), intent(in) ::
     &     modestr
      logical, intent(in) ::
     &     use_s(*)

      integer, intent(inout) ::
     &     iter

      integer, intent(out) ::
     &     nrequest,irectrv(*),irecmvp(*),irecmet(*)

      integer, intent(in) ::
     &     nspecial, nextra, idxspc, nspcfrm

      type(me_list_array), intent(in) ::
     &     me_opt(*), me_dia(*), me_special(nspecial+nextra),
     &     me_mvp(*), me_rhs(*), me_scr(*), me_ext(*)
      type(me_list_array), intent(inout) ::
     &     me_met(*), me_trv(*)
c      type(file_array), intent(in) ::
c     &     ffopt(*), fftrv(*), ffmvp(*), ffmet(*), ffrhs(*), ffdia(*)
      
      type(formula_item), intent(inout) ::
     &     flist
      type(dependency_info) ::
     &     depend
      type(formula_item), intent(in) ::
     &     fspc(nspcfrm)

      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout) ::
     &     opti_stat

      type(orbinf),intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf),intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info

      real(8), intent(inout) ::
     &     xrsnrm(opti_info%nroot*opti_info%nopt),
     &     xeig(opti_info%nroot,2)


*     buffers for incore/out-of-core work:
      integer ::
     &     nincore, nbuf, lenbuf
      real(8), pointer ::
     &     xbuf1(:), xbuf2(:), xbuf3(:)

      character(60) ::
     &     fname
      logical ::
     &     lexit, lconv, triv_s
      integer ::
     &     iprint, iroute, ifree, iopt, idx, iroot, nsub,
     &     irequest
      real(8) ::
     &     de, cpu, sys, wall, cpu0, sys0, wall0


      iprint = max(ntest,iprlvl)
      call atim_csw(cpu0,sys0,wall0)

      if (trim(modestr).ne.'LEQ'.and.trim(modestr).ne.'EVP')
     &     call quit(1,'leq_evp_control','unknown mode - '//
     &     trim(modestr))

* be verbose?
      if (iprint.ge.5) then
        if (modestr(1:3).eq.'LEQ')
     &       call write_title(lulog,wst_section,'LEQ solver')
        if (modestr(1:3).eq.'EVP')
     &       call write_title(lulog,wst_section,'EVP solver')
      end if
      if (ntest.ge.10) then
        write(lulog,*) 'entered leq_evp_control with:'
        write(lulog,*) ' iter: ', iter
      end if

* set iroute:
* 0 conjugate gradient
* 1 subspace expansion
      if (modestr(1:3).eq.'LEQ') iroute = opti_info%mode_leq
* 0 unused
* 1 Davidson        
      if (modestr(1:3).eq.'EVP') iroute = opti_info%mode_evp

      if (iprint.ge.5.or.iter.eq.0) then
        if (modestr(1:3).eq.'LEQ')
     &       write(lulog,*) 'Optimization algorithm:    ',
     &       name_alg_leq(iroute)
        if (modestr(1:3).eq.'EVP')
     &       write(lulog,*) 'Optimization algorithm:    ',
     &       name_alg_evp(iroute)
        write(lulog,'(x,a,i10)')
     &                 'Max. number of iterations: ',opti_info%maxmacit
        write(lulog,'(x,a,e10.2)')
     &                 'Threshold for residual:    ',opti_info%thrgrd(1)
        write(lulog,'(x,a,3i10)')
     &                 'Number of parameters:      ',
     &       opti_info%nwfpar(1:opti_info%nopt)
      end if

      if (ntest.ge.10) then
        write(lulog,*) 'our route is: ',iroute
      end if
*======================================================================*
*     unless this is the initial call:
*     get memory for buffers, decide on in-core/out-of-core
*======================================================================*
      ! set mark for temporarily allocated stuff
      ifree = mem_setmark('leqevpc_temp')
      if (iter.ne.0.or.modestr(1:3).eq.'LEQ') then
        call leqevpc_mem(nincore,lenbuf,
     &                   ifree,opti_info%nwfpar,opti_info%nopt)
c        if (nincore.le.1) then
        if (iter.ne.0) then
          do iopt = 1, opti_info%nopt
            call file_open(opti_stat%ffscr(iopt)%fhand)
            call file_open(opti_stat%ffext(iopt)%fhand)
          end do
        end if
c        end if
      end if

      ! still, we assume this:
      conv = .false.
*======================================================================*
* check iteration number:
* zero --> init everything
*======================================================================*
      if (iter.eq.0) then

        ! set mark for memory allocated until optimization ends
        ifree = mem_setmark('leqevpc_perm')
        call leqevpc_init()

* set task -- mv-product:
        task = 4

        triv_s = .false.  ! needed for quick return for trivial solution (LEQ)

        ! request slaves to put first few trial vectors and
        ! corresponding Mv-product to first nroot records on the
        ! respective files:        
        if (modestr(1:3).eq.'LEQ') then
          call leqc_init(xrsnrm,iroute,
     &       me_opt,me_trv,me_mvp,me_rhs,me_dia,me_met,me_scr,
     &       me_special,nspecial+nextra,
     &       nincore,lenbuf,
     &       xbuf1,xbuf2,xbuf3,
     &       flist,depend,use_s,
     &       fspc,nspcfrm,
     &       opti_info,opti_stat,
     &       orb_info,op_info,str_info,strmap_info)          
          opti_stat%nadd = opti_info%nroot ! ?? <-- check that for LEQ
          triv_s = (opti_info%nroot==1.and.opti_info%nopt==1.and.
     &        xrsnrm(1).lt.opti_info%thrgrd(1))
        else if (modestr(1:3).eq.'EVP') then
          opti_stat%nadd = opti_info%nroot ! <-- check that for LEQ
        else
          call quit(1,'leq_evp_control','this should not happen!')
        end if
        opti_stat%ndel = 0
        nrequest = opti_info%nroot
        do irequest = 1, nrequest
          irectrv(irequest) = irequest
          irecmvp(irequest) = irequest
          irecmet(irequest) = irequest
        end do

        do iopt = 1, opti_info%nopt
          if (opti_stat%ffscr(iopt)%fhand%unit.gt.0) 
     &             call file_close_delete(opti_stat%ffscr(iopt)%fhand)
          if (opti_stat%ffext(iopt)%fhand%unit.gt.0)
     &             call file_close_delete(opti_stat%ffext(iopt)%fhand)
        end do          

        iter = 1

        if (triv_s) then
           write(lulog,'(1x,"DETECTED TRIVIAL SOLUTION")')
           if (luout.ne.lulog)
     &       write(luout,'(1x,"DETECTED TRIVIAL SOLUTION")')
          do iopt = 1, opti_info%nopt
            if (opti_stat%ffscr(iopt)%fhand%unit.gt.0) 
     &             call file_close_delete(opti_stat%ffscr(iopt)%fhand)
            if (opti_stat%ffext(iopt)%fhand%unit.gt.0)
     &             call file_close_delete(opti_stat%ffext(iopt)%fhand)
          end do

          call leqevpc_cleanup()
          ! release all temporary memory
          ifree = mem_flushmark('leqevpc_temp')

          ! flush all other memory allocated by optimizer
          ifree = mem_flushmark('leqevpc_perm')

          task = 8 ! stop it

          return

        end if

        ifree = mem_flushmark('leqevpc_temp')

* return and let the slaves work
        return

*======================================================================*
* here the actual LEQ/EVP solvers come:
*======================================================================*
      else 
        if (ntest.ge.10) then
          write(lulog,*) 'getting new vectors ...'
        end if

        if (modestr(1:3).eq.'LEQ') then
          call leqc_core(iter,
     &         task,iroute,xrsnrm,
     &         use_s,
     &         me_opt,me_trv,me_mvp,me_rhs,me_dia,me_met,me_scr,me_ext,
     &         me_special,nspecial+nextra,
c     &         ffopt,fftrv,ffmvp,ffrhs,ffdia,
     &         nincore,lenbuf,
     &         xbuf1,xbuf2,xbuf3,
     &         flist,depend,
     &         fspc,nspcfrm,
     &         opti_info,opti_stat,
     &         orb_info,op_info,str_info,strmap_info)
        else
          call evpc_core(iter,
     &         task,iroute,xrsnrm,xeig,
     &         use_s,
     &         me_opt,me_trv,me_mvp,me_dia,me_met,me_scr,me_ext,
     &         me_special,nspecial,nextra,idxspc,
c     &         ffopt,fftrv,ffmvp,ffdia,
     &         nincore,lenbuf,
     &         xbuf1,xbuf2,xbuf3,
     &         flist,depend,
     &         fspc,nspcfrm,
     &         opti_info,opti_stat,
     &         orb_info,op_info,str_info,strmap_info)
        end if

      end if

      lexit = .false.
*======================================================================*
* end of macro-iteration:
*  check convergence and max. iterations:
*======================================================================*
      lconv = .true.
      idx = 0
      do iopt = 1, opti_info%nopt
        do iroot = 1, opti_info%nroot
          idx = idx+1
          lconv = lconv.and.xrsnrm(idx).lt.opti_info%thrgrd(iopt)
        end do        
      end do

      if (iprint.ge.5) then
        if (opti_info%norder.eq.1)
     &       write(lulog,*) 'after iteration ',iter
        idx = 0
        do iopt = 1, opti_info%nopt
          do iroot = 1, opti_info%nroot
            idx = idx+1
            write(lulog,'(x,2(a,e10.3),a,l)')
     &                      ' norm of residual:  ', xrsnrm(idx),
     &                           '   threshold:  ',
     &                                   opti_info%thrgrd(iopt),
     &                           '   converged:  ',
     &                   xrsnrm(idx).lt.opti_info%thrgrd(iopt)
          end do
        end do
      end if

      if (lconv)
     &       write(lulog,'(x,a,i5,a)')
     &         'CONVERGED IN ',iter,' ITERATIONS'
      if (lconv.and.luout.ne.lulog.and.iprlvl.ge.5 .and. .not. lmol)
     &       write(luout,'(x,a,i5,a)')
     &         'CONVERGED IN ',iter,' ITERATIONS'
      if (lconv) conv = .true.

      if (.not.lconv) iter = iter + 1

      if (.not.lconv.and.
     &       (iter.gt.opti_info%maxmacit)) then
        call warn('linear solver', 'NO CONVERGENCE OBTAINED')
        !iter = iter - 1
        lexit = .true.
      end if

      if (.not.(lconv.or.lexit)) then
*----------------------------------------------------------------------*
* do some stuff for the next macro-iteration
*----------------------------------------------------------------------*

        do iopt = 1, opti_info%nopt
          if (opti_stat%ffscr(iopt)%fhand%unit.gt.0) 
     &             call file_close_delete(opti_stat%ffscr(iopt)%fhand)
          if (opti_stat%ffext(iopt)%fhand%unit.gt.0)
     &             call file_close_delete(opti_stat%ffext(iopt)%fhand)
        end do

        ! release all temporary memory
        ifree = mem_flushmark('leqevpc_temp')
        
        ! tell the slaves on which records the new test vectors are:
        nsub = opti_stat%ndim_rsbsp
        if (nsub.ne.opti_stat%ndim_vsbsp)
     &       call quit(1,'leq_evp_control',
     &       'different subspace dimensions?')
        nrequest = opti_stat%nadd
        do irequest = 1, nrequest
          irectrv(irequest) =
     &         opti_stat%iord_vsbsp(nsub-nrequest+irequest)
          irecmvp(irequest) =
     &         opti_stat%iord_rsbsp(nsub-nrequest+irequest)
          irecmet(irequest) =
     &         opti_stat%iord_ssbsp(nsub-nrequest+irequest)
        end do
        task = 4
          
      else 
*----------------------------------------------------------------------*
* clean up
*----------------------------------------------------------------------*
        do iopt = 1, opti_info%nopt
          if (opti_stat%ffscr(iopt)%fhand%unit.gt.0) 
     &             call file_close_delete(opti_stat%ffscr(iopt)%fhand)
          if (opti_stat%ffext(iopt)%fhand%unit.gt.0)
     &             call file_close_delete(opti_stat%ffext(iopt)%fhand)
        end do

        call leqevpc_cleanup()

        ! release all temporary memory
        ifree = mem_flushmark('leqevpc_temp')

        ! flush all other memory allocated by optimizer
        ifree = mem_flushmark('leqevpc_perm')

        task = 8 ! stop it

      end if

      if (ntest.ge.10) then
        write(lulog,*) 'at the end of optcont:'
        write(lulog,*) ' task = ',task
        write(lulog,*) ' iter:  ',iter
      end if

      call atim_csw(cpu,sys,wall)
      if (iprint.ge.5)
     &     call prtim(lulog,'time in optimizer',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
*----------------------------------------------------------------------*
*     a few internal subroutines follow
*----------------------------------------------------------------------*

      contains

*----------------------------------------------------------------------*
      subroutine leqevpc_mem(nincore,lenbuf,mem_free,nwfpar,nopt)
*----------------------------------------------------------------------*
*
*     get memory for buffers
*
*     internal subroutine: variable of main routine are global!
*
*----------------------------------------------------------------------*

      integer, intent(in) ::
     &     mem_free, nopt, nwfpar(nopt)
      integer, intent(out) ::
     &     nincore, lenbuf

      integer ::
     &     idx, len1, len2, len3, nmax_per_vec, nbatch,
     &     ifree

      ! currently, nopt==1 is the most often met case
      ! not sure whether this is optimal:
      nmax_per_vec = lblk_da
      do iopt = 1, nopt
        nmax_per_vec = max(nmax_per_vec,nwfpar(iopt))
      end do
      nincore = min(3,mem_free/nmax_per_vec)

      if (nincore.eq.3) then
        nbuf = 3
        len1 = nmax_per_vec
        len2 = nmax_per_vec
        len3 = nmax_per_vec
        lenbuf = nmax_per_vec
      else if (nincore.eq.2) then
        nbuf = 2
        len1 = nmax_per_vec
        len2 = nmax_per_vec
        len3 = 0
        lenbuf = nmax_per_vec
      else if (nincore.eq.1) then
        nbuf = 2
        len1 = nmax_per_vec
        len2 = mem_free-nmax_per_vec
        ! 2nd buffer should at least hold 20% of vector
        if (dble(len2)/dble(nmax_per_vec).lt.0.2d0) then
          nincore = 0
          len1 = mem_free/2
          len2 = mem_free/2
          len3 = 0
        end if
        lenbuf = len2
      else
        nbuf = 2
        len1 = mem_free/2
        len2 = mem_free/2
        len3 = 0
        lenbuf = len2
      end if

      ifree = mem_alloc_real(xbuf1,len1,'buffer_1')
      ifree = mem_alloc_real(xbuf2,len2,'buffer_2')
      if (nbuf.eq.3)
     &     ifree = mem_alloc_real(xbuf3,len3,'buffer_3')

      if (iprint.ge.5) then
        write(lulog,*) ' allocated ',nbuf,' buffers'
        write(lulog,*) ' # incore vectors: ',nincore
        write(lulog,*) ' total size of buffers: ',len1+len2+len3
        write(lulog,*) ' remaining core memory: ',ifree
        if (nincore.le.1) then
          nbatch = nmax_per_vec/lenbuf
          if (nbatch*lenbuf.lt.nmax_per_vec) nbatch = nbatch+1
          write(lulog,*) ' out-of-core routines need ',nbatch,' cycles'
        end if
      end if

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine leqevpc_init()
*----------------------------------------------------------------------*
*
*     initialize counters, open files and allocate some arrays
*
*----------------------------------------------------------------------*
      
      character ::
     &     name*64
      integer ::
     &     lenord, lenmat

* allocate some memory for subspace matrices
      if (iroute.eq.1) then

        allocate(opti_stat%ffrsbsp(opti_info%nopt),
     &       opti_stat%ffvsbsp(opti_info%nopt),
     &       opti_stat%ffssbsp(opti_info%nopt),
     &       opti_stat%ffscr(opti_info%nopt),
     &       opti_stat%ffext(opti_info%nopt))
        do iopt = 1, opti_info%nopt
          opti_stat%ffrsbsp(iopt)%fhand => me_mvp(iopt)%mel%fhand
          opti_stat%ffvsbsp(iopt)%fhand => me_trv(iopt)%mel%fhand
          opti_stat%ffscr(iopt)%fhand => me_scr(iopt)%mel%fhand
          opti_stat%ffext(iopt)%fhand => me_ext(iopt)%mel%fhand
          if (use_s(iopt))
     &         opti_stat%ffssbsp(iopt)%fhand => me_met(iopt)%mel%fhand
        end do

        ! space for Mred, and for each root: xred, RHSred, Sred
        lenord = opti_info%maxsbsp ! includes nopt-factor
        lenmat = 4*lenord**2
      else
        lenmat = 0
        lenord = 0
      end if

      opti_stat%mxdim_sbsp = lenord

      if (lenmat.gt.0) then
        ifree = mem_alloc_real(opti_stat%sbspmat,lenmat,'sbspmat')
        ifree = mem_alloc_int (opti_stat%iord_rsbsp,lenord,'iord_rsbsp')
        ifree = mem_alloc_int (opti_stat%iord_vsbsp,lenord,'iord_vsbsp')
        ifree = mem_alloc_int (opti_stat%iord_ssbsp,lenord,'iord_ssbsp')
        opti_stat%iord_rsbsp(1:lenord) = 0
        opti_stat%iord_vsbsp(1:lenord) = 0
        opti_stat%iord_ssbsp(1:lenord) = 0
      end if

* initialize variables - we start with nroot user-provided guess vectors:
      if (modestr(1:3).eq.'LEQ') then
        opti_stat%ndim_rsbsp = 0 !opti_info%nroot
        opti_stat%ndim_vsbsp = 0 !opti_info%nroot
        opti_stat%ndim_ssbsp = 0 !opti_info%nroot
      else
        opti_stat%ndim_rsbsp = opti_info%nroot
        opti_stat%ndim_vsbsp = opti_info%nroot
        opti_stat%ndim_ssbsp = opti_info%nroot
        do idx = 1, opti_info%nroot
          opti_stat%iord_rsbsp(idx) = idx
          opti_stat%iord_vsbsp(idx) = idx
          opti_stat%iord_ssbsp(idx) = idx
        end do
      end if

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine leqevpc_cleanup()
*----------------------------------------------------------------------*
*     release internal units
*----------------------------------------------------------------------*

      integer ::
     &     iopt

      if (iroute.ge.1) then
        deallocate(opti_stat%ffrsbsp,opti_stat%ffvsbsp,
     &             opti_stat%ffssbsp,opti_stat%ffscr,
     &             opti_stat%ffext)
      end if

      return
      end subroutine


      end


