*----------------------------------------------------------------------*
      subroutine leq_evp_control(modestr,iter,
     &                   task,conv,xrsnrm,xeig,
     &                   use_s,
     &                   nrequest,irectrv,irecmvp,irecmet,
     &                   ffopt,fftrv,ffmvp,ffmet,ffrhs,ffdia,
     &                   opti_info,opti_stat)
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
      include 'def_filinf.h'
      include 'def_file_array.h'
      include 'def_optimize_info.h'
      include 'def_optimize_status.h'
      include 'ifc_memman.h'
      
* parameters
      integer, parameter ::
     &     ntest = 00
      character, parameter ::
     &     name_alg_leq*10(0:1) =
     &     (/"  CONJGRAD","  SUBSPACE"/),
     &     name_alg_evp*10(0:1) =
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

      type(file_array), intent(in) ::
     &     ffopt(*), fftrv(*), ffmvp(*), ffmet(*), ffrhs(*), ffdia(*)
      
      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout) ::
     &     opti_stat

      real(8), intent(inout) ::
     &     xrsnrm(opti_info%nroot),
     &     xeig(opti_info%nroot,2)


*     buffers for incore/out-of-core work:
      integer ::
     &     nincore, nbuf, lenbuf
      real(8), pointer ::
     &     xbuf1(:), xbuf2(:), xbuf3(:)

      type(filinf) ::
     &     ffscr
      logical ::
     &     lexit, lconv
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
     &       call write_title(luout,wst_section,'LEQ solver')
        if (modestr(1:3).eq.'EVP')
     &       call write_title(luout,wst_section,'EVP solver')
      end if
      if (ntest.ge.10) then
        write(luout,*) 'entered leq_evp_control with:'
        write(luout,*) ' iter: ', iter
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
     &       write(luout,*) 'Optimization algorithm:    ',
     &       name_alg_leq(iroute)
        if (modestr(1:3).eq.'EVP')
     &       write(luout,*) 'Optimization algorithm:    ',
     &       name_alg_evp(iroute)
        write(luout,'(x,a,i10)')
     &                 'Max. number of iterations: ',opti_info%maxmacit
        write(luout,'(x,a,e10.2)')
     &                 'Threshold for residual:    ',opti_info%thrgrd(1)
      end if

      if (ntest.ge.10) then
        write(luout,*) 'our route is: ',iroute
      end if
*======================================================================*
*     unless this is the initial call:
*     get memory for buffers, decide on in-core/out-of-core
*======================================================================*
      ! set mark for temporarily allocated stuff
      ifree = mem_setmark('leqevpc_temp')
      if (iter.ne.0.or.modestr(1:3).eq.'LEQ') then
        call leqevpc_mem(nincore,lenbuf,
     &       ifree,opti_info%nwfpar,opti_info%nopt)
c        if (nincore.le.1) then
          call file_init(ffscr,'leqscr.da',ftyp_da_unf,lblk_da)
          call file_open(ffscr)
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

        ! request slaves to put first few trial vectors and
        ! corresponding Mv-product to first nroot records on the
        ! respective files:        
        if (modestr(1:3).eq.'LEQ') then
          call leqc_init(xrsnrm,iroute,
     &       ffopt,fftrv,ffmvp,ffrhs,ffdia,
     &       nincore,lenbuf,ffscr,
     &       xbuf1,xbuf2,xbuf3,
     &       opti_info,opti_stat)          
          opti_stat%nadd = opti_info%nroot ! ?? <-- check that for LEQ
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

        if (ffscr%unit.gt.0) call file_close_delete(ffscr)

        iter = 1

        ifree = mem_flushmark('leqevpc_temp')

* return and let the slaves work
        return

*======================================================================*
* here the actual LEQ/EVP solvers come:
*======================================================================*
      else 
        if (ntest.ge.10) then
          write(luout,*) 'getting new vectors ...'
        end if

        if (modestr(1:3).eq.'LEQ') then
          call leqc_core(iter,
     &         task,iroute,xrsnrm,
     &         ffopt,fftrv,ffmvp,ffrhs,ffdia,
     &         nincore,lenbuf,ffscr,
     &         xbuf1,xbuf2,xbuf3,
     &         opti_info,opti_stat)
        else
          call evpc_core(iter,
     &         task,iroute,xrsnrm,xeig,
     &         use_s,
     &         ffopt,fftrv,ffmvp,ffdia,
     &         nincore,lenbuf,ffscr,
     &         xbuf1,xbuf2,xbuf3,
     &         opti_info,opti_stat)
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
     &       write(luout,*) 'after iteration ',iter
        idx = 0
        do iopt = 1, opti_info%nopt
          do iroot = 1, opti_info%nroot
            idx = idx+1
            write(luout,'(x,2(a,e10.3),a,l)')
     &                      ' norm of residual:  ', xrsnrm(idx),
     &                           '   threshold:  ',
     &                                   opti_info%thrgrd(iopt),
     &                           '   converged:  ',
     &                   xrsnrm(idx).lt.opti_info%thrgrd(iopt)
          end do
        end do
      end if

      if (lconv)
     &       write(luout,'(x,a,i5,a)')
     &         'CONVERGED IN ',iter,' ITERATIONS'
      if (lconv) conv = .true.

      if (.not.lconv) iter = iter + 1

      if (.not.lconv.and.
     &       (iter.gt.opti_info%maxmacit)) then
        write(luout,*) 'NO CONVERGENCE OBTAINED'
        iter = iter - 1
        lexit = .true.
      end if

      if (.not.(lconv.or.lexit)) then
*----------------------------------------------------------------------*
* do some stuff for the next macro-iteration
*----------------------------------------------------------------------*

        if (ffscr%unit.gt.0) call file_close_delete(ffscr)

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
        call leqevpc_cleanup()

        if (ffscr%unit.gt.0) call file_close_delete(ffscr)

        ! release all temporary memory
        ifree = mem_flushmark('leqevpc_temp')

        ! flush all other memory allocated by optimizer
        ifree = mem_flushmark('leqevpc_perm')

        task = 8 ! stop it

      end if

      if (ntest.ge.10) then
        write(luout,*) 'at the end of optcont:'
        write(luout,*) ' task = ',task
        write(luout,*) ' iter:  ',iter
      end if

      call atim_csw(cpu,sys,wall)
      if (iprint.ge.5)
     &     call prtim(luout,'time in optimizer',
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
        write(luout,*) ' allocated ',nbuf,' buffers'
        write(luout,*) ' # incore vectors: ',nincore
        write(luout,*) ' total size of buffers: ',len1+len2+len3
        write(luout,*) ' remaining core memory: ',ifree
        if (nincore.le.1) then
          nbatch = nmax_per_vec/lenbuf
          if (nbatch*lenbuf.lt.nmax_per_vec) nbatch = nbatch+1
          write(luout,*) ' out-of-core routines need ',nbatch,' cycles'
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
     &       opti_stat%ffssbsp(opti_info%nopt))
        do iopt = 1, opti_info%nopt
          opti_stat%ffrsbsp(iopt)%fhand => ffmvp(iopt)%fhand
          opti_stat%ffvsbsp(iopt)%fhand => fftrv(iopt)%fhand
          if (use_s(iopt))
     &         opti_stat%ffssbsp(iopt)%fhand => ffmet(iopt)%fhand
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
     &             opti_stat%ffssbsp)
      end if

      return
      end subroutine


      end


