*----------------------------------------------------------------------*
      subroutine optcont(imacit,imicit,imicit_tot,
     &                   task,conv,
     &                   energy,xngrd,
     &                   me_opt,me_grd,me_dia,
     &                   me_trv,me_h_trv,
     &                   n_states,
     &                   me_special,nspecial,
     &                   fspc,nspcfrm,
     &                   opti_info,opti_stat,
     &                   orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*
* General optimization control routine for non-linear optimization.
*
* The optimization strategy is set in opti_info, see def_optinf.h
* for further comments.
*
* It holds track of the macro/micro-iterations and expects to be
* called with imacit/imicit set to zero at the first time.
*
* optcont() requests information via the itask parameter (bit flags):
*
*     1: calculate energy
*     2: calculate gradient/vectorfunction
*     4: calculate Hessian/Jacobian(H/A) times trialvector product
*     8: exit
*
* Files (now: ME-lists, i.e. files + some more info):
* a) passed to slave routines
*   me_opt   -- new set of wave-function parameters
*
*   for 2nd-order methods:
*   me_trv   -- additional trial-vector for H/A matrix-vector products
*              and dto. for left wave-function and orb-rotations
*
* b) passed from slave routines
*   me_grd   -- gradient or vectorfunction
*
*   2nd order
*   me_h_trv   -- H/A matrix-vector product
*              and dto. for ....
*
* c) own book keeping files
*   ffgrvfold -- previous gradient/vectorfunction
*   ff_h_trv_old  -- previous H/A matrix-vector product
*   ff_newstp -- current new step (scratch)
*   ffst_sbsp -- subspace of previous steps
*   ffgv_sbsp -- subspace of previous gradient/vector function diff.s
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ioparam.h'
      include 'def_file_array.h'
      include 'mdef_operator_info.h'
      include 'def_optimize_info.h'
      include 'def_optimize_status.h'
      include 'ifc_memman.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

* parameters
      integer, parameter ::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="optcont"
      character(len=10), parameter ::
     &     name_alg(0:3) =
     &     (/"Pert. Upd.","      DIIS","  ASSJ/RLE"," 2ND ORDER"/)

      integer, intent(out) ::
     &     task
      logical, intent(out) ::
     &     conv
      integer, intent(in) ::
     &     nspcfrm

      type(formula_item), intent(in) ::
     &     fspc(nspcfrm)

      type(orbinf),intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf),intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info

      integer, intent(inout) ::
     &     imacit, imicit, imicit_tot

      real(8), intent(inout) ::
     &     energy(0:*), xngrd(*)

      integer, intent(in) ::
     &     nspecial, n_states

      type(me_list_array), intent(in) ::
     &     me_opt(*), me_grd(*), me_dia(*), me_special(nspecial),
     &     me_trv(*), me_h_trv(*)

      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout), target ::
     &     opti_stat

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
     &     iprint, iroute, ifree, iopt, i_state
      real(8) ::
     &     de, cpu, sys, wall, cpu0, sys0, wall0

      iprint = max(ntest,iprlvl)
      call atim_csw(cpu0,sys0,wall0)

* be verbose?
      if (iprint.ge.5) then
        write(lulog,'(/,x,a,/,x,a)')
     &         'Optimization control',
     &         '===================='
      end if
      if (ntest.ge.10) then
        write(lulog,*) 'entered optcont with:'
        write(lulog,*) ' imacit, imicit, imicit_tot: ',
     &         imacit, imicit, imicit_tot
      end if

* set iroute:
* 0 simple perturbation (for testing, not recommended)
* 1 DIIS
* 2 RLE (ASSJ)
* 3 2ND-ORDER optimization
      !iroute = opti_info%mode_nleq
      iroute = 0

      if (iprint.ge.5.or.imacit.eq.0) then
        write(lulog,*) 'Optimization algorithm:    ',name_alg(iroute)
        write(lulog,'(x,a,i10)')
     &                 'Max. number of iterations: ',opti_info%maxmacit
        write(lulog,'(x,a,e10.2)')
     &                 'Threshold for residual:    ',opti_info%thrgrd(1)
      end if

      if (ntest.ge.10) then
        write(lulog,*) 'our route is: ',iroute
      end if
*======================================================================*
*     unless this is the initial call:
*     get memory for buffers, decide on in-core/out-of-core
*======================================================================*
      ! set mark for temporarily allocated stuff
      if (imacit.ne.0) then
        ifree = mem_setmark('optc_temp')
        call optc_mem(nincore,lenbuf,
     &       ifree,opti_info%nwfpar,opti_info%nopt,opti_info%max_incore)
        call file_init(ffscr,'optscr.da',ftyp_da_unf,lblk_da)
        if (nincore.le.1) then
          call file_open(ffscr)
        end if
      end if

      ! still, we assume this:
      conv = .false.
*======================================================================*
* check iteration number:
* zero --> init everything
*======================================================================*
      if (imacit.eq.0) then

        ! set mark for memory allocated until optimization ends
        ifree = mem_setmark('optc_perm')
        call optc_init()

* set task -- we want the energy in any way
        task = 1
* ... and the gradient/vector function
        task = task + 2

        imacit = 1

* return and let the slaves work
        return

*======================================================================*
* beginning of a macro iteration
*======================================================================*
      else if (imicit.eq.0) then
        if (ntest.ge.10) then
          write(lulog,*) 'macro iteration part entered'
        end if

        call optc_macit(imacit,imicit,imicit_tot,
     &       task,iroute,opti_info%nopt,
     &       me_opt,me_grd,me_dia,
     &       me_trv,me_h_trv,
     &       n_states,
     &       me_special,nspecial,
     &       nincore,lenbuf,ffscr,
     &       xbuf1,xbuf2,xbuf3,
     &       fspc,nspcfrm,energy,xngrd,
     &       opti_info,opti_stat,
     &       orb_info,op_info,str_info,strmap_info)

        de = opti_stat%energy_last - energy(0)

      else
*======================================================================*
* micro-iteration
*======================================================================*
        if (ntest.ge.10) then
          write(lulog,*) 'micro iteration part entered'
        end if

        call optc_micit()

      end if

      lexit = .false.
      lconv = .false.
      ! convergence check:
      !   end of iteration for 1st-order methods
      !   before first micro-iteration for 2nd-order methods
      if ((opti_info%norder.eq.1.and.imicit.eq.0) .or.
     &    (opti_info%norder.eq.2.and.imicit.eq.1)) then
*======================================================================*
* end of macro-iteration:
*  check convergence and max. iterations:
*======================================================================*
        lexit = .false.

        lconv = .true.
        do iopt = 1, opti_info%nopt
          lconv = lconv.and.xngrd(iopt).lt.opti_info%thrgrd(iopt)
        end do

        if (iprint.ge.5) then
          if (opti_info%norder.eq.1)
     &         write(lulog,*) 'after iteration ',imacit
          if (opti_info%norder.eq.2)
     &         write(lulog,*) 'after macro-iteration ',imacit
          do iopt = 1, opti_info%nopt
            write(lulog,'(x,2(a,e10.3),a,l)')
     &                      ' norm of gradient:  ', xngrd(iopt),
     &                           '   threshold:  ',
     &                                   opti_info%thrgrd(iopt),
     &                           '   converged:  ',
     &                    xngrd(iopt).lt.opti_info%thrgrd(iopt)
          end do
        end if

        if (lconv.and.opti_info%norder.eq.1)
     &         write(lulog,'(x,a,i5,a)')
     &         'CONVERGED IN ',imacit,' ITERATIONS'
        if (lconv.and.opti_info%norder.eq.2) then
          imicit_tot = imicit_tot-1
          write(lulog,'(x,a,i5,a,i6,a)')
     &         'CONVERGED IN ',imacit,' MACRO-ITERATIONS (',imicit_tot,
     &         ' MICRO-ITERATIONS)'
        end if
        if (lconv) conv = .true.
        if (lconv) imicit = 0

        if (.not.lconv) imacit = imacit + 1

        if (.not.lconv.and.
     &       (imacit.gt.opti_info%maxmacit.or.
     &       imicit_tot.gt.opti_info%maxmicit)) then
           call warn('non-linear solver','NO CONVERGENCE OBTAINED')
c          imacit = imacit - 1
c          imicit_tot = imicit_tot - 1
          imicit = 0
          lexit = .true.
        end if
      end if

      if (imicit.eq.0) then

        if (.not.(lconv.or.lexit)) then
*----------------------------------------------------------------------*
* do some stuff for the next macro-iteration
*----------------------------------------------------------------------*

c          if (opti_info%norder.eq.2.and.imicit.eq.0) itransf = 1
c          call optc_prepnext()
          opti_stat%energy_last = energy(0)
          task = 3

        else
*----------------------------------------------------------------------*
* clean up
*----------------------------------------------------------------------*
          call optc_cleanup()

          if (ffscr%unit.gt.0) call file_close_delete(ffscr)

          ! release all temporary memory
          ifree = mem_flushmark('optc_temp')

          ! flush all other memory allocated by optimizer
          ifree = mem_flushmark('optc_perm')

          task = 8 ! stop it

          ! that was it ....
          return

        end if

        if (ntest.ge.10) then
          write(lulog,*) 'at the end of optcont:'
          write(lulog,*) ' task = ',task
          write(lulog,*) ' imacit,imicit,imicit_tot: ',
     &         imacit,imicit,imicit_tot
        end if

      end if ! end-of-macro-iteration part

      if (ffscr%unit.gt.0) call file_close_delete(ffscr)

      ! release all temporary memory
      ifree = mem_flushmark('optc_temp')

      call atim_csw(cpu,sys,wall)
      if (iprint.ge.1)
     &     call prtim(lulog,'time in optimizer',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
*----------------------------------------------------------------------*
*     a few internal subroutines follow
*----------------------------------------------------------------------*

      contains

*----------------------------------------------------------------------*
      subroutine optc_mem(nincore,lenbuf,mem_free,nwfpar,nopt,
     &                    max_incore)
*----------------------------------------------------------------------*
*
*     get memory for buffers
*
*     internal subroutine: variable of main routine are global!
*     max_incore --> for checking only
*
*----------------------------------------------------------------------*
c      implicit none

      integer, intent(in) ::
     &     mem_free, nopt, nwfpar(nopt), max_incore
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
      if (max_incore.gt.0) then
        nincore = min(nincore,max_incore)
        write(lulog,*) ' restricting nincore to ',nincore
      end if

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
      if (nbuf.ne.3) xbuf3 => null()

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
      subroutine optc_init()
*----------------------------------------------------------------------*
*
*     initialize counters, open files and allocate some arrays
*
*----------------------------------------------------------------------*

c      implicit none

c      include 'stdunit.h'

      character ::
     &     name*64
      integer ::
     &     lenord, lenmat
      type(optimize_status), pointer ::
     &     c_opti_stat

      c_opti_stat => opti_stat

* open internal files
      do i_state = 1,n_states

      allocate(c_opti_stat%next_state)
      c_opti_stat%i_state = i_state

      if (iroute.ge.1.and.i_state.eq.1) then
        allocate(c_opti_stat%ffrsbsp(opti_info%nopt),
     &       c_opti_stat%ffvsbsp(opti_info%nopt))
        do iopt = 1, opti_info%nopt
          allocate(c_opti_stat%ffrsbsp(iopt)%fhand)
          allocate(c_opti_stat%ffvsbsp(iopt)%fhand)
          write(name,'("resd_sbsp",i2.2,".da")') iopt
          call file_init(c_opti_stat%ffrsbsp(iopt)%fhand,
     &         trim(name),ftyp_da_unf,lblk_da)
          write(name,'("vect_sbsp",i2.2,".da")') iopt
          call file_init(c_opti_stat%ffvsbsp(iopt)%fhand,
     &         trim(name),ftyp_da_unf,lblk_da)
          call file_open(c_opti_stat%ffrsbsp(iopt)%fhand)
          call file_open(c_opti_stat%ffvsbsp(iopt)%fhand)
        end do
      end if

* allocate some memory for subspace matrices
      if (iroute.eq.1) then
        lenmat = (opti_info%maxsbsp+1)*(opti_info%maxsbsp+2)/2
        lenord = opti_info%maxsbsp
      else if (iroute.eq.2) then
        lenmat = 2*opti_info%maxsbsp*opti_info%maxsbsp
        lenord = opti_info%maxsbsp
      else if (iroute.eq.3) then
        ! space for hred, cred, gred
        lenmat = (2*opti_info%micifac)**2 + 2*2*opti_info%micifac
        lenord = 2*opti_info%micifac
      else
        lenmat = 0
        lenord = 0
      end if

      c_opti_stat%mxdim_sbsp = lenord

      if (lenmat.gt.0) then
        ifree = mem_alloc_real(c_opti_stat%sbspmat,lenmat,'sbspmat')
        ifree = mem_alloc_int
     &       (c_opti_stat%iord_rsbsp,lenord,'iord_rsbsp')
        ifree = mem_alloc_int
     &       (c_opti_stat%iord_vsbsp,lenord,'iord_vsbsp')
        c_opti_stat%iord_rsbsp(1:lenord) = 0
        c_opti_stat%iord_vsbsp(1:lenord) = 0
      end if

* initialize variables
      c_opti_stat%ndim_rsbsp = 0
      c_opti_stat%ndim_vsbsp = 0

      c_opti_stat%trrad = opti_info%trini
      c_opti_stat%energy_last = 0d0

* init 2nd-order solver
      if (opti_info%norder.eq.2) then
        ! begin with Newton-eigenvector method ...
        c_opti_stat%mode_2nd = 2
        ! ... and a gamma of 1
        c_opti_stat%gamma = 1d0
      end if

      c_opti_stat => c_opti_stat%next_state

      end do

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine optc_cleanup()
*----------------------------------------------------------------------*
*     release internal units
*----------------------------------------------------------------------*
c      implicit none

      integer ::
     &     iopt

      type(optimize_status), pointer ::
     &     c_opti_stat, c_opti_stat_2

      if (iroute.ge.1) then
        do iopt = 1, opti_info%nopt
          call file_delete(opti_stat%ffvsbsp(iopt)%fhand)
          call file_delete(opti_stat%ffrsbsp(iopt)%fhand)
          deallocate(opti_stat%ffrsbsp(iopt)%fhand)
          deallocate(opti_stat%ffvsbsp(iopt)%fhand)
        end do
        deallocate(opti_stat%ffrsbsp,opti_stat%ffvsbsp)
      end if

      c_opti_stat => opti_stat
      do i_state=1,n_states
       c_opti_stat_2 => c_opti_stat%next_state
       deallocate(c_opti_stat)
       c_opti_stat => c_opti_stat_2
      end do
      deallocate(c_opti_stat)
      nullify(c_opti_stat_2)

      return
      end subroutine


      end


