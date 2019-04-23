*----------------------------------------------------------------------*
!>     solve eigenvalue problem  Mx = lambda x
!!
!!     the formula with label "label_form" describes how to calculate
!!     the matrix trial-vector products and the r.h.s.
!!
!!     @para nopt                  number of x operators to be solved for
!!                           in case of coupled equations
!!     @param nroots                number of roots per x operator
!!     @param targ_root             the target root per x operator
!!
!!     @param label_opt(nopt)       label of solution vectors
!!     @param label_prc(nopt)       label of preconditioners
!!     @param label_op_mvp(nopt)    label operators describing Mx-products
!!     @param label_op_met(nopt)    label operators describing Sx-products
!!                                  if S is unity, pass label of operator
!!                                   associated with ME-list label_opt
!!
!!     the latter two are used to initilize temporary ME-lists
!!
!!     @param op_info:   operator/ME-list definitions
!!     @param form_info: formula definitions
!!     @param str_info: string information (to be passed to subroutines)
!!     @param strmap_info: string mappings (to be passed to subroutines)
!!     @param orb_info: orbital space information (to be passed)
!!
!!     @param thr_suggest: allows weaker convergence threshold
!!
!!     'choice' is used to prefer a particular operator over other for the guess
!!     choice = 0 (default), consider all the operators involved
!!     choice = iopt, consider operator iopt, iopt=1,nopt, nopt is the number
!!              of operators involved
!!     init: it is used to initialise guesses for each of list_opt(s)
!!           from existing files instead of going through 'init_guess'.
!!           The same can also be done using the keyword 'calculate.solve.eigen.resume',
!!           But that would be imposed to all the calls of this solver during the process.
!----------------------------------------------------------------------*
      subroutine solve_evp2(mode_str,
     &     nopt,nroots,targ_root,label_opt,label_prc,label_op_mvp,
     &     label_op_met,label_form,
     &     label_special,nspecial,label_spcfrm,nspcfrm,thr_suggest,
     &     choice_opt,init,
     &     op_info,form_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
      implicit none             ! for sure
!     DAVIDSON algorithm: solves M*v+x*S*v=0 M,S are matrices, v is an (unknown) vector x is a scalar(the root)

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'def_optimize_info.h'
      include 'def_file_array.h'
      include 'def_optimize_status.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'mdef_formula_info.h'
      include 'def_dependency_info.h'
      include 'ifc_memman.h'
      include 'mdef_target_info.h'
      include 'ifc_input.h'
      include 'def_davidson_subspace.h'
      include "def_guess_gen.h"


      integer, parameter ::
     &     ntest = 000
      character(len=*),parameter::
     &     i_am="solve_evp2"

      integer, intent(in) ::
     &     nopt, nroots, nspecial, nspcfrm, targ_root,choice_opt
      character(*), intent(in) ::
     &     mode_str,
     &     label_opt(nopt),
     &     label_prc(nopt),
     &     label_op_mvp(nopt),
     &     label_op_met(nopt),
     &     label_special(nspecial),
     &     label_spcfrm(nspcfrm),
     &     label_form
      real(8), intent(in) ::
     &     thr_suggest
      type(formula_info) ::
     &     form_info
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf) ::
     &     orb_info
      logical, intent(in) ::
     &     init

      logical ::
     &     conv, use_s_t, use_s(nopt), trafo(nopt),use_init_guess(nopt),
     &     home_in, restart
      character(len_opname) ::
     &     label
      integer ::
     &     iter, iprint, task, ifree, iopt, jopt, nintm, irequest,
     &     nrequest, nvectors, iroot, idx, ierr, idxmel, nout,
     &     jdx,
     &     lenbuf, nincore,
     &     nlists
      real(8) ::
     &     trf_nrm,
     &     xresmax, xdum, xnrm,
     &     xeig(nroots,2),reig(nroots,2),
     &     xresnrm(nroots,nopt),
     &     xnrm2(nroots)
      type(me_list_array), pointer ::
     &     me_opt(:), me_dia(:),
     &     me_trv(:), me_mvp(:), me_mvpprj(:),
     &     me_mvort(:), me_vort(:),
     &     me_special(:), me_scr(:), me_home(:),
     &     me_met(:),me_metort(:), me_res(:),
     &     me_pnt_list(:)
      type(file_array), pointer ::
     &     ffdia(:), ff_trv(:),
     &     ffopt(:), ff_mvp(:), ff_met(:), ffspecial(:), ff_scr(:),
     &     ffhome(:), ff_ext(:)
      type(me_list), pointer ::
     &     me_pnt
      type(dependency_info) ::
     &     depend
      type(optimize_info) ::
     &     opti_info
      type(optimize_status) ::
     &     opti_stat
      type(formula), pointer ::
     &     form_mvp
      type(formula_item) ::
     &     fl_mvp, fl_spc(nspcfrm)

      type(guess_generator)::
     &     guess_gen
      integer, pointer ::
     &     irecmvp(:), irectrv(:), irecmet(:)
      real(8), pointer::
     &     xret(:), xbuf1(:), xbuf2(:), xbuf3(:), xoverlap(:),
     &     xrsnrm(:,:)

      character ::
     &     fname*256
      type(davidson_subspace_t)::
     &     dvdsbsp

      logical, external ::
     &     file_exists, generate_guess
      integer, external ::
     &     idx_formlist, idx_mel_list, idx_xret
      real(8), external ::
     &     da_ddot


      ifree = mem_setmark('solve_evp')

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'entered solve_evp')
        write(lulog,*) 'nopt   = ',nopt
        write(lulog,*) 'nroots = ',nroots
        write(lulog,*) 'targ_root = ',targ_root
      end if
      call print_solver_header_h(luout)
      if (lulog .ne. luout)  call print_solver_header_h(lulog)




      idx = idx_formlist(label_form,form_info)
      if (idx.le.0)
     &     call quit(1,i_am,
     &     'did not find formula '//trim(label_form))
      form_mvp => form_info%form_arr(idx)%form

      allocate(me_opt(nopt),
     &     me_dia(nopt),
     &     me_mvp(nopt),me_trv(nopt),me_mvpprj(nopt),
     &     me_mvort(nopt),me_vort(nopt),
     &     me_special(nspecial),me_scr(nopt),
     &     me_met(nopt),me_metort(nopt),me_res(nopt))
      allocate(ffopt(nopt),ffdia(nopt),
     &     ff_trv(nopt),ff_mvp(nopt),ff_met(nopt),ffspecial(nspecial),
     &     ff_scr(nopt),ff_ext(nopt) )
      allocate(xrsnrm(nroots,nopt))

      do iopt = 1, nopt
         ! lists to be optimized
         me_opt(iopt)%mel   => get_mel_h(label_opt(iopt), op_info)
         if (.not.mel_has_file( me_opt(iopt)%mel))
     &        call quit(1,i_am,
     &       'no file associated with list '//trim(label_opt(iopt)))
         ! preconditioner (diag because they might be inverted diagonals of Hamiltonian)
         me_dia(iopt)%mel   => get_mel_h(label_prc(iopt),op_info)
         if (.not.mel_has_file( me_dia(iopt)%mel))
     &        call quit(1,i_am,
     &        'no file associated with list '//trim(label_prc(iopt)))

      end do

      ! special lists needed?
      do idx = 1, nspecial
         me_special(idx)%mel => get_mel_h(label_special(idx), op_info)
         if (.not.mel_has_file( me_special(idx)%mel))
     &        call quit(1,i_am,
     &        'no file associated with list '//trim(label_special(idx)))
      end do



      do jdx = 1, nspcfrm
        idx = idx_formlist(label_spcfrm(jdx),form_info)
        if (idx.le.0)
     &       call quit(1,i_am,
     &       'did not find formula '//trim(label_spcfrm(jdx)))
        call read_form_list(form_info%form_arr(idx)%form%fhand,
     &                      fl_spc(jdx),.true.)
      end do

      call set_opti_info(opti_info,3,nopt,nroots,me_opt,mode_str)

      nvectors = opti_info%maxsbsp
      use_s_t = .false.
      trafo=.false.
      do iopt=1,nopt
        ! weaker convergence threshold requested?
         opti_info%thrgrd(iopt)=max(opti_info%thrgrd(iopt),thr_suggest)
        ! which operators use a metric?
         use_s(iopt) = trim(label_op_met(iopt)).ne.
     &        trim(me_opt(iopt)%mel%op%name)
         ! does any operator use a metric
         if (use_s(iopt)) use_s_t = .true.
      end do

      do iopt=1,nopt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     intermediate lists in normal space
         ! mvp = matrix vector product
         write(fname,'("mvp_",i3.3)') iopt
         me_mvp(iopt)%mel => me_from_template(
     &        fname, label_op_mvp(iopt), me_opt(iopt)%mel,
     &        nvectors,
     &        op_info,  orb_info, str_info, strmap_info )

         select case (opti_info%typ_prc(iopt))
         case(optinf_prc_spinp)
            me_mvpprj(iopt)%mel => me_special(1)%mel
         case(optinf_prc_prj)
           me_mvpprj(iopt)%mel => me_mvp(iopt)%mel
        case(optinf_prc_spinrefp)
           me_mvpprj(iopt)%mel => me_special(1)%mel
        case default
           me_mvpprj(iopt)%mel => me_mvp(iopt)%mel
        end select

        ! trv = trialvector
        write(fname,'("trv_",i3.3)') iopt
        me_trv(iopt)%mel => me_from_template(
     &       fname, me_opt(iopt)%mel%op%name, me_opt(iopt)%mel,
     &       nvectors,
     &       op_info,  orb_info, str_info, strmap_info)



! use of metric requested?
        if (use_s(iopt)) then
             ! get a ME list for metric-times-vector products
! (have same symmtry properties as result!)
!S * vector product -- svp
           write(fname,'("svp_",i3.3)') iopt
           me_met(iopt)%mel => me_from_template(
     &          fname, label_op_met(iopt) , me_trv(iopt)%mel,
     &          nvectors,
     &          op_info,  orb_info, str_info, strmap_info)
           ff_met(iopt)%fhand => me_met(iopt)%mel%fhand
        else
           me_met(iopt)%mel => me_trv(iopt)%mel
       end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     intermediate lists in possibly transformed space

! vort is used as template
! in case of ab-sym braking trafo, get sym props from special list
       if (opti_info%typ_prc(iopt).eq.optinf_prc_traf
     &      .and.nspecial.ge.3) then
          trafo(iopt)=.true.
          me_vort(iopt)%mel => me_special(1)%mel
       elseif (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
          trafo(iopt)=.true.
          me_vort(iopt)%mel => me_special(4)%mel
       else
          trafo(iopt)=.false.
          me_vort(iopt)%mel => me_opt(iopt)%mel
       end if

         ! scr = scratch used for mvp is orthogonal space
         write(fname,'("scr_",i3.3)') iopt
         me_scr(iopt)%mel => me_from_template(
     &        fname, me_opt(iopt)%mel%op%name, me_vort(iopt)%mel,
     &        nvectors,
     &        op_info,  orb_info, str_info, strmap_info)

        ! can use me_scr/me_mvp
        if (trafo(iopt))then
           me_mvort(iopt)%mel => me_scr(iopt)%mel
        else
           me_mvort(iopt)%mel => me_mvp(iopt)%mel
        end if

        ! metric in orthogonal space
        write(fname,'("metort_",i3.3)') iopt
        me_metort(iopt)%mel => me_from_template(
     &       fname, me_opt(iopt)%mel%op%name,me_vort(iopt)%mel,
     &       nvectors,
     &       op_info,  orb_info, str_info, strmap_info)

        ! residual (in orthogonal space)
        write(fname,'("res_",i3.3)') iopt
        me_res(iopt)%mel => me_from_template(
     &       fname, me_opt(iopt)%mel%op%name, me_vort(iopt)%mel,
     &       nvectors,
     &       op_info,  orb_info, str_info, strmap_info)

        ! assign me_vort to the list actually used in the calculation
        if (opti_info%typ_prc(iopt).eq.optinf_prc_traf
     &       .and.nspecial.ge.3) then
           trafo(iopt)=.true.
           me_vort(iopt)%mel => me_special(1)%mel
        elseif (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
           trafo(iopt)=.true.
           me_vort(iopt)%mel => me_special(4)%mel !redundant
        else
           trafo(iopt)=.false.
           me_vort(iopt)%mel => me_trv(iopt)%mel !changing
        end if
      end do ! iopt

      call dvdsbsp_init(dvdsbsp, opti_info%maxsbsp, me_vort, nopt,
     &     use_s)

      ! write information to opti_info about signs which occur
      ! in trv*mvp or trv*met  multiplication
      ! relevant when trv is njoined=1 op. but mvp (met) are njoined=2 op's
      call set_opti_info_signs(opti_info,3,nopt,
     &     me_trv,me_mvp,me_met,me_met,use_s)

      call print_settings(lulog, opti_info)

      do iopt=1, nopt
         call assign_me_list(me_opt(iopt)%mel%label,
     &        me_opt(iopt)%mel%op%name,op_info)
      end do
      ! read formula
      call read_form_list(form_mvp%fhand,fl_mvp,.true.)
      ! number of info values returned on xret
      nout = depend%ntargets
      allocate(xret(nout))

      ! records with trial vectors and Mv-products, needed in evp_control:
      ifree = mem_alloc_int(irectrv,nroots,'rectrv')
      ifree = mem_alloc_int(irecmvp,nroots,'recmvp')
      ifree = mem_alloc_int(irecmet,nroots,'recmet')

      use_init_guess = .true.
      home_in = .false.
      do iopt = 1, nopt
        ! open result vector file(s)
! if file already open, use as initial guess (if requested)!
        me_opt(iopt)%mel%fhand => me_opt(iopt)%mel%fhand
        if (me_opt(iopt)%mel%fhand%unit.gt.0.and.opti_info%resume
     &       .and.nroots.eq.1) then ! switched off due to frequent
! problems for nroots>1
          write(lulog,'(a,i4,a)')
     &         'Using last vector as initial guess (iopt =',iopt,')'
          use_init_guess(iopt) = .false.
        end if
        if (me_opt(iopt)%mel%fhand%unit.gt.0.and.nroots.gt.1) then
! copy this root so that we may home in on it later
          if (nopt.ne.1) call quit(1,i_am,
     &         'homing in available only for one opt. vector yet')
          home_in = .true.
          allocate(me_home(1),ffhome(1))
          if (opti_info%typ_prc(iopt).eq.optinf_prc_traf
     &         .and.nspecial.ge.3) then
            me_pnt => me_special(1)%mel
          elseif (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
            me_pnt => me_special(4)%mel
          else
            me_pnt => me_opt(iopt)%mel
          end if

          fname = "home"
          me_home(1)%mel => me_from_template(
     &         fname, me_opt(iopt)%mel%op%name,me_pnt,
     &         1,
     &         op_info,  orb_info, str_info, strmap_info)
          ffhome(1)%fhand => me_home(1)%mel%fhand
          call file_open(me_home(1)%mel%fhand)
          call assign_me_list(me_trv(iopt)%mel%label,
     &         me_opt(iopt)%mel%op%name,op_info)

          call list_copy(me_opt(iopt)%mel,me_home(1)%mel,.false.)
        end if
        call file_ensure_open(me_opt(iopt)%mel%fhand)
        call file_ensure_open(me_scr(iopt)%mel%fhand)
        call file_ensure_open(me_trv(iopt)%mel%fhand)
! open corresponding matrix vector products ...
        call file_ensure_open(me_mvp(iopt)%mel%fhand)
        call file_ensure_open(me_mvort(iopt)%mel%fhand)
        call file_ensure_open(me_metort(iopt)%mel%fhand)
        call file_ensure_open(me_res(iopt)%mel%fhand)
        if (associated(me_met(iopt)%mel))
     &       call file_ensure_open(me_met(iopt)%mel%fhand)
        ! ... and corresponding preconditioner(s)
        call file_ensure_open(me_dia(iopt)%mel%fhand)
      end do

      do idx = 1, nspecial
         call file_ensure_open(me_special(idx)%mel%fhand)
      end do
      do iopt=1, nopt
        call assign_me_list(me_trv(iopt)%mel%label,
     &       me_opt(iopt)%mel%op%name,op_info)
      end do

      ! set dependency info for submitted formula list
      call set_formula_dependencies(depend,fl_mvp,op_info)

      do iopt=1,nopt
         if (opti_info%typ_prc(iopt).eq.optinf_prc_traf
     &        .and.nspecial.ge.3) then
            trafo(iopt)=.true.
            me_vort(iopt)%mel => me_special(1)%mel
         elseif (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
            trafo(iopt)=.true.
            me_vort(iopt)%mel => me_special(4)%mel
         else
            trafo(iopt)=.false.
            me_vort(iopt)%mel => me_trv(iopt)%mel
         end if
      end do
      if (init) then
      ! get the initial amplitudes from files
        do iopt = 1,nopt
          call file_ensure_open(me_opt(iopt)%mel%fhand)
          inquire(file=trim(me_opt(iopt)%mel%fhand%name),exist=restart)
          if (.not.restart) call warn(i_am,
     &         'No amplitude file found for restart! Setting to zero.')
          if (restart) then
            write(lulog,'(x,a,i1,a)')
     &           'Using old amplitude file for vector ',iopt,'!'
            do iroot = 1, nroots
              call switch_mel_record(me_trv(iopt)%mel,iroot)
              call switch_mel_record(me_opt(iopt)%mel,iroot)
              call list_copy(me_opt(iopt)%mel,me_trv(iopt)%mel,.false.)
            enddo
          else
            do iroot = 1, nroots
              call switch_mel_record(me_trv(iopt)%mel,iroot)
              call zeroop(me_trv(iopt)%mel)
            enddo
          endif
        enddo

      else
! get initial amplitudes
        do iopt=1,nopt
          call assign_me_list(me_trv(iopt)%mel%label,
     &          me_opt(iopt)%mel%op%name,op_info)
        end do
        call init_guess_gen(guess_gen, max(1000,4*nroots),
     &       me_dia, me_vort, nopt,
     &       op_info,str_info,strmap_info, orb_info, opti_info)
        iroot = 1
        trial_gen_loop: do while(iroot .le.nroots)
        do iopt = 1,nopt
          call switch_mel_record(me_vort(iopt)%mel,iroot)
          call zeroop(me_vort(iopt)%mel)
          call switch_mel_record(me_trv(iopt)%mel,iroot)
          call touch_file_rec(me_trv(iopt)%mel%fhand)
          call zeroop(me_trv(iopt)%mel)
        end do
        if (.not.generate_guess(guess_gen, me_vort,
     &       iopt,nopt,choice_opt))
     &       call quit(1,i_am,
     &       'Could not find enough guess vectors')
        if (trafo(iopt))then
           if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
              call set_blks(me_vort(iopt)%mel,
     &             "P,H|P,V|V,H|V,V",0d0)
           end if
           call transform_back_wrap(fl_mvp,depend,
     &         me_special,me_vort(iopt)%mel,me_trv(iopt)%mel,
     &         trf_nrm,
     &         iopt,nspecial,
     &         me_trv(iopt)%mel,
     &         op_info, str_info, strmap_info, orb_info, opti_info)
! guess vectors of wrong spin symmetry and twinned guess vectors will be discarded
           if (should_discard_vector_h(trf_nrm, nopt,nroots,
     &           iroot,xnrm2,me_trv, opti_info) )then
              cycle trial_gen_loop
           else
              xnrm2(iroot) = trf_nrm**2
           end if
        end if

        if (opti_info%typ_prc(iopt).eq.optinf_prc_spinp.or.
     &          opti_info%typ_prc(iopt).eq.optinf_prc_prj.or.
     &          opti_info%typ_prc(iopt).eq.optinf_prc_spinrefp) then
             call  spin_proj_h(opti_info%typ_prc(iopt),
     &            opti_info%nwfpar(iopt),me_trv(iopt)%mel,me_special,
     &            nspecial,0,1,
     &            fl_spc, xnrm,
     &            opti_info,orb_info, op_info, str_info, strmap_info)
             if (xnrm**2.lt.1d-12) then
               if (iprlvl.ge.5) write(lulog,*)
     &              'Discarding guess vector due to projection.'
               me_trv(iopt)%mel%fhand%last_mod(iroot) = -1
               cycle trial_gen_loop
             else if (iroot.gt.1) then
               call orthogonalize_roots_h(iroot, opti_info%nwfpar(iopt),
     &              xnrm2,xnrm**2,me_trv(iopt)%mel)
             else
               xnrm2(iroot) = xnrm**2
             end if
        end if
! Assumption: If the controlflow reaches here, a new guess vector was created
        iroot = iroot +1
       end do trial_gen_loop
       call del_guess_gen(guess_gen)
      endif
      deallocate(xret)
      ! number of info values returned on xret
      nout = depend%ntargets
      allocate(xret(nout))

      do iroot=1,nroots
         do iopt=1,nopt
            call switch_mel_record(me_trv(iopt)%mel,iroot)
            call touch_file_rec(me_trv(iopt)%mel%fhand)
         end do
      end do

      call init_buffers(opti_info%nwfpar, nopt,
     &     xbuf1,xbuf2,xbuf3,nincore, lenbuf)

      iter = 0
      task = 4
      nrequest=nroots
      xrsnrm=0d0
      xeig=0d0
      reig=0d0

      opt_loop: do while(task.lt.8)
         iter=iter+1
         if (iter.gt.1) then
            call print_step_results(iter-1,xrsnrm, xeig,
     &           nroots, nopt)
         end if


! 4 - get residual
        if (iand(task,4).eq.4) then
!   outside loop over requested Mv-products
           do irequest = 1, nrequest
              do iopt = 1, nopt

                 call assign_me_list(me_trv(iopt)%mel%label,
     &                me_trv(iopt)%mel%op%name,op_info)
                 call assign_me_list(me_mvp(iopt)%mel%label,
     &                label_op_mvp(iopt),op_info)

                 call switch_mel_record(me_trv(iopt)%mel,irequest)
                 call switch_mel_record(me_mvp(iopt)%mel,irequest)
                 if (use_s(iopt))
     &                call switch_mel_record(me_met(iopt)%mel,
     &                irequest)

              ! enforce MS-combination symmetry of trial vectors
! (if requested)
                 if (me_trv(iopt)%mel%absym.ne.0)then
                     call sym_ab_list(
     &                0.5d0,me_trv(iopt)%mel,me_trv(iopt)%mel,
     &                xdum,.false.,
     &                op_info,str_info,strmap_info,orb_info)
                  end if
                 call touch_file_rec(me_trv(iopt)%mel%fhand)
              end do



              call frm_sched(xret,fl_mvp,depend,0,0,
     &             .true.,.false.,op_info,str_info,strmap_info,orb_info)

            ! apply sign-fix (if needed)
            do iopt = 1, nopt
               call optc_fix_signs2(me_mvp(iopt)%mel%fhand,
     &              irequest,
     &              opti_info,iopt,
     &              opti_info%nwfpar(iopt),xbuf1)
               if (use_s(iopt))
     &              call optc_fix_signs2(me_met(iopt)%mel%fhand,
     &              irequest,
     &              opti_info,iopt,
     &              opti_info%nwfpar(iopt),xbuf1)
            end do
c dbg
c            write(lulog,*) 'output for request: ',irequest
c            call wrt_mel_file(lulog,5,me_mvp(1)%mel,
c     &           1,me_mvp(1)%mel%op%n_occ_cls,
c     &           str_info,orb_info)
c dbg

! enforce MS-combination symmetry of Mv-products:
! (if requested)
            do iopt = 1, nopt
               if (me_mvp(iopt)%mel%absym.ne.0)
     &              call sym_ab_list(
     &              0.5d0,me_mvp(iopt)%mel,me_mvp(iopt)%mel,
     &              xdum,.false.,
     &              op_info,str_info,strmap_info,orb_info)

              ! project out spin contaminations?
               if (opti_info%typ_prc(iopt).eq.optinf_prc_spinp.or.
     &              opti_info%typ_prc(iopt).eq.optinf_prc_prj.or.
     &              opti_info%typ_prc(iopt).eq.optinf_prc_spinrefp) then
                  if (nincore .lt. 2)
     &                 call quit(0,i_am, "not enough memory")
!     assign op. with list containing the mvp vector
                  call assign_me_list(me_mvp(iopt)%mel%label,
     &                 me_opt(iopt)%mel%op%name,op_info)
                  if (opti_info%typ_prc(iopt).eq.optinf_prc_spinp) then
                 call spin_project(me_mvp(iopt)%mel,me_mvpprj(iopt)%mel,
     &                    fl_spc(1),opti_info%nwfpar(iopt),
     &                    xbuf1,xbuf2,.false.,xnrm,
     &                    opti_info,orb_info,
     &                    op_info,str_info,strmap_info)
                  elseif (opti_info%typ_prc(iopt).eq.
     &                  optinf_prc_spinrefp) then
                 call spin_project(me_mvp(iopt)%mel,me_mvpprj(iopt)%mel,
     &                    fl_spc(2),opti_info%nwfpar(iopt),
     &                    xbuf1,xbuf2,.false.,xnrm,
     &                    opti_info,orb_info,
     &                    op_info,str_info,strmap_info)
                     call evaluate2(fl_spc(1),.false.,.false.,
     &                  op_info,str_info,strmap_info,orb_info,
     &                  xnrm,.false.)
                  else
                     call evaluate2(fl_spc(1),.false.,.false.,
     &                    op_info,str_info,strmap_info,orb_info,
     &                    xnrm,.false.)
                  end if
!     reassign lists to correct ops
                  call assign_me_list(me_trv(iopt)%mel%label,
     &                 me_opt(iopt)%mel%op%name,op_info)
                  call assign_me_list(me_mvp(iopt)%mel%label,
     &                 label_op_mvp(iopt),op_info)
               end if
            end do
         end do
      end if
      !
      call transform_forward_h(
     &     me_mvp, me_mvort,
     &     me_met, me_metort,
     &     me_trv,
     &     nopt, nrequest,
     &     me_special, nspecial,
     &     fl_mvp, depend,
     &     xrsnrm,
     &     trafo, use_s,
     &     op_info, str_info, strmap_info, orb_info, opti_info)

        call davidson_driver(
     &     dvdsbsp,
     &     iter, task, nrequest,
     &     opti_info%nroot, nopt,
     &     trafo, use_s,
     &     xrsnrm , xeig, reig,
     &     me_opt,me_dia,
     &     me_metort,
     &     me_scr,me_res,
     &     me_trv,me_mvp,me_vort,me_mvort,
     &     me_special, nspecial,
     &     xbuf1,xbuf2, xbuf3, nincore,lenbuf,
     &     fl_mvp,depend,
     &     fl_spc,nspcfrm,
     &     opti_info, opti_stat,
     &     orb_info, op_info, str_info,strmap_info
     &     )
        if (iand(task,8).eq.8)then
          me_pnt_list=> me_opt
          nlists = min(nroots,mel_get_maxrec_h(me_opt(1)%mel))
        else
          me_pnt_list=>me_trv
          nlists = nrequest
        end if
       call transform_back_h(fl_mvp,depend,
     &      trafo,
     &      me_special,
     &      me_vort, me_pnt_list,    !vort -> trv !new_trialvector created
     &      nopt, nlists,
     &      nspecial,
     &      me_trv,
     &      op_info,str_info,
     &      strmap_info, orb_info,
     &      opti_info)



        if (iand(task,8).eq.8)
     &       call print_step_results(iter,
     &       xrsnrm, xeig,nroots, nopt)

      end do opt_loop


      if(nrequest.eq.0)then
         write(lulog,'(x,a,i5,a)')
     &        'CONVERGED IN ',iter,' ITERATIONS'
         if (luout.ne.lulog.and.iprlvl.ge.5)
     &        write(luout,'(x,a,i5,a)')
     &        'CONVERGED IN ',iter,' ITERATIONS'
      else
         write(lulog,'(x,a,i5,a)') "Stopping after",iter,"iterations"
         call warn('linear solver', 'NO CONVERGENCE OBTAINED')
      end if

      do iopt = 1, nopt



!     remove the temporary lists
         write(fname,'("mvp_",i3.3)') iopt
         call me_ensure_deleted(fname,op_info)
         write(fname,'("trv_",i3.3)') iopt
         call me_ensure_deleted(fname,op_info)
         write(fname,'("svp_",i3.3)') iopt
         call me_ensure_deleted(fname,op_info)
         write(fname,'("res_",i3.3)') iopt
         call me_ensure_deleted(fname,op_info)
         write(fname,'("scr_",i3.3)') iopt
         call me_ensure_deleted(fname,op_info)
         write(fname,'("metort_",i3.3)') iopt
         call me_ensure_deleted(fname,op_info)

         call assign_me_list(label_opt(iopt),
     &        me_opt(iopt)%mel%op%name,op_info)




c        ! solution vector has been updated (if we had some iteration)

         call touch_file_rec(me_opt(iopt)%mel%fhand)

         if (home_in) then
!     home in on root with largest overlap with prior solution
            ifree = mem_setmark('solve_evp.home_in')
            ifree = mem_alloc_real(xoverlap,nroots,'xoverlap')
            ifree = mem_alloc_real(xbuf1,opti_info%nwfpar(iopt),'xbuf1')
            ifree = mem_alloc_real(xbuf2,opti_info%nwfpar(iopt),'xbuf2')
            xresmax = 0d0
            do iroot = 1, nroots
               xoverlap(iroot) = da_ddot(me_home(1)%mel%fhand,1,
     &              me_opt(iopt)%mel%fhand,iroot,
     &                                opti_info%nwfpar(iopt),
     &                                xbuf1,xbuf2,
     &              opti_info%nwfpar(iopt))
               xoverlap(iroot) = abs(xoverlap(iroot))
               if (xoverlap(iroot).gt.xresmax) then
                  idx = iroot
                  xresmax = xoverlap(iroot)
               end if
c     dbg
c     print *,'root / overlap: ',iroot,xoverlap(iroot)
c dbgend
            end do
            ifree = mem_flushmark()
            if (idx.ne.targ_root) then
               write(lulog,'(a,i4,a,f8.4)')
     &              'Homing in on root ',idx,' with overlap ',xresmax
! Interchange this record and the current record
! and leave everything else unchanged (a bit dirty)
            call switch_mel_record(me_opt(iopt)%mel,targ_root)
            call list_copy(me_opt(iopt)%mel,me_home(1)%mel,.false.)
            call switch_mel_record(me_opt(iopt)%mel,idx)
            call list_copy(me_home(1)%mel,me_opt(iopt)%mel,.true.)
            call switch_mel_record(me_opt(iopt)%mel,targ_root)
            call list_copy(me_home(1)%mel,me_opt(iopt)%mel,.false.)
         end if
         call del_me_list(me_home(1)%mel%label,op_info)
         deallocate(me_home,ffhome)
      end if

      end do

      do idx = 1, nspecial
         if (me_special(idx)%mel%fhand%unit.gt.0)
     &        call file_close_keep(me_special(idx)%mel%fhand)
      end do

      call print_roots(lulog, xrsnrm, nroots, nopt, xeig)
      if (lulog.ne.luout.and.iprlvl.ge.10)
     &     call print_roots(luout, xrsnrm, nroots, nopt, xeig)

      ! switch to target root if possible

      if (targ_root.ge.me_opt(1)%mel%fhand%active_records(1).and.
     &     targ_root.le.me_opt(1)%mel%fhand%active_records(2))
     &     call switch_mel_record(me_opt(1)%mel,targ_root)

      call clean_formula_dependencies(depend)

      ! note that only the pointer array ffopt (but not the entries)
      ! is deallocated:
      deallocate(me_opt,me_dia,me_trv,me_mvp,me_met,me_special,me_scr,
     &     me_vort, me_mvort, me_res)
      deallocate(ff_trv,ff_mvp,ffdia,ffopt,ff_met,xret,ffspecial,
     &     ff_scr,ff_ext)

!     not freeing xrsnrm is a memory leak.
!     but deallocating it results in a double free error.
!     I suspect an implicit reallocation somwhere with an incorrect intent(out)

!     deallocate(xrsnrm)

      call dealloc_formula_list(fl_mvp)
      do jdx = 1, nspcfrm
         call dealloc_formula_list(fl_spc(jdx))
      end do

      call dvdsbsp_del(dvdsbsp)
      ifree = mem_flushmark('solve_evp')

      return

      contains
*----------------------------------------------------------------------*
!>
!!
*----------------------------------------------------------------------*
      subroutine print_step_results(iter,xrsnrm,xeig,nroots,nopt)
      implicit none
      ! lulog,luout and iprlvl from include from parent
      integer, intent(in)::
     &     nroots,
     &     nopt,
     &     iter


      real(8),Dimension(nroots,nopt),intent(in)::
     &     xrsnrm
      real(8)::
     &     xeig(nroots,2)

      real(8)::
     &     xresmax

      real(8),external::
     &     fndmnx

      integer::
     &     iroot,
     &     idx


      xresmax = fndmnx(RESHAPE(xrsnrm, (/nroots*nopt/)),
     &     nroots*nopt,2)

      write(lulog,'("E>>",i3,24x,x,g10.4)') iter,xresmax
      if (lulog.ne.luout)
     &     write(luout,'("   ",i3,24x,x,g10.4)') iter,xresmax
      if (iprlvl.gt.0) then
         do iroot = 1, nroots
            write(lulog,'(" E>",3x,f24.12,x,3g10.4)')
     &           xeig(iroot,1),(xrsnrm(iroot,idx),
     &           idx = 1, nopt)
            if (lulog.ne.luout.and.iprlvl.ge.5)
     &           write(luout,'("   ",3x,f24.12,x,3g10.4)')
     &           xeig(iroot,1),(xrsnrm(iroot,idx),
     &           idx = 1, nopt)
         end do
      end if
      return
      end subroutine

*----------------------------------------------------------------------*
!>
*----------------------------------------------------------------------*
      subroutine print_roots(lu, xresnrm, nroots, nopt, xeig)
*----------------------------------------------------------------------*
      integer, intent(in) :: lu, nroots, nopt
      real(8)::
     &     xresnrm(nroots,nopt),
     &     xeig(nroots,2)

      call write_title(lu,wst_title,
     &     'Results for '//trim(label_opt(1)))
      write(lu,'("E>>",66("="))')
      write(lu,'("E>>",2x,'//
     &     '"root     eigenvalue (real)       eigenvalue (img.)'//
     &     '  |residual|")')
      write(lu,'("E>>",66("-"))')
      do iroot = 1, nroots
         write(lu,'("E>>",2x,i3,x,f22.12,2x,f22.12,2x,x,g10.4)')
     &        iroot,xeig(iroot,1),xeig(iroot,2),xresnrm(iroot,1)

c dbg
c        if (xresnrm(iroot).lt.opti_info%thrgrd(iopt)
c    &                  .and.abs(xeig(iroot,2)).lt.1d-12) then
c         do iopt=1,nopt
c           call switch_mel_record(me_opt(iopt)%mel,iroot)
c           call wrt_mel_file(lulog,5,me_opt(iopt)%mel,
c    &             1,me_opt(iopt)%mel%op%n_occ_cls,
c    &             str_info,orb_info)
c         enddo
c        end if
c dbg
      end do
      write(lu,'("E>>",66("="))')

      end subroutine
*----------------------------------------------------------------------*
!>    create an me_list from a template list
!!  returns a pointer to the created me_list
*----------------------------------------------------------------------*
      function me_from_template(label, label_op, me_template, nrec,
     &     op_info, orb_info, str_info, strmap_info)
*----------------------------------------------------------------------*
      implicit none

      character(*), intent(in) ::
     &     label, label_op
      integer,intent(in)::
     &     nrec
      type(me_list),intent(in)::
     &     me_template
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(inout) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(me_list),pointer::
     &     me_from_template
      integer,external::
     &     idx_mel_list


      call define_me_list(label,label_op,
     &     me_template%absym,
     &     me_template%casym,
     &     me_template%gamt,
     &     me_template%s2,
     &     me_template%mst,.false.,
     &     -1,1,nvectors,0,0,0,
     &     op_info,orb_info,str_info,strmap_info)
      me_from_template   => get_mel_h(label, op_info)

      return
      end function
*----------------------------------------------------------------------*
!!    resolves a label to a pointer to the matrix element list object
!!
*----------------------------------------------------------------------*
      function get_mel_h(label, op_info)
*----------------------------------------------------------------------*
      implicit none
      character(*),intent(in) ::
     &     label
      type(operator_info), intent(inout) ::
     &     op_info

      type(me_list),pointer::
     &     get_mel_h
      integer,external::
     &     idx_mel_list
      integer::
     &     idxmel
      idxmel = idx_mel_list(label,op_info)
      if (idxmel.le.0)
     &       call quit(1,i_am,
     &       'did not find list '//trim(label_opt(iopt)))
      get_mel_h=> op_info%mel_arr(idxmel)%mel
      return
      end function
      pure  function mel_has_file(mel)
      implicit none
      type(me_list),intent(in)::
     &     mel
      logical::
     &     mel_has_file
      mel_has_file=associated(mel%fhand)
      return
      end function
      subroutine file_ensure_open(fhand)

      type(filinf)::
     &     fhand
      if(fhand%unit.le.0)then
         call file_open(fhand)
      end if
      return
      end subroutine

*----------------------------------------------------------------------*
*
*     get memory for buffers
*
*
*
*----------------------------------------------------------------------*
      subroutine init_buffers(nwfpar,nopt,
     &     xbuf1,xbuf2,xbuf3, nincore, lenbuf)
      implicit none

      integer,parameter::
     &     ntest=00

      character(len=*),parameter::
     &     i_am="init_buffers"

      integer, intent(in) ::
     &     nopt, nwfpar(nopt)
      integer, intent(out) ::
     &     nincore, lenbuf

      real(8),Dimension(:),pointer,intent(out)::
     &   xbuf1,xbuf2,xbuf3

      integer ::
     &     idx, len1, len2, len3, nmax_per_vec, nbatch,
     &     ifree, mem_free, nbuf

      ifree=mem_setmark('solve_evp_buffer')
      mem_free=ifree*0.9        !allowed to use 90% of memory

      nmax_per_vec = lblk_da    !blocklength in direct access file (from include)
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
      end subroutine
      subroutine del_buffers(ifree)
*--------------------------------------------------------------*
      integer,intent(out)::
     &     ifree
      ifree=mem_flushmark('solve_evp_buffer')
      end subroutine
*--------------------------------------------------------------*
!> ensure that a me-list is deleted
!!
!!
*--------------------------------------------------------------*

      subroutine me_ensure_deleted(label, op_info )
      implicit none
      character(len=*),intent(in)::
     &     label
      type(operator_info), intent(inout) ::
     &     op_info
      type(me_list),pointer::
     &     me_pnt
      integer,external::
     &     idx_mel_list
      integer::
     &     idxmel

      idxmel = idx_mel_list(label,op_info)
      if (idxmel .gt. 0 ) then
            me_pnt => get_mel_h(label,op_info)
         if (has_perm_file(me_pnt, me_special,nspecial) ! Attention these are accessed by host association
     &        .or. has_perm_file(me_pnt, me_opt, nopt) )then
            call detach_mel_file(me_pnt,.false.)
         end if
         call del_me_list(label,op_info)
      end if
      end subroutine
      function has_perm_file(mel, mels, nmels)
      implicit none
      type(me_list) ,intent(in)::
     &     mel
      type(me_list_array),intent(in)::
     &     mels(*)
      integer,intent(in)::nmels
      logical::has_perm_file
      integer::ii
      has_perm_file=.false.
      do ii=1,nmels
         has_perm_file=has_perm_file
     &        .and. .not.
     &        ( mel%fhand%name
     &          .eq.mels(ii)%mel%fhand%name)
      end do
      end function
*----------------------------------------------------------------------*
      subroutine print_settings(lu,opti_info)
*----------------------------------------------------------------------*
      implicit none

      type(optimize_info), intent(in) ::
     &     opti_info
      integer,intent(in)::lu
      character(len=*), parameter ::
     &     name_alg_evp(0:1) =
     &     (/"  xxxxxxxx","  DAVIDSON"/)
      !TODO list the mode
      write(lu,*) 'Optimization algorithm:    ',
     &     name_alg_evp(opti_info%mode_evp)
      write(lu,'(x,a,i10)')
     &     'Max. number of iterations: ',opti_info%maxmacit
      write(lu,'(x,a,e10.2)')
     &     'Threshold for residual:    ',opti_info%thrgrd(1)
      write(lu,'(x,a,3i10)')
     &     'Number of parameters:      ',
     &     opti_info%nwfpar(1:opti_info%nopt)
      write(lu,'(x,a,3i10)')
     &     'modes      ',
     &     opti_info%typ_prc(1:opti_info%nopt)


      end subroutine
*----------------------------------------------------------------------*
!!   writes a header for the solver
!
*----------------------------------------------------------------------*
      subroutine print_solver_header_h(lu)
*----------------------------------------------------------------------*
      integer, intent(in)::
     &     lu

      write(lu,*) "This is the 'NEW' solver by Arne Bargholz 2017"
      end subroutine

*----------------------------------------------------------------------*
!!    transforms mvp and metric into orthogonal state
*----------------------------------------------------------------------*
      subroutine transform_forward_h(
     &     me_mvp, me_mvort,
     &     me_met, me_metort,
     &     me_trv,
     &     nopt, nnew,
     &     me_special, nspecial,
     &     flist, depend,
     &     xrsnrm,
     &     trafo, use_s,
     &     op_info, str_info, strmap_info, orb_info, opti_info)
*----------------------------------------------------------------------*
      implicit none

      integer,intent(in)::
     &     nnew, nopt,
     &     nspecial

      logical, intent(in)::
     &     trafo(nopt), use_s(nopt)
      type(me_list_array), intent(inout) ::
     &     me_mvp(nopt), me_mvort(nopt),
     &     me_met(nopt), me_metort(nopt),
     &     me_special(nspecial),
     &     me_trv(nopt)                  ! needed as target for transformation formula

      real(8),intent(inout)::
     &     xrsnrm(nnew,nopt)
      type(formula_item), intent(inout) ::
     &     flist

      type(dependency_info) ,intent(in)::
     &     depend

      type(optimize_info), intent(in) ::
     &     opti_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in)::
     &     strmap_info

      integer ::
     &     iroot,iopt
      do iroot = 1, nnew
         do iopt = 1, nopt
            call switch_mel_record(me_mvp(iopt)%mel,iroot)
            call switch_mel_record(me_mvort(iopt)%mel,iroot)
            call switch_mel_record(me_trv(iopt)%mel,iroot)
            if (trafo(iopt)) then
               call transform_forward_wrap(flist,depend,
     &              me_special,me_mvp(iopt)%mel,me_mvort(iopt)%mel, !mvp-> mvort
     &              xrsnrm(iroot,iopt),
     &              iopt, nspecial,
     &              me_trv(iopt)%mel,
     &              op_info, str_info, strmap_info, orb_info, opti_info)
               if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
                call set_blks(me_Mvort(iopt)%mel,"P,H|P,V|V,H|V,V",0d0)
               endif

               if (use_s(iopt))then
                  call switch_mel_record(me_met(iopt)%mel,iroot)
                  call switch_mel_record(me_metort(iopt)%mel,iroot)
                  call transform_forward_wrap(flist,depend,
     &                 me_special,me_met(iopt)%mel,me_metort(iopt)%mel, !met-> metort
     &                 xrsnrm(iopt,iroot),
     &                 iopt, nspecial,
     &                 me_trv(iopt)%mel,
     &                op_info,str_info,strmap_info, orb_info, opti_info)
                 if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
                    call set_blks(me_metort(iopt)%mel,
     &                   "P,H|P,V|V,H|V,V",0d0)
                 end if
               else
                  me_metort(iopt)%mel=> null()
               end if
           else
              if (use_s(iopt))then
                 me_metort(iopt)%mel=> me_met(iopt)%mel
              else
                 me_metort(iopt)%mel=>null()
              end if
           end if
        end do !iopt
      end do !iroot

      return
      end subroutine

*----------------------------------------------------------------------*
!>    orthogonalizes the newest root (iroot) to all previous roots
*----------------------------------------------------------------------*
      subroutine orthogonalize_roots_h(iroot,len_op,xnrm2,xnrm,me_root)
*----------------------------------------------------------------------*
      implicit none
      integer, intent(in)::
     &     len_op
      integer, intent(inout)::
     &     iroot
      type(me_list)::
     &     me_root

      integer::
     &     ifree,jroot,len_buf
      real(8)::
     &     xnrm2(*),xnrm

      real(8),pointer::
     &     xbuf1(:),xbuf2(:)
      real(8)::
     &     xover,fac1,fac2

      ifree = mem_setmark('init_guess.check_guess')
      len_buf = len_op
      ifree = mem_alloc_real(xbuf1,len_buf,
     &     'xbuf1')
      ifree = mem_alloc_real(xbuf2,len_buf,
     &     'xbuf2')

      do jroot = 1,iroot-1
         xover = da_ddot(
     &        me_root%fhand,jroot,
     &        me_root%fhand,iroot,
     &        len_op,
     &        xbuf1, xbuf2,
     &        len_buf)
! if the two are to similar then |j * i|/|j|*|i|=1   ==>   (j * i)^2/(j^2 * i^2) = 1    ==>    (j * i)^2/j^2 - i^2 = 0
! xover =  j* i  , xnrm2(j) = j^2 ,   xnrm = i^2
         if(abs(xover**2/xnrm2(jroot)-xnrm).lt.1d-6)then
            if (iprlvl.ge.5) write(lulog,*)
     &           'Discarding redundant guess vector.'
            me_root%fhand%last_mod(iroot) = -1
            iroot = iroot - 1
            exit
         else if (abs(xover).ge.1d-12) then
! remove this component using norm-conserving factors
            if (iprlvl.ge.5) write(lulog,*)
     &           'Removing overlap to root',jroot,'(',xover,')'
            fac1 = 1d0/sqrt(1d0-xover**2/(xnrm2(jroot)*xnrm))
            fac2 = -sign(1d0/sqrt((xnrm2(jroot)/xover)**2
     &           -xnrm2(jroot)/xnrm),xover)
            call da_vecsum(me_root%fhand,iroot,
     &           me_root%fhand,iroot,fac1,
     &           me_root%fhand,jroot,fac2,
     &           len_op,
     &           xbuf1,xbuf2,
     &           len_buf)
c                  xnrm = xnrm - xover**2/xnrm2(jroot) ! new norm**2
         end if
         if (jroot.eq.iroot-1) xnrm2(iroot) = xnrm
      end do
      ifree = mem_flushmark()
      end subroutine

*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      function should_discard_vector_h(trf_nrm, nopt, nroots,iroot,
     &     xnrm2,me_trv,opti_info)
*----------------------------------------------------------------------*
      implicit none

      logical::
     &     should_discard_vector_h
      real(8),intent(in)::
     &     trf_nrm,xnrm2(nroots)

      integer,intent(in)::
     &     nopt,nroots,iroot

      type(me_list_array)::
     &     me_trv(nopt)
      type(optimize_info), intent(in)::
     &     opti_info

      if (trf_nrm.lt.1d-12) then
        if (iprlvl.ge.5) write(lulog,*)
     &       'Discarding guess vector with wrong spin symmetry.'
        should_discard_vector_h =.true.
      else
        if(check_last_guess_h(nopt,iopt,nroots,iroot,
     &       (trf_nrm**2),
     &       xnrm2, me_trv,opti_info) ) then
          should_discard_vector_h =.false.
        else
          if (iprlvl.ge.5) write(lulog,*)
     &         'Discarding twin guess vector.'
          should_discard_vector_h =.true.
        end if

      end if
      return
      end function
*----------------------------------------------------------------------*
!>    check the last guess for twinning
*----------------------------------------------------------------------*
      function check_last_guess_h( nopt, iopt,nroots,iroot,guess_norm2,
     &     xnrm2, me_trv,opti_info)
*----------------------------------------------------------------------*
      implicit none
      logical ::
     &     check_last_guess_h

      integer,intent(in)::
     &     nopt, iopt,
     &     nroots,iroot
      real(8),intent(in)::
     &     guess_norm2,xnrm2(nroots)

      type(me_list_array)::
     &     me_trv(nopt)

      type(optimize_info), intent(in)::
     &     opti_info


      integer::
     &     ifree
      real(8)::
     &     xover
      real(8),pointer::
     &     xbuf1(:), xbuf2(:)
      real(8),external::
     &     da_ddot

      if (iroot .eq.1) then
         check_last_guess_h = .true.
         return
      end if
      ifree = mem_setmark('init_guess.check_guess')
      ifree = mem_alloc_real(xbuf1,opti_info%nwfpar(iopt),
     &     'xbuf1')
      ifree = mem_alloc_real(xbuf2,opti_info%nwfpar(iopt),
     &     'xbuf2')
      xover = da_ddot(me_trv(iopt)%mel%fhand,iroot-1,
     &     me_trv(iopt)%mel%fhand,iroot,
     &     opti_info%nwfpar(iopt),
     &     xbuf1, xbuf2,
     &     opti_info%nwfpar(iopt))
      ifree = mem_flushmark()
      if (abs(
     &     xover**2/xnrm2(iroot-1)-guess_norm2).lt.1d-6) then
         check_last_guess_h = .false.
      else
         check_last_guess_h = .true.
      end if
      return
      end function


*----------------------------------------------------------------------*
!>    projects out spin parts
*----------------------------------------------------------------------*
      subroutine spin_proj_h( typ_prc,len_list,me_root,me_special,
     &     nspecial, nextra, idxspc,
     &     fl_spc, xnrm,
     &     opti_info, orb_info, op_info,str_info, strmap_info)
*----------------------------------------------------------------------*
      implicit none
      integer,intent(in)::
     &     typ_prc,len_list, nspecial,nextra,idxspc
      real(8),intent(inout)::
     &     xnrm

      type(me_list)::
     &     me_root
      type(me_list_array)::
     &     me_special(nspecial+nextra)
      type(formula_item)::
     &     fl_spc(*)
      type(optimize_info)::
     &     opti_info
      type(orbinf)::
     &     orb_info
      type(operator_info)::
     &     op_info
      type(strinf)::
     &     str_info
      type(strmapinf)::
     &     strmap_info

      real(8),pointer::
     &     xbuf1(:),xbuf2(:)

      integer::
     &     ifree, len_buf
      len_buf = len_list
      ifree = mem_setmark('solve_evp2.spin_proj')
      ifree = mem_alloc_real(xbuf1,len_buf,'xbuf1')
      ifree = mem_alloc_real(xbuf2,len_buf,'xbuf2')
      select case(typ_prc)
      case(optinf_prc_spinp)
         call spin_project(me_root,me_special(idxspc)%mel,
     &        fl_spc(1),len_buf,
     &        xbuf1,xbuf2,.true.,xnrm,
     &        opti_info,orb_info,
     &        op_info,str_info,strmap_info)
      case(optinf_prc_spinrefp)
         call spin_project(me_root,me_special(idxspc)%mel,
     &        fl_spc(2),len_buf,
     &        xbuf1,xbuf2,.true.,xnrm,
     &        opti_info,orb_info,
     &        op_info,str_info,strmap_info)
         call evaluate2(fl_spc(1),.false.,.false.,
     &        op_info,str_info,strmap_info,orb_info,
     &        xnrm,.true.)
      case(optinf_prc_prj)
         call evaluate2(fl_spc(1),.false.,.false.,
     &        op_info,str_info,strmap_info,orb_info,
     &        xnrm,.true.)
      end select
      ifree = mem_flushmark()
      return
      end subroutine
*----------------------------------------------------------------------*
!!     helper routine to transform the trialvector back into non-orthogonal space
!
*----------------------------------------------------------------------*
      subroutine transform_back_h(
     &     flist, depend,
     &     trafo,
     &     me_special,
     &     me_vort,me_trv,
     &     nopt, maxvec,
     &     nspecial,
     &     me_tgt,
     &     op_info,str_info,
     &     strmap_info,orb_info,
     &     opti_info)
*----------------------------------------------------------------------*
      implicit none
      integer, intent(in)::
     &     nspecial,nopt, maxvec
      logical, intent(in)::
     &     trafo(nopt)
      type(formula_item),intent(in)::
     &     flist
      type(dependency_info),intent(in)::
     &     depend
      type(me_list_array)::
     &     me_special(nspecial),
     &     me_vort(nopt),me_trv(nopt),
     &     me_tgt(nopt)

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(optimize_info)::
     &     opti_info

      integer::
     &     iroot,iopt
      real(8)::
     &     xdummy

      do iroot=1,maxvec
        do iopt=1,nopt
          call switch_mel_record(me_vort(iopt)%mel,iroot)
          call switch_mel_record(me_trv(iopt)%mel,iroot)
          if (trafo(iopt) ) then
            call transform_back_wrap(flist,depend,
     &           me_special,me_vort(iopt)%mel,me_trv(iopt)%mel, !vort -> opt !
     &           xdummy,
     &           iopt, nspecial,
     &           me_tgt(iopt)%mel,
     &           op_info, str_info, strmap_info,
     &           orb_info, opti_info)
!     else me_vort => me_trv => me_opt
          end if
        end do
      end do
      end subroutine
*----------------------------------------------------------------------*
!!    return the maximum record of ME-list
*----------------------------------------------------------------------*
      pure function mel_get_maxrec_h(mel)
      implicit none
      integer:: mel_get_maxrec_h
      type(me_list),intent(in)::mel
      mel_get_maxrec_h=mel%fhand%active_records(2)
      end function
*----------------------------------------------------------------------*

      end

