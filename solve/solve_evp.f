*----------------------------------------------------------------------*
      subroutine solve_evp(mode_str,
     &     nopt,nroots,targ_root,label_opt,label_prc,label_op_mvp,
     &     label_op_met,label_form,
     &     label_special,nspecial,label_spcfrm,nspcfrm,thr_suggest,
     &     op_info,form_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     solve eigenvalue problem  Mx = lambda x
*
*     the formula with label "label_form" describes how to calculate 
*     the matrix trial-vector products and the r.h.s.
*
*     nopt                  number of x operators to be solved for
*                           in case of coupled equations
*     nroots                number of roots per x operator
*     targ_root             the target root per x operator
*
*     label_opt(nopt)       label of solution vectors
*     label_prc(nopt)       label of preconditioners
*     label_op_mvp(nopt)    label operators describing Mx-products
*     label_op_met(nopt)    label operators describing Sx-products
*                           if S is unity, pass label of operator
*                           associated with ME-list label_opt
*
*     the latter two are used to initilize temporary ME-lists
*
*     op_info:   operator/ME-list definitions
*     form_info: formula definitions
*     str_info: string information (to be passed to subroutines)
*     strmap_info: string mappings (to be passed to subroutines)
*     orb_info: orbital space information (to be passed)
*
*     thr_suggest: allows weaker convergence threshold
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
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_file_array.h'
      include 'mdef_formula_info.h'
      include 'def_dependency_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nopt, nroots, nspecial, nspcfrm, targ_root
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

      logical ::
     &     conv, use_s_t, use_s(nopt), trafo, init(nopt), home_in
      character(len_opname) ::
     &     label
      integer ::
     &     iter, iprint, task, ifree, iopt, jopt, nintm, irequest,
     &     nrequest, nvectors, iroot, idx, ierr, idxmel, nout,
     &     jdx
      real(8) ::
     &     xresmax, xdum, xnrm,
     &     xeig(nroots,2), xresnrm(nroots*nopt)
      type(me_list_array), pointer ::
     &     me_opt(:), me_dia(:), me_trv(:), me_mvp(:), me_met(:),
     &     me_special(:), me_scr(:), me_home(:)
      type(file_array), pointer ::
     &     ffdia(:), ff_trv(:),
     &     ffopt(:), ff_mvp(:), ff_met(:), ffspecial(:), ff_scr(:),
     &     ffhome(:)
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

      integer, pointer ::
     &     irecmvp(:), irectrv(:), irecmet(:)
      real(8), pointer ::
     &      xret(:), xbuf1(:), xbuf2(:), xoverlap(:)

      character ::
     &     fname*256

      logical, external ::
     &     file_exists
      integer, external ::
     &     idx_formlist, idx_mel_list, idx_xret
      real(8), external ::
     &     fndmnx, da_ddot

      ifree = mem_setmark('solve_evp')

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'entered solve_evp')
        write(luout,*) 'nopt   = ',nopt
        write(luout,*) 'nroots = ',nroots
        write(luout,*) 'targ_root = ',targ_root
      end if

      idx = idx_formlist(label_form,form_info)
      if (idx.le.0)
     &     call quit(1,'solve_evp',
     &     'did not find formula '//trim(label_form))
      form_mvp => form_info%form_arr(idx)%form

      allocate(me_opt(nopt),me_dia(nopt),me_trv(nopt),me_mvp(nopt),
     &         me_met(nopt),me_special(nspecial),me_scr(nopt))
      allocate(ffopt(nopt),ffdia(nopt),
     &     ff_trv(nopt),ff_mvp(nopt),ff_met(nopt),ffspecial(nspecial),
     &     ff_scr(nopt))
      do iopt = 1, nopt
        ! pointer array for operators:
        ierr = 1
        jopt = iopt
        idxmel = idx_mel_list(label_opt(iopt),op_info)
        if (idxmel.le.0) exit
        ierr = 2
        me_opt(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        ffopt(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand
        if (.not.associated(ffopt(iopt)%fhand)) exit
        ierr = 3
        jopt = iopt
        idxmel = idx_mel_list(label_prc(iopt),op_info)
        if (idxmel.le.0) exit
        me_dia(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        ierr = 4
        ffdia(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand
        if (.not.associated(ffdia(iopt)%fhand)) exit
        ierr = 0
      end do

      ! special lists needed?
      if (ierr.eq.0) then
        do idx = 1, nspecial
          jopt = idx
          idxmel = idx_mel_list(label_special(idx),op_info)
          ierr = 5
          if (idxmel.le.0) exit
          me_special(idx)%mel  => op_info%mel_arr(idxmel)%mel
          ffspecial(idx)%fhand => op_info%mel_arr(idxmel)%mel%fhand
          ierr = 6
          if (.not.associated(ffspecial(idx)%fhand)) exit
          ierr = 0
        end do
      end if

      ! error handling
      if (ierr.gt.0) then
        if (ierr.eq.1.or.ierr.eq.2) label = label_opt(jopt)
        if (ierr.eq.3.or.ierr.eq.4) label = label_prc(jopt)
        if (ierr.eq.5.or.ierr.eq.6) label = label_special(jopt)
        if (mod(ierr,2).eq.1)
     &       call quit(1,'solve_evp',
     &       'did not find list '//trim(label))
        if (mod(ierr,2).eq.0)
     &       call quit(1,'solve_evp',
     &       'no file associated with list '//trim(label))
      end if

      ! special formulae
      do jdx = 1, nspcfrm
        idx = idx_formlist(label_spcfrm(jdx),form_info)
        if (idx.le.0)
     &       call quit(1,'solve_evp',
     &       'did not find formula '//trim(label_spcfrm(jdx)))
        ! read formula
        call read_form_list(form_info%form_arr(idx)%form%fhand,
     &                      fl_spc(jdx),.true.)
      end do

      call set_opti_info(opti_info,3,nopt,nroots,me_opt,mode_str)

      nvectors = opti_info%maxsbsp
      use_s_t = .false.

      do iopt = 1, nopt
        ! weaker convergence threshold requested?
        opti_info%thrgrd(iopt)=max(opti_info%thrgrd(iopt),thr_suggest)

        ! get a ME-list for scratch trial-vectors
        ! in case of ab-sym braking trafo, get sym props from special list
        if (opti_info%typ_prc(iopt).eq.optinf_prc_traf
     &      .and.nspecial.eq.3) then
          me_pnt => me_special(1)%mel
        else
          me_pnt => me_opt(iopt)%mel
        end if
        write(fname,'("scr_",i3.3)') iopt
        call define_me_list(fname,me_opt(iopt)%mel%op%name,
     &       me_pnt%absym,me_pnt%casym,
     &       me_pnt%gamt,me_pnt%s2,
     &       me_pnt%mst,.false.,
     &       -1,1,nvectors,0,0,0,
     &       op_info,orb_info,str_info,strmap_info)
        idxmel = idx_mel_list(fname,op_info)
        me_scr(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        ff_scr(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand

        ! get a ME-list for trial-vectors
        write(fname,'("trv_",i3.3)') iopt
        call define_me_list(fname,me_opt(iopt)%mel%op%name,
     &       me_opt(iopt)%mel%absym,me_opt(iopt)%mel%casym,
     &       me_opt(iopt)%mel%gamt,me_opt(iopt)%mel%s2,
     &       me_opt(iopt)%mel%mst,.false.,
     &       -1,1,nvectors,0,0,0,
     &       op_info,orb_info,str_info,strmap_info)
        idxmel = idx_mel_list(fname,op_info)
        me_trv(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        ff_trv(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand

        ! get a ME list for matrix-vector products
        ! (have same symmtry properties as result!)
        write(fname,'("mvp_",i3.3)') iopt
        call define_me_list(fname,label_op_mvp(iopt),
     &       me_opt(iopt)%mel%absym,me_opt(iopt)%mel%casym,
     &       me_opt(iopt)%mel%gamt,me_opt(iopt)%mel%s2,
     &       me_opt(iopt)%mel%mst,.false.,
     &       -1,1,nvectors,0,0,0,
     &       op_info,orb_info,str_info,strmap_info)
        idxmel = idx_mel_list(fname,op_info)
        me_mvp(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        ff_mvp(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand

        ! use of metric requested?
        use_s(iopt) = trim(label_op_met(iopt)).ne.
     &       trim(me_opt(iopt)%mel%op%name)
        if (use_s(iopt)) then
          use_s_t = .true.
          ! get a ME list for metric-times-vector products
          ! (have same symmtry properties as result!)
          write(fname,'("svp_",i3.3)') iopt
          call define_me_list(fname,label_op_met(iopt),
     &         me_opt(iopt)%mel%absym,me_opt(iopt)%mel%casym,
     &         me_opt(iopt)%mel%gamt,me_opt(iopt)%mel%s2,
     &         me_opt(iopt)%mel%mst,.false.,
     &         -1,1,nvectors,0,0,0,
     &         op_info,orb_info,str_info,strmap_info)
          idxmel = idx_mel_list(fname,op_info)
          me_met(iopt)%mel   => op_info%mel_arr(idxmel)%mel
          ff_met(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand          
        end if

      end do

      ! write information to opti_info about signs which occur
      ! in trv*mvp or trv*met  multiplication
      ! relevant when trv is njoined=1 op. but mvp (met) are njoined=2 op's
      call set_opti_info_signs(opti_info,3,nopt,
     &                         me_trv,me_mvp,me_met,me_met,use_s)

      ! read formula
      call read_form_list(form_mvp%fhand,fl_mvp,.true.)

      ! set dependency info for submitted formula list
      call set_formula_dependencies(depend,fl_mvp,op_info)

      ! number of info values returned on xret
      nout = depend%ntargets
      allocate(xret(nout))

      ! records with trial vectors and Mv-products, needed in evp_control:
      ifree = mem_alloc_int(irectrv,nroots,'rectrv')
      ifree = mem_alloc_int(irecmvp,nroots,'recmvp')
      ifree = mem_alloc_int(irecmet,nroots,'recmet')

      init = .true.
      home_in = .false.
      do iopt = 1, nopt
        ! open result vector file(s)
        ! if file already open, use as initial guess (if requested)!
        if (ffopt(iopt)%fhand%unit.gt.0.and.opti_info%resume) then
          write(luout,'(a,i4,a)')
     &          'Using last vector as initial guess (iopt =',iopt,')'
          init(iopt) = .false.
        end if
        if (ffopt(iopt)%fhand%unit.gt.0.and.nroots.gt.1) then
          ! copy this root so that we may home in on it later
          if (nopt.ne.1) call quit(1,'solve_evp',
     &         'homing in available only for one opt. vector yet')
          home_in = .true.
          allocate(me_home(1),ffhome(1))
          if (opti_info%typ_prc(iopt).eq.optinf_prc_traf
     &        .and.nspecial.eq.3) then
            me_pnt => me_special(1)%mel
          else
            me_pnt => me_opt(iopt)%mel
          end if
          call define_me_list('home',me_opt(iopt)%mel%op%name,
     &         me_pnt%absym,me_pnt%casym,
     &         me_pnt%gamt,me_pnt%s2,
     &         me_pnt%mst,.false.,
     &         -1,1,1,0,0,0,
     &         op_info,orb_info,str_info,strmap_info)
          idxmel = idx_mel_list('home',op_info)
          me_home(1)%mel   => op_info%mel_arr(idxmel)%mel
          ffhome(1)%fhand => op_info%mel_arr(idxmel)%mel%fhand
          call assign_me_list(me_trv(iopt)%mel%label,
     &                        me_opt(iopt)%mel%op%name,op_info)
          call file_open(ffhome(1)%fhand)
          call list_copy(me_opt(iopt)%mel,me_home(1)%mel,.false.)
c dbg
c          print *,'preparing for homing in later. Saved vector:'
c          call wrt_mel_file(luout,5,
c     &         me_home(1)%mel,
c     &         1,me_trv(iopt)%mel%op%n_occ_cls,
c     &         str_info,orb_info)
c dbgend
        else if (ffopt(iopt)%fhand%unit.le.0) then
          call file_open(ffopt(iopt)%fhand)
        end if
        call file_open(ff_scr(iopt)%fhand)
        call file_open(ff_trv(iopt)%fhand)
        ! open corresponding matrix vector products ...
        call file_open(ff_mvp(iopt)%fhand)
        if (use_s(iopt))
     &       call file_open(ff_met(iopt)%fhand)
        ! ... and corresponding preconditioner(s)
        if (ffdia(iopt)%fhand%unit.le.0)
     &       call file_open(ffdia(iopt)%fhand)
      end do

      do idx = 1, nspecial
        if (ffspecial(idx)%fhand%unit.le.0)
     &       call file_open(ffspecial(idx)%fhand)
      end do

      ! get initial amplitudes
      call init_guess(nopt,init,nroots,
     &                me_opt,me_trv,me_dia,me_special,nspecial,
     &                fl_mvp,depend,fl_spc,nspcfrm,
     &                opti_info,orb_info,op_info,str_info,strmap_info)

      ! start optimization loop
      iter = 0
      task = 0
      opt_loop: do while(task.lt.8)

        call leq_evp_control
     &       ('EVP',iter,
     &       task,conv,xresnrm,xeig,
     &       use_s,
     &       nrequest,irectrv,irecmvp,irecmet,
     &       me_opt,me_scr,me_trv,me_mvp,me_met,me_dia,me_dia,
     &       me_special,nspecial,
c     &       ffopt,ff_trv,ff_mvp,ff_met,ffdia,ffdia,  ! #5 is dummy
     &       fl_mvp,depend,
     &       fl_spc,nspcfrm,
     &       opti_info,opti_stat,
     &       orb_info,op_info,str_info,strmap_info)

        if (iter.gt.1) then
          xresmax = fndmnx(xresnrm,nroots*nopt,2)
          write(luout,'("E>>",i3,24x,x,g10.4)') iter-1,xresmax
          if (iprlvl.gt.0) then
            do iroot = 1, nroots
              if (xeig(iroot,2).eq.0d0) then
                write(luout,'(" E>",3x,f24.12,x,3g10.4)')
     &               xeig(iroot,1),(xresnrm(iroot+idx*nroots),
     &                              idx = 0, nopt-1)
              else
                write(luout,
     &               '(" E>",3x,f24.12,x,g10.4," (img=",g24.12,")")')
     &               xeig(iroot,1),xresnrm(iroot),xeig(iroot,2)
              end if
            end do
          end if
        end if

        ! 4 - get residual
        if (iand(task,4).eq.4) then
          ! preliminary solution: 
          !   outside loop over requested Mv-products
          do irequest = 1, nrequest
            do iopt = 1, nopt
              call switch_mel_record(me_trv(iopt)%mel,irectrv(irequest))
              call switch_mel_record(me_mvp(iopt)%mel,irecmvp(irequest))
              if (use_s(iopt))
     &             call switch_mel_record(me_met(iopt)%mel,
     &                                                irecmet(irequest))
              
              ! enforce MS-combination symmetry of trial vectors
              ! (if requested)
c              if (me_trv(iopt)%mel%absym.ne.0)
c dbg
c        write(luout,*) 'current trial vector (before):'
c        call wrt_mel_file(luout,5,
c     &       me_trv(iopt)%mel,
c     &       1,me_trv(iopt)%mel%op%n_occ_cls,
c     &       str_info,orb_info)
c dbgend
              if (iter.gt.1.and.me_trv(iopt)%mel%absym.ne.0)
     &             call sym_ab_list(
     &             0.5d0,me_trv(iopt)%mel,me_trv(iopt)%mel,
     &             xdum,.false.,
     &             op_info,str_info,strmap_info,orb_info)

              ! here?
              call touch_file_rec(me_trv(iopt)%mel%fhand)
            end do

c dbg
c            write(luout,*) 'input for request: ',irequest
c            call wrt_mel_file(luout,5,me_trv(1)%mel,
c     &           1,me_trv(1)%mel%op%n_occ_cls,
c     &           str_info,orb_info)
c dbg

            call frm_sched(xret,fl_mvp,depend,0,0,
     &           .true.,.false.,op_info,str_info,strmap_info,orb_info)

c dbg
c            write(luout,*) 'output for request: ',irequest
c            call wrt_mel_file(luout,5,me_mvp(1)%mel,
c     &           1,me_mvp(1)%mel%op%n_occ_cls,
c     &           str_info,orb_info)
c dbg
              ! enforce MS-combination symmetry of Mv-products:
              ! (if requested)
            do iopt = 1, nopt
              if (me_mvp(iopt)%mel%absym.ne.0)
     &             call sym_ab_list(
     &             0.5d0,me_mvp(iopt)%mel,me_mvp(iopt)%mel,
     &             xdum,.false.,
     &             op_info,str_info,strmap_info,orb_info)

              ! project out spin contaminations?
              if (opti_info%typ_prc(iopt).eq.optinf_prc_spinp.or.
     &            opti_info%typ_prc(iopt).eq.optinf_prc_prj) then
                ifree = mem_setmark('solve_evp.spin_proj_res')
                ifree = mem_alloc_real(xbuf1,opti_info%nwfpar(iopt),
     &                                 'xbuf1')
                ifree = mem_alloc_real(xbuf2,opti_info%nwfpar(iopt),
     &                                 'xbuf2')
                ! assign op. with list containing the mvp vector
                call assign_me_list(me_mvp(iopt)%mel%label,
     &                              me_opt(iopt)%mel%op%name,op_info)
                if (opti_info%typ_prc(iopt).eq.optinf_prc_spinp) then
                  call spin_project(me_mvp(iopt)%mel,me_special(1)%mel,
     &                             fl_spc(1),opti_info%nwfpar(iopt),
     &                             xbuf1,xbuf2,.false.,xnrm,
     &                             opti_info,orb_info,
     &                             op_info,str_info,strmap_info)
                else
                  call evaluate2(fl_spc(1),.false.,.false.,
     &                           op_info,str_info,strmap_info,orb_info,
     &                           xnrm,.false.)
                end if
                ! reassign lists to correct ops
                call assign_me_list(me_trv(iopt)%mel%label,
     &                              me_opt(iopt)%mel%op%name,op_info)
                call assign_me_list(me_mvp(iopt)%mel%label,
     &                              label_op_mvp(iopt),op_info)
                ifree = mem_flushmark()
              end if
            end do

          end do
        end if

      end do opt_loop

      do iopt = 1, nopt

        ! remove the temporary lists
        call del_me_list(me_scr(iopt)%mel%label,op_info)
        call del_me_list(me_trv(iopt)%mel%label,op_info)
        call del_me_list(me_mvp(iopt)%mel%label,op_info)
        if (use_s(iopt))
     &       call del_me_list(me_met(iopt)%mel%label,op_info)

        ! make sure that the operator is now associated with
        ! the list containing the solution vector
        call assign_me_list(label_opt(iopt),
     &                      me_opt(iopt)%mel%op%name,op_info)


c        ! solution vector has been updated (if we had some iteration)
c        if (iter.gt.1) call touch_file_rec(me_opt(iopt)%mel%fhand)
        ! solution vector has been updated
        call touch_file_rec(me_opt(iopt)%mel%fhand)

        if (home_in) then
          ! home in on root with largest overlap with prior solution
          ifree = mem_setmark('solve_evp.home_in')
          ifree = mem_alloc_real(xoverlap,nroots,'xoverlap')
          ifree = mem_alloc_real(xbuf1,opti_info%nwfpar(iopt),'xbuf1')
          ifree = mem_alloc_real(xbuf2,opti_info%nwfpar(iopt),'xbuf2')
          xresmax = 0d0
          do iroot = 1, nroots
            xoverlap(iroot) = da_ddot(ffhome(1)%fhand,1,1,
     &                                ffopt(iopt)%fhand,iroot,1,
     &                                opti_info%nwfpar(iopt),
     &                                xbuf1,xbuf2,
     &                                opti_info%nwfpar(iopt))
            xoverlap(iroot) = abs(xoverlap(iroot))
            if (xoverlap(iroot).gt.xresmax) then
              idx = iroot
              xresmax = xoverlap(iroot)
            end if
c dbg
c            print *,'root / overlap: ',iroot,xoverlap(iroot)
c dbgend
          end do
          ifree = mem_flushmark()
          if (idx.ne.targ_root) then
            write(luout,'(a,i4,a,f8.4)') 
     &            'Homing in on root ',idx,' with overlap ',xresmax
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
        if (ffspecial(idx)%fhand%unit.gt.0)
     &       call file_close_keep(ffspecial(idx)%fhand)
      end do

      ! print results
      call write_title(luout,wst_title,
     &     'Results for '//trim(label_opt(1)))
      write(luout,'("E>>",66("="))')
      write(luout,'("E>>",2x,'//
     &     '"root     eigenvalue (real)       eigenvalue (img.)'//
     &     '  |residual|")')
      write(luout,'("E>>",66("-"))') 
      do iroot = 1, nroots
        if (xeig(iroot,2).eq.0d0) then
          write(luout,'("E>>",2x,i3,x,f22.12,20x,"---",2x,x,g10.4)')
     &         iroot,xeig(iroot,1),xresnrm(iroot)
        else
          write(luout,
     &         '("E>>",3x,i2,x,f22.12,x,g24.12,x,g10.4)')
     &         iroot,xeig(iroot,1:2),xresnrm(iroot)
        end if
c dbg
c         do iopt=1,nopt
c           call switch_mel_record(me_opt(iopt)%mel,iroot)
c           call wrt_mel_file(luout,5,me_opt(iopt)%mel,
c    &             1,me_opt(iopt)%mel%op%n_occ_cls,
c    &             str_info,orb_info)
c         enddo
c dbg     
      end do
      write(luout,'("E>>",66("="))') 

      ! switch to target root if possible
!      ! (we assume that nroots has been chosen for this reason,
!      !  otherwise a new keyword must be set up for this purpose)
!      if (nroots.ge.me_opt(1)%mel%fhand%active_records(1).and.
!     &    nroots.le.me_opt(1)%mel%fhand%active_records(2))
!     &           call switch_mel_record(me_opt(1)%mel,nroots)
!
      if (targ_root.ge.me_opt(1)%mel%fhand%active_records(1).and.
     &    targ_root.le.me_opt(1)%mel%fhand%active_records(2))
     &           call switch_mel_record(me_opt(1)%mel,targ_root)

      call clean_formula_dependencies(depend)

      ! note that only the pointer array ffopt (but not the entries)
      ! is deallocated:
      deallocate(me_opt,me_dia,me_trv,me_mvp,me_met,me_special,me_scr)
      deallocate(ff_trv,ff_mvp,ffdia,ffopt,ff_met,xret,ffspecial,ff_scr)
      call dealloc_formula_list(fl_mvp)
      do jdx = 1, nspcfrm
        call dealloc_formula_list(fl_spc(jdx))
      end do

      ifree = mem_flushmark()

      return
      end


