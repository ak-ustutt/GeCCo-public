*----------------------------------------------------------------------*
      subroutine solve_leq(mode_str,
     &     nopt,nroots,label_opt,label_prc,label_op_mvp,label_op_met,
     &     label_op_rhs,xrhsnorm,
     &     label_form,
     &     label_special,nspecial,label_spcfrm,nspcfrm,thr_suggest,
     &     op_info,form_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     solve linear equations  Mx = -g
*
*     the formula with label "label_form" describes how to calculate 
*     the matrix trial-vector products and the r.h.s.
*
*     nopt                  number of x operators to be solved for
*                           in case of coupled equations
*     nroots                number of roots per x operator
*
*     label_opt(nopt)       label of solution vectors
*     label_prc(nopt)       label of preconditioners
*     label_op_mvp(nopt)    label operators describing Mx-products
*     label_op_rhs(nopt)    label of r.h.s. operators (g)
*     label_op_met(nopt)    label operators describing Sx-products
*                           if S is unity, pass label of operator
*                           associated with ME-list label_opt
*     xrhsnorm              return the norm of the RHS (just in case)
*
*     the latter two are used to initilize temporary ME-lists
*
*     op_info:   operator/ME-list definitions
*     form_info: formula definitions
*     str_info: string information (to be passed to subroutines)
*     strmap_info: string mappings (to be passed to subroutines)
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
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_file_array.h'
      include 'mdef_formula_info.h'
      include 'def_dependency_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nopt, nroots, nspecial, nspcfrm
      character(*), intent(in) ::
     &     mode_str,
     &     label_opt(nopt),
     &     label_prc(nopt),
     &     label_op_mvp(nopt),
     &     label_op_met(nopt),
     &     label_special(nspecial),
     &     label_op_rhs(nopt),
     &     label_spcfrm(nspcfrm),
     &     label_form
      real(8), intent(in) ::
     &     thr_suggest
      real(8), intent(out) ::
     &     xrhsnorm
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
     &     conv, use_s_t, use_s(nopt)
      character(len_opname) ::
     &     label
      integer ::
     &     iter, iprint, task, ifree, iopt, jopt, nintm, irequest,
     &     nrequest, nvectors, iroot, idx, ierr, idxmel, nout, idxrhs,
     &     nselect, jdx
      real(8) ::
     &     energy, xresnrm(nroots,nopt), xdum, xresmax
      type(me_list_array), pointer ::
     &     me_opt(:), me_trv(:), me_mvp(:), me_rhs(:), me_dia(:),
     &     me_met(:), me_special(:), me_scr(:), me_ext(:)
      type(file_array), pointer ::
     &     ffdia(:), ff_rhs(:), ff_trv(:),ff_ext(:),
     &     ffopt(:), ff_mvp(:), ff_met(:), ffspecial(:), ff_scr(:)
      type(me_list), pointer ::
     &     me_pnt
      type(dependency_info) ::
     &     depend
      type(optimize_info) ::
     &     opti_info
      type(optimize_status) ::
     &     opti_stat
      type(formula), pointer ::
     &     form_rhs_mvp
      type(formula_item) ::
     &     fl_rhs_mvp, fl_spc(nspcfrm)

      integer, pointer ::
     &     irecmvp(:), irectrv(:), irecmet(:), idxselect(:)
      real(8), pointer ::
     &      xret(:), xbuf1(:)

      character ::
     &     fname*256

      integer, external ::
     &     idx_formlist, idx_mel_list, idx_oplist2, idx_xret
      real(8), external ::
     &     fndmnx

      ifree = mem_setmark('solve_leq')

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'entered solve_leq')
        write(lulog,*) 'nopt   = ',nopt
        write(lulog,*) 'nroots = ',nroots
      end if

c      if (nopt.gt.1)
c     &     call quit(1,'solve_leq','did not yet consider coupled LEQs')

      idx = idx_formlist(label_form,form_info)
      if (idx.le.0)
     &     call quit(1,'solve_leq',
     &     'did not find formula '//trim(label_form))
      form_rhs_mvp => form_info%form_arr(idx)%form

      use_s(1:nopt) = .false.

      allocate(me_opt(nopt),me_rhs(nopt),me_trv(nopt),me_mvp(nopt),
     &     me_dia(nopt),me_met(nopt),me_special(nspecial),
     &     me_scr(nopt),me_ext(nopt))
      allocate(ffopt(nopt),ffdia(nopt),
     &     ff_trv(nopt),ff_mvp(nopt),ff_rhs(nopt),ff_ext(nopt),
     &     ff_met(nopt),ffspecial(nspecial),ff_scr(nopt))
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
        ierr = 4
        me_dia(iopt)%mel  => op_info%mel_arr(idxmel)%mel
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
     &       call quit(1,'solve_leq',
     &       'did not find list '//trim(label))
        if (mod(ierr,2).eq.0)
     &       call quit(1,'solve_leq',
     &       'no file associated to list '//trim(label))
      end if

      ! special formulae
      do jdx = 1, nspcfrm
        idx = idx_formlist(label_spcfrm(jdx),form_info)
        if (idx.le.0)
     &       call quit(1,'solve_leq',
     &       'did not find formula '//trim(label_spcfrm(jdx)))
        ! read formula
c dbg
        print *,'reading special form: ',trim(label_spcfrm(jdx))
c dbg
        call read_form_list(form_info%form_arr(idx)%form%fhand,
     &                      fl_spc(jdx),.true.)
      end do      

      call set_opti_info(opti_info,2,nopt,nroots,me_opt,mode_str)

      nvectors = opti_info%maxsbsp
      use_s_t = .false.

      do iopt = 1, nopt
        ! weaker convergence threshold requested?
        opti_info%thrgrd(iopt)=max(opti_info%thrgrd(iopt),thr_suggest)

        ! get a ME-list for scratch vectors
        ! in case of ab-sym braking trafo, get sym props from special list
        if (opti_info%typ_prc(iopt).eq.optinf_prc_traf
     &      .and.nspecial.eq.4) then
          me_pnt => me_special(2)%mel
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

        ! Here is a new ME-list that will be fed in the 
        ! routine 'optc_minspace'. Previously ff_scr was
        ! used there but the use was erroneous
        write(fname,'("ext_",i3.3)') iopt
        call define_me_list(fname,me_opt(iopt)%mel%op%name,
     &       me_pnt%absym,me_pnt%casym,
     &       me_pnt%gamt,me_pnt%s2,
     &       me_pnt%mst,.false.,
     &       -1,1,nvectors,0,0,0,
     &       op_info,orb_info,str_info,strmap_info)
        idxmel = idx_mel_list(fname,op_info)
        me_ext(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        ff_ext(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand

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

        ! get a ME list for RHS
        write(fname,'("rhs_",i3.3)') iopt
        call define_me_list(fname,label_op_rhs(iopt),
     &       me_opt(iopt)%mel%absym,me_opt(iopt)%mel%casym,
     &       me_opt(iopt)%mel%gamt,me_opt(iopt)%mel%s2,
     &       me_opt(iopt)%mel%mst,.false.,
     &       -1,1,nvectors,0,0,0,
     &       op_info,orb_info,str_info,strmap_info)
        idxmel = idx_mel_list(fname,op_info)
        me_rhs(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        ff_rhs(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand

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
      ! in trv*mvp or trv*met or trv*rhs multiplication
      ! relevant when trv is njoined=1 op. but mvp (met) are njoined=2 op's
      call set_opti_info_signs(opti_info,2,nopt,
     &                         me_trv,me_mvp,me_met,me_rhs,use_s)

      ! read formula
      call read_form_list(form_rhs_mvp%fhand,fl_rhs_mvp,.true.)

      ! set dependency info for submitted formula list
      call set_formula_dependencies(depend,fl_rhs_mvp,op_info)

      ! number of info values returned on xret
      nout = depend%ntargets
      allocate(xret(nout),idxselect(nout))

      ! records with trial vectors and Mv-products, needed in leq_control:
      ifree = mem_alloc_int(irectrv,nroots,'rectrv')
      ifree = mem_alloc_int(irecmvp,nroots,'recmvp')
      ifree = mem_alloc_int(irecmet,nroots,'recmet')

      do iopt = 1, nopt
        ! open result vector file(s)
        call file_open(ffopt(iopt)%fhand)
        call file_open(ff_scr(iopt)%fhand)
        call file_open(ff_ext(iopt)%fhand)
        call file_open(ff_trv(iopt)%fhand)
        ! open corresponding matrix vector products ...
        call file_open(ff_mvp(iopt)%fhand)
        if (use_s(iopt))
     &       call file_open(ff_met(iopt)%fhand)
        ! right hand sides ...
        call file_open(ff_rhs(iopt)%fhand)
        ! ... and corresponding preconditioner(s)
        if (ffdia(iopt)%fhand%unit.le.0)
     &       call file_open(ffdia(iopt)%fhand)
      end do

      do idx = 1, nspecial
        if (ffspecial(idx)%fhand%unit.le.0)
     &       call file_open(ffspecial(idx)%fhand)
      end do

      ! get initial amplitudes
      do iopt = 1, nopt
        ! get idx of RHS operator
        nselect = 0
        call select_formula_target(idxselect,nselect,
     &       me_rhs(iopt)%mel%label,depend,op_info)

        do iroot = 1, nroots
          
          call switch_mel_record(me_rhs(iopt)%mel,iroot)

          call frm_sched(xret,fl_rhs_mvp,depend,idxselect,nselect,
     &         .true.,.false.,op_info,str_info,strmap_info,orb_info)

          call touch_file_rec(me_rhs(iopt)%mel%fhand)

          ! store norm of RHS
          xresnrm(iroot,iopt) = xret(idxselect(1))

          ! apply sign-fix (if needed)
          ifree = mem_setmark('solve_leq.fix_sign')
          ifree = mem_alloc_real(xbuf1,opti_info%nwfpar(iopt),'xbuf1')
          call optc_fix_signs2(me_rhs(iopt)%mel%fhand,
     &                        iroot,opti_info,iopt,
     &                        opti_info%nwfpar(iopt),xbuf1)
          ifree = mem_flushmark()
        end do
      end do

      if (ntest.ge.1000.and.nroots.eq.1) then
        do iopt = 1, nopt
          write(lulog,*) 'dump of '//trim(me_rhs(iopt)%mel%label)
          write(lulog,*) 'iopt = ',iopt
          call analyze_list_core(
     &         me_rhs(iopt),me_rhs(iopt),1,1,'norm',
     &         orb_info,str_info)
          call wrt_mel_file(lulog,5,
     &         me_rhs(iopt)%mel,
     &         1,me_rhs(iopt)%mel%op%n_occ_cls,
     &         str_info,orb_info)
        end do
      end if

      ! start optimization loop
      iter = 0
      task = 0
      opt_loop: do while(task.lt.8)

        call leq_evp_control
     &       ('LEQ',iter,
     &       task,conv,xresnrm,xdum,
     &       use_s,
     &       nrequest,irectrv,irecmvp,irecmet, 
     &       me_opt,me_scr,me_trv,me_mvp,me_met,me_rhs,me_dia,
     &       me_ext,
     &       me_special,nspecial,0,1,
c     &       ffopt,ff_trv,ff_mvp,ff_mvp,ff_rhs,ffdia, ! dto.
     &       fl_rhs_mvp,depend,
     &       fl_spc,nspcfrm,
     &       opti_info,opti_stat,
     &       orb_info,op_info,str_info,strmap_info)

        do iopt = 1, nopt
          xresmax = fndmnx(xresnrm(1:nroots,iopt),nroots,2)
          if (conv) then
            write(lulog,'("L>> conv.",21x,x,g10.4)') xresmax
            if (lulog.ne.luout)
     &         write(luout,'("    conv.",21x,x,g10.4)') xresmax
          else if (iter.eq.1) then
            write(lulog,'("L>> |rhs|",21x,x,g10.4)') xresmax
            if (lulog.ne.luout)   
     &        write(luout,'("    |rhs|",21x,x,g10.4)') xresmax
            xrhsnorm = xresmax
          else
            write(lulog,'("L>>",i3,24x,x,g10.4)')iter-1,xresmax
            if (lulog.ne.luout) 
     &         write(luout,'("   ",i3,24x,x,g10.4)')iter-1,xresmax
          end if
        end do

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
            end do

            if (ntest.ge.1000) then
              do iopt = 1, nopt
                write(lulog,*) 'dump of '//trim(me_trv(iopt)%mel%label)
                write(lulog,*) 'iopt = ',iopt
                call analyze_list_core(
     &            me_trv(iopt),me_trv(iopt),1,1,'norm',
     &                orb_info,str_info)
                call wrt_mel_file(lulog,5,
     &               me_trv(iopt)%mel,
     &               1,me_trv(iopt)%mel%op%n_occ_cls,
     &               str_info,orb_info)
              end do
            end if

            call frm_sched(xret,fl_rhs_mvp,depend,0,0,
     &           .true.,.false.,op_info,str_info,strmap_info,orb_info)

            do iopt = 1, nopt
              call touch_file_rec(me_trv(iopt)%mel%fhand)
            end do

            ! apply sign-fix (if needed)
            do iopt = 1, nopt
c             write(lulog,*) 'Fixing signs of residual+metric,iopt=',iopt
              ifree = mem_setmark('solve_leq.fix_sign')
              ifree = mem_alloc_real(xbuf1,opti_info%nwfpar(iopt),
     &                                         'xbuf1')
              call optc_fix_signs2(me_mvp(iopt)%mel%fhand,
     &                            irecmvp(irequest),
     &                            opti_info,iopt,
     &                            opti_info%nwfpar(iopt),xbuf1)
              if (use_s(iopt))
     &           call optc_fix_signs2(me_met(iopt)%mel%fhand,
     &                            irecmet(irequest),
     &                            opti_info,iopt,
     &                            opti_info%nwfpar(iopt),xbuf1)
              ifree = mem_flushmark()
            end do

            if (ntest.ge.1000) then
              do iopt = 1, nopt
                write(lulog,*) 'dump of '//
     &               trim(me_mvp(iopt)%mel%label)
               call analyze_list_core(
     &            me_mvp(iopt),me_mvp(iopt),1,1,'norm',
     &                orb_info,str_info)
                call wrt_mel_file(lulog,5,
     &               me_mvp(iopt)%mel,
     &               1,me_mvp(iopt)%mel%op%n_occ_cls,
     &               str_info,orb_info)
              end do
            end if

          end do
        end if

      end do opt_loop

      do iopt = 1, nopt

        ! remove the temporary lists
        call del_me_list(me_scr(iopt)%mel%label,op_info)
        call del_me_list(me_ext(iopt)%mel%label,op_info)
        call del_me_list(me_trv(iopt)%mel%label,op_info)
        call del_me_list(me_mvp(iopt)%mel%label,op_info)
        call del_me_list(me_rhs(iopt)%mel%label,op_info)
        if (use_s(iopt))
     &       call del_me_list(me_met(iopt)%mel%label,op_info)

        ! make sure that the operator is now associated with
        ! the list containing the solution vector
        call assign_me_list(label_opt(iopt),
     &                      me_opt(iopt)%mel%op%name,op_info)

        ! close the output file
        call file_close_keep(ffopt(iopt)%fhand)

      end do

      do idx = 1, nspecial
        if (ffspecial(idx)%fhand%unit.gt.0)
     &       call file_close_keep(ffspecial(idx)%fhand)
      end do

      if (ntest.ge.1000) then
        do iopt = 1, nopt
          write(lulog,*) 'dump of final '//trim(me_opt(iopt)%mel%label)
          write(lulog,*) 'iopt = ',iopt
          call analyze_list_core(
     &            me_opt(iopt),me_opt(iopt),1,1,'norm',
     &            orb_info,str_info)
          call wrt_mel_file(lulog,5,
     &         me_opt(iopt)%mel,
     &         1,me_opt(iopt)%mel%op%n_occ_cls,
     &         str_info,orb_info)
        end do
      end if

      call clean_formula_dependencies(depend)

      ! note that only the pointer array ffopt (but not the entries)
      ! is deallocated:
      deallocate(me_opt,me_trv,me_rhs,me_mvp,me_dia,me_met,me_special,
     &           me_scr,me_ext)
      deallocate(ff_trv,ff_rhs,ff_mvp,ffdia,ffopt,ff_met,ffspecial,
     &     xret,idxselect,ff_scr,ff_ext)
      call dealloc_formula_list(fl_rhs_mvp)
      do jdx = 1, nspcfrm
        call dealloc_formula_list(fl_spc(jdx))
      end do

      ifree = mem_flushmark()

      return
      end


