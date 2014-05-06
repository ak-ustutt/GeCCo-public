*----------------------------------------------------------------------*
      subroutine solve_nleq(mode_str,
     &     nopt,label_opt,label_res,label_prc,label_en,
     &     label_form,
     &     label_special,nspecial,           !<- eg. for R12
     &     label_spcfrm,nspcfrm,             !<- eg. for MRCC
     &     op_info,form_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     solve non-linear equations
*
*     the formula with label "label_form" describes how to calculate 
*     energy and residual
*
*     nopt               number of operators to be simultaneously optimized
*     label_opt(1..nop_opt) labels of those operators
*     label_opt(1..nop_opt) labels of preconditioners
*     label_res(nop_opt+1,..)   residuals
*     label_en                  energy
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
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_file_array.h'
      include 'mdef_formula_info.h'
      include 'def_dependency_info.h'
      include 'ifc_memman.h'
c dbg
      include 'ifc_input.h'
c dbgend
      include 'routes.h'

      integer, parameter ::
     &     ntest = 000

      integer, intent(in) ::
     &     nopt, nspecial, nspcfrm
      character(*), intent(in) ::
     &     mode_str,
     &     label_opt(nopt),
     &     label_res(nopt),
     &     label_prc(nopt),
     &     label_special(nspecial),
     &     label_spcfrm(nspcfrm),
     &     label_en,
     &     label_form
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
     &     conv, restart, traf
      character(len_opname) ::
     &     label, dia_label
      integer ::
     &     imacit, imicit, imicit_tot, iprint, task, ifree, iopt, jopt,
     &     idx, idxmel, ierr, nout, idx_en_xret,
     &     ndx, idx_res_xret(nopt), jdx, it_print
      real(8) ::
     &     energy, xresnrm(nopt), xdum, thr_suggest
      real(8), pointer ::
     &     xret(:)
      type(dependency_info) ::
     &     depend
      type(filinf), target ::
     &     ff_unit_dummy
      type(me_list_array), pointer ::
     &     me_opt(:), me_grd(:), me_dia(:), me_special(:),
     &     me_trv(:), me_h_trv(:)   ! not yet needed
      type(file_array), pointer ::
     &     ffopt(:), ffgrd(:), ffdia(:), ffspecial(:),
     &     ff_trv(:), ff_h_trv(:)   ! not yet needed
      type(optimize_info) ::
     &     opti_info
      type(optimize_status) ::
     &     opti_stat
      type(formula), pointer ::
     &     form_en_res
      type(formula_item) ::
     &     fl_en_res, fl_spc(nspcfrm)

      integer, external ::
     &     idx_formlist, idx_mel_list, idx_xret
      logical, external ::
     &     file_exists
      real(8), external ::
     &     xnormop

      ifree = mem_setmark('solve_nleq')

      if (iprlvl.ge.5) then
        write(lulog,*) 'formula: ',trim(label_form)
        if (nopt.gt.1)
     &       write(lulog,*) 'solving for ',nopt,
     &       ' operators simultaneously'
        write(lulog,*)   'solving for: ',trim(label_opt(1))
        do iopt = 2, nopt
          write(lulog,*) '             ',trim(label_opt(iopt))
        end do
      end if

      idx = idx_formlist(label_form,form_info)
      if (idx.le.0)
     &     call quit(1,'solve_nleq',
     &     'did not find formula '//trim(label_form))
      form_en_res => form_info%form_arr(idx)%form

      allocate(ffopt(nopt),ffdia(nopt),ffgrd(nopt),ffspecial(nspecial),
     &     me_opt(nopt),me_dia(nopt),me_grd(nopt),me_special(nspecial))

      do iopt = 1, nopt
        jopt = iopt
        idxmel = idx_mel_list(label_opt(iopt),op_info)
        ierr = 1
        if (idxmel.le.0) exit
        me_opt(iopt)%mel =>  op_info%mel_arr(idxmel)%mel
        ffopt(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand
        ierr = 2
        if (.not.associated(ffopt(iopt)%fhand)) exit
        idxmel = idx_mel_list(label_res(iopt),op_info)
        ierr = 3
        if (idxmel.le.0) exit
        me_grd(iopt)%mel =>  op_info%mel_arr(idxmel)%mel
        ffgrd(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand
        ierr = 4
        if (.not.associated(ffgrd(iopt)%fhand)) exit
        idxmel = idx_mel_list(label_prc(iopt),op_info)
        ierr = 5
        if (idxmel.le.0) exit
        me_dia(iopt)%mel =>  op_info%mel_arr(idxmel)%mel
        ffdia(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand
        ierr = 6
        if (.not.associated(ffdia(iopt)%fhand)) exit
        ierr = 0
      end do

      ! special lists needed?
      if (ierr.eq.0) then
        do idx = 1, nspecial
          jopt = idx
          idxmel = idx_mel_list(label_special(idx),op_info)
          ierr = 7
          if (idxmel.le.0) exit
          me_special(idx)%mel  => op_info%mel_arr(idxmel)%mel
          ffspecial(idx)%fhand => op_info%mel_arr(idxmel)%mel%fhand
          ierr = 8
          if (.not.associated(ffspecial(idx)%fhand)) exit
          ierr = 0
        end do
      end if

      ! error handling
      if (ierr.gt.0) then
        if (ierr.eq.1.or.ierr.eq.2) label = label_opt(jopt)
        if (ierr.eq.3.or.ierr.eq.4) label = label_res(jopt)
        if (ierr.eq.5.or.ierr.eq.6) label = label_prc(jopt)
        if (ierr.eq.7.or.ierr.eq.8) label = label_special(jopt)

        if (mod(ierr,2).eq.1)
     &       call quit(1,'solve_nleq',
     &       'did not find list "'//trim(label)//'"')
        if (mod(ierr,2).eq.0)
     &       call quit(1,'solve_nleq',
     &       'no file associated to list '//trim(label))
      end if

      ! special formulae
      do jdx = 1, nspcfrm
        idx = idx_formlist(label_spcfrm(jdx),form_info)
        if (idx.le.0)
     &       call quit(1,'solve_nleq',
     &       'did not find formula '//trim(label_spcfrm(jdx)))
        ! read formula
        call read_form_list(form_info%form_arr(idx)%form%fhand,
     &                      fl_spc(jdx),.true.)
      end do

      ! for safety reasons, we allocate the two guys
      allocate(me_trv(1),me_h_trv(1))
      call get_argument_value('calculate.solve.non_linear',
     &                        'restart',lval=restart)

      do iopt = 1, nopt
        ! open result vector file(s)
cmh     if file already open, use as initial guess!
        if (ffopt(iopt)%fhand%unit.gt.0) then
          write(lulog,'(x,a,i1,a)')
     &         'Using existing amplitudes as initial guess for vector ',
     &         iopt,'!'
        else
          call file_open(ffopt(iopt)%fhand)
          ! hard restart? (just use old amplitude file as initial guess)
          if (restart) then
            inquire(file=trim(ffopt(iopt)%fhand%name),exist=restart)
            if (.not.restart) call warn('solve_nleq',
     &         'No amplitude file found for restart! Setting to zero.')
          end if
          if (restart) then
            write(lulog,'(x,a,i1,a)')
     &         'Using old amplitude file for vector ',iopt,'!'
          else
            call zeroop(me_opt(iopt)%mel)
          end if
        end if
        ! open corresponding residuals ...
        call file_open(ffgrd(iopt)%fhand)
        ! ... and corresponding preconditioner(s)
        if (ffdia(iopt)%fhand%unit.le.0)
     &       call file_open(ffdia(iopt)%fhand)
      end do

      do idx = 1, nspecial
        if (ffspecial(idx)%fhand%unit.le.0)
     &       call file_open(ffspecial(idx)%fhand)
      end do
      
cmh      ! get initial amplitudes
cmh      do iopt = 1, nopt
cmhc        if (.not.file_exists(me_opt(iopt)%mel%fhand)) then
cmh          call zeroop(me_opt(iopt)%mel)
cmhc dbg DBGDBG!!!
cmhc        else
cmhc          call warn('solve_nleq','debug version of restart active')
cmhc        end if
cmh      end do

      ! use mode_str to set special preconditioning, e.g. for R12

      call set_opti_info(opti_info,1,nopt,1,me_opt,mode_str)

      ! write information to opti_info about signs which might occur
      ! in preconditioning step
      ! relevant when amp is njoined=1 op. but grd is njoined=2 op.
      call set_opti_info_signs(opti_info,1,nopt,
     &                         me_opt,me_grd,me_grd,me_grd,.false.)

      ! read formula
      call read_form_list(form_en_res%fhand,fl_en_res,.true.)

      ! set dependency info for submitted formula list
      call set_formula_dependencies(depend,fl_en_res,op_info)

      if (opti_info%skip_resx) then
        ! exclude residuals which are not needed from dependency info
        ! (these residuals are also defined as up to date)
        idxmel = idx_mel_list(label_res(nopt),op_info)
        call trunc_formula_dependencies(depend,idxmel,op_info)
      end if

      ! number of info values returned on xret
      nout = depend%ntargets
      allocate(xret(nout))

      ! find out, which entries of xret are the ones that we need
      idx_en_xret = idx_xret(label_en,op_info,depend)
c dbg
c      print *,'idx_en_xret: ',idx_en_xret
c dbg
      if (idx_en_xret.le.0)
     &     call quit(1,'solve_nleq',
     &     'formula does not provide an update for the energy')
      do iopt = 1, nopt
        idx_res_xret(iopt) = idx_xret(label_res(iopt),op_info,depend)
        if (idx_res_xret(iopt).le.0)
     &       call quit(1,'solve_nleq',
     &       'formula does not provide an update for all residuals')
      end do

      ! start optimization loop
      imacit = 0
      imicit = 0
      imicit_tot = 0
      task = 0
      opt_loop: do !while(task.lt.8)

        call optcont
     &       (imacit,imicit,imicit_tot,
     &       task,conv,
     &       energy,xresnrm,
     &       me_opt,me_grd,me_dia,
     &       me_trv,me_h_trv,
     &       me_special, nspecial,! <- R12: pass B, X, H here
c     &       ffopt,ffgrd,ffdia,ffmet, ! <- R12: pass X here (metric)
c     &       ff_trv,ff_h_trv,
     &       fl_spc,nspcfrm,
     &       opti_info,opti_stat,
     &       orb_info,op_info,str_info,strmap_info)

        ! output
        it_print = imacit-1
        if (conv) it_print = imacit
        if (luout.ne.lulog) then
          if (nopt.eq.1.and.imacit.gt.1)
     &       write(luout,'(3x,i3,f24.12,x,g10.4)')
     &       it_print,energy,xresnrm(1)
          if (nopt.eq.2.and.imacit.gt.1)
     &       write(luout,'(3x,i3,f24.12,2(x,g10.4))')
     &       it_print,energy,xresnrm(1:2)
          if (nopt.eq.3.and.imacit.gt.1)
     &       write(luout,'(3x,i3,f24.12,3(x,g10.4))')
     &       it_print,energy,xresnrm(1:3)
          if (.not.conv.and.task.ge.8) then
            write(luout,'("    NOT CONVERGED!")')
          else if (task.ge.8) then
            write(luout,'("    CONVERGED")')
C            write(luout,'(">>> final energy:",f24.12," <<<")')
C     &       energy
          end if
        end if

        if (nopt.eq.1.and.imacit.gt.1)
     &     write(lulog,'(">>>",i3,f24.12,x,g10.4)')
     &     it_print,energy,xresnrm(1)
        if (nopt.eq.2.and.imacit.gt.1)
     &     write(lulog,'(">>>",i3,f24.12,2(x,g10.4))')
     &     it_print,energy,xresnrm(1:2)
        if (nopt.eq.3.and.imacit.gt.1)
     &     write(lulog,'(">>>",i3,f24.12,3(x,g10.4))')
     &     it_print,energy,xresnrm(1:3)
        if (.not.conv.and.task.ge.8) then
          write(lulog,'(">>> NOT CONVERGED! <<<")')
          exit opt_loop
        else if (task.ge.8) then
          write(lulog,'(">>> CONVERGED <<<")')
          write(lulog,'(">>> final energy:",f24.12," <<<")')
     &     energy
           exit opt_loop
        end if

        ! quick and dirty (for experimental use):
        ! do C0 optimization if requested
        if (opti_info%optref.eq.-3.and.
     &      (imacit.gt.1.or.restart.or.opti_info%skip_resx)
     &      .and..not.conv.and.nspcfrm.gt.1) then
          call get_argument_value('method.MR','ciroot',
     &       ival=idx)
          call get_argument_value('method.MR','maxroot',
     &       ival=ndx)
          if(ndx.le.0) ndx=idx
          call me_list_label(dia_label,'DIA',orb_info%lsym,
     &                       0,0,0,.false.)
          dia_label = trim(dia_label)//'C0'
          ! use weaker convergence threshold for micro-iterations
          thr_suggest = min(xresnrm(1)*opti_info%mic_ahead,1d-4)
          if (spinadapt.gt.0) then
            call solve_evp('SPP',1,ndx,idx,
     &                 'ME_C0',trim(dia_label),'A_C0',
     &                 'C0','FOPT_OMG_C0',
     &                 'ME_C0_sp',1,
     &                 'FOPT_C0_sp',1,thr_suggest,
     &                 op_info,form_info,str_info,strmap_info,orb_info)
          else
            call solve_evp('DIA',1,ndx,idx,
     &                 'ME_C0',trim(dia_label),'A_C0',
     &                 'C0','FOPT_OMG_C0','-',0,'-',0,thr_suggest,
     &                 op_info,form_info,str_info,strmap_info,orb_info)
          end if
c dbg
c          idx = idx_mel_list('ME_C0',op_info)
c          print *,'current C0 vector: '
c          call wrt_mel_file(lulog,1000,op_info%mel_arr(idx)%mel,
c     &       1,op_info%mel_arr(idx)%mel%op%n_occ_cls,
c     &       str_info,orb_info)
c dbgend
        end if

c dbg
c        if (opti_info%typ_prc(1).eq.optinf_prc_traf
c     &      .and.nopt.ge.2) then
c          print *,'current C0 vector: '
c          call wrt_mel_file(lulog,1000,me_opt(2)%mel,
c     &       1,me_opt(2)%mel%op%n_occ_cls,
c     &       str_info,orb_info)
c        end if
c dbgend

        if (ntest.ge.1000) then
          do iopt = 1, nopt
            write(lulog,*) 'dump of '//trim(me_opt(iopt)%mel%label)
            write(lulog,*) 'iopt = ',iopt
            call wrt_mel_file(lulog,5,
     &           me_opt(iopt)%mel,
     &           1,me_opt(iopt)%mel%op%n_occ_cls,
     &           str_info,orb_info)
          end do
        end if

        ! here?
        do iopt = 1, nopt
          if (opti_info%typ_prc(iopt).eq.optinf_prc_norm.and.
     &        opti_info%optref.ne.0) cycle !no unneeded metric update
          call touch_file_rec(ffopt(iopt)%fhand)
        end do
c dbg
c        if (opti_info%typ_prc(1).eq.optinf_prc_traf.and.
c     &      nspecial.ge.7)
c     &      call touch_file_rec(me_special(7)%mel%fhand)
c        if (opti_info%typ_prc(1).eq.optinf_prc_traf.and.
c     &      nspecial.ge.6.and.nopt.eq.1) then
c          ! very dirty trick: prevent calc. of Res.#2
c            call touch_file_rec(op_info%mel_arr(
c     &           idx_mel_list('ME_A_C0',op_info))%mel%fhand)
c          ! also delete information that Res.#2 depends on energy
c          if (imacit.eq.1) then
c            do idx = 1, depend%ndepend
c              if (depend%depends_on_idxlist(idx,3).eq.
c     &            depend%idxlist(1))
c     &           depend%depends_on_idxlist(idx,3) = 0
c            end do
c          end if
c        end if
c dbgend

        ! 1 - get energy
        ! 2 - get residual
        if (iand(task,1).eq.1.or.iand(task,2).eq.2) then
          call frm_sched(xret,fl_en_res,depend,0,0,
     &         .true.,.false.,op_info,str_info,strmap_info,orb_info)
          ! intermediates should be generated first, energy
          ! is expected to be the last "intermediate"
          energy =  xret(idx_en_xret)
c dbg
c          print *,'xret : ',xret
c dbg

          if (ntest.ge.1000) then
            do iopt = 1, nopt
              write(lulog,*) 'dump of residual '//
     &             trim(me_grd(iopt)%mel%label)
              call wrt_mel_file(lulog,5,
     &           me_grd(iopt)%mel,
     &           1,me_grd(iopt)%mel%op%n_occ_cls,
     &           str_info,orb_info)
            end do
          end if
c dbg:
c test
c          if (imacit.ge.3) then
c            print *,'>> put O to zero'
c            do idx=1,20
c              print *,'   put O to zero'
c            end do
c            call scale_opblk(xdum,0d0,me_grd(1)%mel,me_grd(1)%mel,
c     &           1,1,orb_info)
c          end if
c test

          do iopt = 1, nopt
            xresnrm(iopt) = abs(xret(idx_res_xret(iopt)))
            traf = traf.or.opti_info%typ_prc(iopt).eq.optinf_prc_traf
     &                 .or.opti_info%typ_prc(iopt).eq.optinf_prc_invH0
          end do
        end if

        ! report untransformed residual
        if (traf) 
     &   write(lulog,'(x,"norm of untransformed residual ",3(x,g10.4))')
     &   xresnrm(1:opti_info%nopt)

        ! another quick and dirty call to the linear solver
        ! for advanced preconditioning
        if (opti_info%optref.eq.-3.and.
     &      opti_info%typ_prc(1).eq.optinf_prc_invH0.and.
     &      .not.conv) then
          thr_suggest = min(xresnrm(1)*opti_info%mic_ahead,1d-4)

          idx = idx_mel_list('ME_C0',op_info)
          if (opti_info%optref.ne.0.and.
     &         op_info%mel_arr(idx)%mel%fhand%last_mod(
     &         op_info%mel_arr(idx)%mel%fhand%current_record).gt.
     &         me_special(2)%mel%fhand%last_mod(1)) then
            call update_metric(me_dia(1)%mel,me_special,nspecial,
     &           fl_spc,nspcfrm,orb_info,op_info,str_info,strmap_info,
     &           opti_info%update_prc.gt.0.and.
     &           mod(imacit,max(opti_info%update_prc,1)).eq.0)
          end if

          call solve_leq('TRF',
     &                 1,1,'ME_DlT',label_prc(1),'H0_DlT','DlT',!'S_DlT',
     &                 'OMGprj',xdum,'FOPT_H0INV',
     &            (/'ME_Tout  ','ME_Ttr   ','ME_Dtr   ','ME_Dtrdag'/),4,
     &                 'FOPT_Ttr_GEN',1,  thr_suggest,
     &                 op_info,form_info,str_info,strmap_info,orb_info)
          ! overwrite the gradient with preconditioned gradient
          call scale_copy_op(label_res(1),'ME_DlT',1d0,1,'--',0,
     &                  op_info,orb_info,str_info)
          ! get norm of projected gradient and report as gradient norm
          xresnrm(1) = xdum
        end if

      end do opt_loop

      call clean_formula_dependencies(depend)

      ! close files
      do iopt = 1, nopt
        ! open result vector file(s)
c dbg - for R12:
c        if (iopt.eq.2) then
c          write(lulog,*) ' iopt = ',iopt
c          call wrt_mel_file(lulog,1000,me_opt(iopt)%mel,
c     &       1,me_opt(iopt)%mel%op%n_occ_cls,
c     &       str_info,orb_info)        
c        end if
c dbg
c dbg
        ! very dirty: don't close file if needed for following opt.
        if ((opti_info%typ_prc(1).ne.optinf_prc_traf.and.
     &       opti_info%typ_prc(1).ne.optinf_prc_invH0).or.
     &      opti_info%optref.eq.0.or.nopt.ne.1) then
c dbgend
        call file_close_keep(ffopt(iopt)%fhand)
c dbg
        end if
c dbgend
        ! open corresponding residuals ...
        call file_close_keep(ffgrd(iopt)%fhand)
        ! ... and corresponding preconditioner(s)
        if (ffdia(iopt)%fhand%unit.gt.0)
     &       call file_close_keep(ffdia(iopt)%fhand)
      end do

      do idx = 1, nspecial
        if (ffspecial(idx)%fhand%unit.gt.0)
     &       call file_close_keep(ffspecial(idx)%fhand)
      end do

      deallocate(ffopt,ffdia,ffgrd,ffspecial,
     &     me_opt,me_dia,me_grd,me_special,me_trv,me_h_trv,xret)
      call dealloc_formula_list(fl_en_res)
      do jdx = 1, nspcfrm
        call dealloc_formula_list(fl_spc(jdx))
      end do
      ifree = mem_flushmark()

      return
      end

