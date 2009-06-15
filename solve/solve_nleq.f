*----------------------------------------------------------------------*
      subroutine solve_nleq(mode_str,
     &     nopt,label_opt,label_res,label_prc,label_en,
     &     label_form,
     &     label_special,nspecial,           !<- eg. for R12
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

      integer, parameter ::
     &     ntest = 000

      integer, intent(in) ::
     &     nopt, nspecial
      character(*), intent(in) ::
     &     mode_str,
     &     label_opt(nopt),
     &     label_res(nopt),
     &     label_prc(nopt),
     &     label_special(nspecial),
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
     &     conv
      character(len_opname) ::
     &     label
      integer ::
     &     imacit, imicit, imicit_tot, iprint, task, ifree, iopt, jopt,
     &     idx, idxmel, ierr, nout, idx_en_xret, idx_res_xret(nopt)
      real(8) ::
     &     energy, xresnrm(nopt), xdum
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
     &     fl_en_res

      integer, external ::
     &     idx_formlist, idx_mel_list, idx_xret
      logical, external ::
     &     file_exists
      real(8), external ::
     &     xnormop

      ifree = mem_setmark('solve_nleq')

      if (iprlvl.ge.5) then
        write(luout,*) 'formula: ',trim(label_form)
        if (nopt.gt.1)
     &       write(luout,*) 'solving for ',nopt,
     &       ' operators simultaneously'
        write(luout,*)   'solving for: ',trim(label_opt(1))
        do iopt = 2, nopt
          write(luout,*) '             ',trim(label_opt(iopt))
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

      ! for safety reasons, we allocate the two guys
      allocate(me_trv(1),me_h_trv(1))

      do iopt = 1, nopt
        ! open result vector file(s)
        call file_open(ffopt(iopt)%fhand)
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
      
      ! get initial amplitudes
      do iopt = 1, nopt
c        if (.not.file_exists(me_opt(iopt)%mel%fhand)) then
          call zeroop(me_opt(iopt)%mel)
c dbg DBGDBG!!!
c        else
c          call warn('solve_nleq','debug version of restart active')
c        end if
      end do

      ! use mode_str to set special preconditioning, e.g. for R12

      call set_opti_info(opti_info,1,nopt,1,me_opt,mode_str)

      ! read formula
      call read_form_list(form_en_res%fhand,fl_en_res)

      ! set dependency info for submitted formula list
      call set_formula_dependencies(depend,fl_en_res,op_info)

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
      opt_loop: do while(task.lt.8)

        call optcont
     &       (imacit,imicit,imicit_tot,
     &       task,conv,
     &       energy,xresnrm,
     &       me_opt,me_grd,me_dia,
     &       me_trv,me_h_trv,
     &       me_special, nspecial,! <- R12: pass B, X, H here
c     &       ffopt,ffgrd,ffdia,ffmet, ! <- R12: pass X here (metric)
c     &       ff_trv,ff_h_trv,
     &       opti_info,opti_stat,
     &       orb_info,op_info,str_info,strmap_info)

        if (ntest.ge.1000) then
          do iopt = 1, nopt
            write(luout,*) 'dump of '//trim(me_opt(iopt)%mel%label)
            write(luout,*) 'iopt = ',iopt
            call wrt_mel_file(luout,5,
     &           me_opt(iopt)%mel,
     &           1,me_opt(iopt)%mel%op%n_occ_cls,
     &           str_info,orb_info)
          end do
        end if

        ! here?
        do iopt = 1, nopt
          call touch_file_rec(ffopt(iopt)%fhand)
        end do

        ! 1 - get energy
        ! 2 - get residual
        if (iand(task,1).eq.1.or.iand(task,2).eq.2) then
          call frm_sched(xret,fl_en_res,depend,0,0,
     &         op_info,str_info,strmap_info,orb_info)
          ! intermediates should be generated first, energy
          ! is expected to be the last "intermediate"
          energy =  xret(idx_en_xret)
c dbg
c          print *,'xret : ',xret
c dbg

          if (ntest.ge.1000) then
            do iopt = 1, nopt
              write(luout,*) 'dump of residual '//
     &             trim(me_grd(iopt)%mel%label)
              call wrt_mel_file(luout,5,
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
          end do
        end if

        if (.not.conv.and.task.lt.8) then
          if (nopt.eq.1)
     &       write(luout,'(">>>",i3,f24.12,x,g10.4)')
     &       imacit,energy,xresnrm(1)
          if (nopt.eq.2)
     &       write(luout,'(">>>",i3,f24.12,2(x,g10.4))')
     &       imacit,energy,xresnrm(1:2)
          if (nopt.eq.3)
     &       write(luout,'(">>>",i3,f24.12,3(x,g10.4))')
     &       imacit,energy,xresnrm(1:3)
        else if (.not.conv) then
          write(luout,'(">>> NOT CONVERGED! <<<")')
        else
          write(luout,'(">>> final energy:",f24.12," <<<")')
     &       energy
        end if

      end do opt_loop

      call clean_formula_dependencies(depend)

      ! close files
      do iopt = 1, nopt
        ! open result vector file(s)
c dbg - for R12:
        if (iopt.eq.2) then
          write(luout,*) ' iopt = ',iopt
          call wrt_mel_file(luout,1000,me_opt(iopt)%mel,
     &       1,me_opt(iopt)%mel%op%n_occ_cls,
     &       str_info,orb_info)        
        end if
c dbg
        call file_close_keep(ffopt(iopt)%fhand)
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
      ifree = mem_flushmark()

      return
      end

