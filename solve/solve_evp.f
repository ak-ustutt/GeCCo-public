*----------------------------------------------------------------------*
      subroutine solve_evp(mode_str,
     &     nopt,nroots,label_opt,label_prc,label_op_mvp,
     &     label_form,
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
*
*     label_opt(nopt)       label of solution vectors
*     label_prc(nopt)       label of preconditioners
*     label_op_mvp(nopt)    label operators describing Mx-products
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
     &     nopt, nroots
      character(*), intent(in) ::
     &     mode_str,
     &     label_opt(nopt),
     &     label_prc(nopt),
     &     label_op_mvp(nopt),
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
     &     iter, iprint, task, ifree, iopt, jopt, nintm, irequest,
     &     nrequest, nvectors, iroot, idx, ierr, idxmel, nout,
     &     idxlist(2*nroots)
      real(8) ::
     &     xresmax,
     &     xeig(nroots,2), xresnrm(nroots), xlist(2*nroots)
      type(me_list_array), pointer ::
     &     me_opt(:), me_dia(:), me_trv(:), me_mvp(:)
      type(file_array), pointer ::
     &     ffdia(:), ff_trv(:),
     &     ffopt(:), ff_mvp(:)
      type(dependency_info) ::
     &     depend
      type(optimize_info) ::
     &     opti_info
      type(optimize_status) ::
     &     opti_stat
      type(formula), pointer ::
     &     form_mvp
      type(formula_item) ::
     &     fl_mvp

      integer, pointer ::
     &     irecmvp(:), irectrv(:)
      real(8), pointer ::
     &      xret(:)

      character ::
     &     fname*256

      integer, external ::
     &     idx_formlist, idx_mel_list, idx_xret
      real(8), external ::
     &     fndmnx

      ifree = mem_setmark('solve_leq')

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'entered solve_evp')
        write(luout,*) 'nopt   = ',nopt
        write(luout,*) 'nroots = ',nroots
      end if

      if (nopt.gt.1)
     &     call quit(1,'solve_evp','did not yet consider coupled EVPs')

      idx = idx_formlist(label_form,form_info)
      if (idx.le.0)
     &     call quit(1,'solve_evp',
     &     'did not find formula '//trim(label_form))
      form_mvp => form_info%form_arr(idx)%form


      allocate(me_opt(nopt),me_dia(nopt),me_trv(nopt),me_mvp(nopt))
      allocate(ffopt(nopt),ffdia(nopt),
     &     ff_trv(nopt),ff_mvp(nopt))
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

      ! error handling
      if (ierr.gt.0) then
        if (ierr.eq.1.or.ierr.eq.2) label = label_opt(jopt)
        if (ierr.eq.3.or.ierr.eq.4) label = label_prc(jopt)
        if (mod(ierr,2).eq.1)
     &       call quit(1,'solve_evp',
     &       'did not find list '//trim(label))
        if (mod(ierr,2).eq.0)
     &       call quit(1,'solve_evp',
     &       'no file associated to list '//trim(label))
      end if

      call set_opti_info(opti_info,3,nopt,nroots,me_opt,mode_str)

      nvectors = opti_info%maxsbsp

      do iopt = 1, nopt
        ! get a ME-list for trial-vectors
        write(fname,'("trv_",i3.3)') iopt
        call define_me_list(fname,me_opt(iopt)%mel%op%name,
     &       me_opt(iopt)%mel%absym,me_opt(iopt)%mel%casym,
     &       me_opt(iopt)%mel%gamt,me_opt(iopt)%mel%s2,
     &       me_opt(iopt)%mel%mst,
     &       1,nvectors,
     &       op_info,orb_info,str_info,strmap_info)
        idxmel = idx_mel_list(fname,op_info)
        me_trv(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        ff_trv(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand

        ! get a ME list for matrix-vector products
        ! (have same symmtry properties as result!)
        write(fname,'("mvp_",i3.3)') iopt
        call define_me_list(fname,label_op_mvp,
     &       me_opt(iopt)%mel%absym,me_opt(iopt)%mel%casym,
     &       me_opt(iopt)%mel%gamt,me_opt(iopt)%mel%s2,
     &       me_opt(iopt)%mel%mst,
     &       1,nvectors,
     &       op_info,orb_info,str_info,strmap_info)
        idxmel = idx_mel_list(fname,op_info)
        me_mvp(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        ff_mvp(iopt)%fhand => op_info%mel_arr(idxmel)%mel%fhand

      end do

      ! read formula
      call read_form_list(form_mvp%fhand,fl_mvp)

      ! set dependency info for submitted formula list
      call set_formula_dependencies(depend,fl_mvp,op_info)

      ! number of info values returned on xret
      nout = depend%ntargets
      allocate(xret(nout))

      ! records with trial vectors and Mv-products, needed in evp_control:
      ifree = mem_alloc_int(irectrv,nroots,'rectrv')
      ifree = mem_alloc_int(irecmvp,nroots,'recmvp')

      do iopt = 1, nopt
        ! open result vector file(s)
        call file_open(ffopt(iopt)%fhand)
        call file_open(ff_trv(iopt)%fhand)
        ! open corresponding matrix vector products ...
        call file_open(ff_mvp(iopt)%fhand)
        ! ... and corresponding preconditioner(s)
        if (ffdia(iopt)%fhand%unit.le.0)
     &       call file_open(ffdia(iopt)%fhand)
      end do

      ! get initial amplitudes
      do iopt = 1, nopt
        call find_nmin_list(xlist,idxlist,2*nroots,me_dia(iopt)%mel)
        do iroot = 1, nroots
          
          call switch_mel_record(me_trv(iopt)%mel,iroot)
          call diag_guess(me_trv(iopt)%mel,
     &         xlist,idxlist,2*nroots,iroot,0)

        end do
      end do

      ! start optimization loop
      iter = 0
      task = 0
      opt_loop: do while(task.lt.8)

        call leq_evp_control
     &       ('EVP',iter,
     &       task,conv,xresnrm,xeig,
     &       nrequest,irectrv,irecmvp,
     &       ffopt,ff_trv,ff_mvp,ffdia,ffdia,  ! #4 is dummy
     &       opti_info,opti_stat)

        if (iter.gt.1) then
          xresmax = fndmnx(xresnrm,nroots,2)
          write(luout,'(">>>",i3,24x,x,g10.4)') iter-1,xresmax
          if (iprlvl.gt.0) then
            do iroot = 1, nroots
              if (xeig(iroot,2).eq.0d0) then
                write(luout,'(" >>",3x,f24.12,x,g10.4)')
     &               xeig(iroot,1),xresnrm(iroot)
              else
                write(luout,
     &               '(" >>",3x,f24.12,x,g10.4," (img=",g24.12,")")')
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
     &           op_info,str_info,strmap_info,orb_info)

c dbg
c            write(luout,*) 'output for request: ',irequest
c            call wrt_mel_file(luout,5,me_mvp(1)%mel,
c     &           1,me_mvp(1)%mel%op%n_occ_cls,
c     &           str_info,orb_info)
c dbg

          end do
        end if

      end do opt_loop

      do iopt = 1, nopt

        ! remove the temporary lists
        call del_me_list(me_trv(iopt)%mel%label,op_info)
        call del_me_list(me_mvp(iopt)%mel%label,op_info)

        ! make sure that the operator is now associated with
        ! the list containing the solution vector
        call assign_me_list(label_opt(iopt),
     &                      me_opt(iopt)%mel%op%name,op_info)

      end do


      ! note that only the pointer array ffopt (but not the entries)
      ! is deallocated:
      deallocate(me_opt,me_dia,me_trv,me_mvp)
      deallocate(ff_trv,ff_mvp,ffdia,ffopt,xret)

      ifree = mem_flushmark()

      return
      end


