*----------------------------------------------------------------------*
      subroutine init_guess(nopt,init,nroots,
     &                me_opt,me_trv,me_dia,me_special,nspecial,
     &                fl_mvp,depend,fl_spc,nspcfrm,
     &                opti_info,orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*
* get suitable initial guess vectors for the EVP solver
*
* matthias, nov 2012
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
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_file_array.h'
      include 'mdef_formula_info.h'
      include 'def_dependency_info.h'
      include 'ifc_memman.h'

      integer, intent(in) ::
     &     nopt, nroots, nspecial, nspcfrm
      logical, intent(in) ::
     &     init(nopt)
      type(me_list_array), intent(inout) :: 
     &     me_opt(*), me_trv(*), me_dia(*), me_special(nspecial)
      type(formula_item), intent(inout) ::
     &     fl_mvp, fl_spc(nspcfrm)
      type(dependency_info), intent(in) ::
     &     depend
      type(optimize_info), intent(in) ::
     &     opti_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info

      integer ::
     &     ifree, iopt, iroot, ntrials, nout,
     &     idxlist(max(100,4*nroots)), nselect, iguess
      logical ::
     &     trafo
      real(8) ::
     &     xnrm, xretlast, xover, xlist(max(100,4*nroots))
      type(me_list), pointer ::
     &     me_pnt

      integer, pointer ::
     &     idxselect(:)
      real(8), pointer ::
     &      xret(:), xbuf1(:), xbuf2(:)

      real(8), external ::
     &     da_ddot

      ntrials = max(100,4*nroots)
      nout = depend%ntargets
      allocate(xret(nout))

      do iopt = 1, nopt
        if (.not.init(iopt)) then
          do iroot = 1, nroots
            call switch_mel_record(me_trv(iopt)%mel,iroot)
            call switch_mel_record(me_opt(iopt)%mel,iroot)
            call list_copy(me_opt(iopt)%mel,me_trv(iopt)%mel,.false.)
          end do
          cycle
        end if
        ! preliminary solution: set only component 1, rest is zero
        if (iopt.gt.1) then
          do iroot = 1, nroots
            call switch_mel_record(me_trv(iopt)%mel,iroot)
            call zeroop(me_trv(iopt)%mel)
          end do
          cycle
        end if
        ! transformed preconditioner => transformed initial guess vector
        if (opti_info%typ_prc(iopt).eq.optinf_prc_traf) then
          me_pnt => me_special(1)%mel
          trafo = .true.
        else
          me_pnt => me_trv(iopt)%mel
          trafo = .false.
        end if
        call find_nmin_list(xlist,idxlist,ntrials,me_dia(iopt)%mel)
        iroot = 0
        do iguess = 1, ntrials
          iroot = iroot + 1    

          call switch_mel_record(me_pnt,iroot)
          call diag_guess(me_pnt,
     &         xlist,idxlist,ntrials,iguess,me_pnt%absym,
     &         op_info,str_info,strmap_info,orb_info)

          ! if requested, back-transformation of initial guess vector
          if (trafo) then
            ! use non-daggered transformation matrix if requested
            if (nspecial.eq.3)
     &         call assign_me_list(me_special(2)%mel%label,
     &                             me_special(2)%mel%op%name,op_info)
            ! do the transformation
            allocate(idxselect(nout))
            nselect = 0
            call select_formula_target(idxselect,nselect,
     &                  me_trv(iopt)%mel%label,depend,op_info)
            call switch_mel_record(me_trv(iopt)%mel,iroot)
            call frm_sched(xret,fl_mvp,depend,idxselect,nselect,
     &             .true.,.false.,op_info,str_info,strmap_info,orb_info)
            ! guess vectors of wrong spin symmetry will be discarded
            if (abs(xret(idxselect(1))).lt.1d-12) then
              if (iprlvl.ge.5) write(luout,*)
     &           'Discarding guess vector with wrong spin symmetry.'
              me_trv(iopt)%mel%fhand%last_mod(iroot) = -1
              iroot = iroot - 1
            else 
              if (abs(xret(idxselect(1))).lt.1d0-1d-12)
     &            call warn('init_guess','guess vector not normalized')
                
              if (iroot.gt.1) then
                ! Due to symmetrization we might get same guess twice
                ifree = mem_setmark('init_guess.check_guess')
                ifree = mem_alloc_real(xbuf1,opti_info%nwfpar(iopt),
     &                                 'xbuf1')
                ifree = mem_alloc_real(xbuf2,opti_info%nwfpar(iopt),
     &                                 'xbuf2')
                xover = da_ddot(me_trv(iopt)%mel%fhand,iroot-1,1,
     &                            me_trv(iopt)%mel%fhand,iroot,1,
     &                            opti_info%nwfpar(iopt),
     &                            xbuf1,xbuf2,
     &                            opti_info%nwfpar(iopt))
                if (abs(abs(xover)-xretlast).lt.1d-6) then
                  if (iprlvl.ge.5) write(luout,*)
     &                'Discarding twin guess vector.'
                  me_trv(iopt)%mel%fhand%last_mod(iroot) = -1
                  iroot = iroot - 1
                end if  
                ifree = mem_flushmark()
              end if
              xretlast = xret(idxselect(1))**2
            end if
            deallocate(idxselect)
          end if

          ! project out spin contaminations?
          if (opti_info%typ_prc(iopt).eq.optinf_prc_spinp) then
            ifree = mem_setmark('init_guess.spin_proj')
            ifree = mem_alloc_real(xbuf1,opti_info%nwfpar(iopt),'xbuf1')
            ifree = mem_alloc_real(xbuf2,opti_info%nwfpar(iopt),'xbuf2')
            call spin_project(me_trv(iopt)%mel,me_special(1)%mel,
     &                        fl_spc(1),opti_info%nwfpar(iopt),
     &                        xbuf1,xbuf2,.true.,xnrm,
     &                        opti_info,orb_info,
     &                        op_info,str_info,strmap_info)
            ifree = mem_flushmark()
            if (xnrm.lt.1d-12) then
              if (iprlvl.ge.5) write(luout,*)
     &           'Discarding guess vector with wrong spin symmetry.'
              me_trv(iopt)%mel%fhand%last_mod(iroot) = -1
              iroot = iroot - 1
            end if
          end if

          if (iroot.eq.nroots) exit

c dbg
c          if (file_exists(me_opt(iopt)%mel%fhand)) then
c            print *,' *** RESTARTING ***'
c            call switch_mel_record(me_opt(iopt)%mel,iroot)
c            call list_copy(me_opt(iopt)%mel,me_trv(iopt)%mel,.false.)
c          end if
c dbg
        end do
        if (iroot.ne.nroots) call quit(1,'init_guess',
     &        'Could not find enough guess vectors')
      end do

      deallocate(xret)

      return
      end
