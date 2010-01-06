*----------------------------------------------------------------------*
      subroutine process_me_lists(rule,
     &     form_info,op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     evaluate/import elements of ME-list according to action "rule"
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_target.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_strmapinf.h'
      include 'par_actions.h'     
      include 'ifc_input.h'

      type(action), intent(in) ::
     &     rule
      type(formula_info), intent(inout) ::
     &     form_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(inout) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info

      integer, parameter ::
     &     maxfac = 20, maximum_order = 10
      real(8) ::
     &     fac(maxfac), freq
      integer ::
     &     idxblk(maxfac), minblk, maxblk,
     &     idx, jdx, ioff, nfac, nblk, nspecial, imode,
     &     absym,casym,gamma,s2,ms,nopt,nroots,ndens,rank
      logical ::
     &     ms_fix, form_test
      type(me_list), pointer ::
     &     mel_pnt
      character(len_command_par) ::
     &     title, env_type, list_type, mode
      character(len_command_par), allocatable ::
     &     label_met(:)

      integer ::
     &     idx_formlist, order, dummy, loopdum
      integer, allocatable ::
     &     ifreq_temp(:), ifreq(:)

      integer, external ::
     &     idx_mel_list

      if (rule%type.ne.ttype_opme)
     &     call quit(1,'process_me_lists',
     &               'called for wrong target type')
      if (rule%n_labels.lt.1)
     &     call quit(1,'process_me_lists',
     &               'at least one label expected')

      ! form_test = true skips time consuming steps
      call get_argument_value('general','form_test',lval=form_test)
      loop:do loopdum = 1,1

      select case(trim(rule%command))
      case(DEF_ME_LIST)
        if (rule%n_labels.lt.2)
     &     call quit(1,'process_me_lists','two labels expected')

        call me_list_parameters(+1,rule%parameters,
     &       absym,casym,gamma,s2,ms,ms_fix)

        call define_me_list(rule%labels(1),rule%labels(2),
     &       absym,casym,gamma,s2,ms,ms_fix,
     &       -1,-1,
     &       op_info,orb_info,str_info,strmap_info)

      case(UNITY)

        idx = idx_mel_list(rule%labels(1),op_info)
        if(idx.lt.0)
     &       call quit(1,'process_me_lists','Label not on list: "'//
     &       trim(rule%labels(1))//'"')

        mel_pnt => op_info%mel_arr(idx)%mel
        call dens_parameters(+1,rule%parameters,
     &               minblk,maxblk,imode)

        if (imode.ne.1) then
          call zeroop(mel_pnt)
          minblk = 1
          maxblk = mel_pnt%op%n_occ_cls
        end if

        do idx = minblk, maxblk
          call add_unity(1d0,mel_pnt,idx,orb_info,str_info)
        end do
c dbg
c        write(luout,*)'writing unity'
c        call wrt_mel_file(luout,5,mel_pnt,1,1,
c     &       str_info,orb_info)
c dbg

      case(ASSIGN_ME2OP)
        
        if (rule%n_labels.lt.2)
     &     call quit(1,'process_me_lists','two labels expected for '
     &       //trim(ASSIGN_ME2OP))

        call assign_me_list(rule%labels(1),rule%labels(2),op_info)

      case(RES_ME_LIST)
        
        if (rule%n_labels.ne.1)
     &     call quit(1,'process_me_lists','one label expected for '
     &       //trim(RES_ME_LIST))

        call reset_me_list(rule%labels(1),op_info)

      case(DELETE_ME_LIST)
        
        if (rule%n_labels.ne.1)
     &     call quit(1,'process_me_lists','one label expected for '
     &       //trim(DELETE_ME_LIST))

        call del_me_list(rule%labels(1),op_info)

      case(IMPORT)

        if (form_test) exit loop
        call import_parameters(+1,rule%parameters,
     &       list_type,env_type)

        call import_op_el(rule%labels(1),
     &       list_type,env_type,
     &       op_info,str_info,strmap_info,orb_info)

      case(PRECONDITIONER)

        if (form_test) exit loop
        if (rule%n_labels.lt.2)
     &     call quit(1,'process_me_lists',
     &       'at least two labels expected for '
     &       //trim(PRECONDITIONER))

        ! quick fix:
        if (rule%n_labels.eq.2) then
          mode(1:len(mode)) = ' '
          mode = 'dia-F'
          ! very quick and dirty now:
          if (rule%n_parameter_strings.eq.2) mode(5:5) = 'H'
          if (rule%n_parameter_strings.eq.3) mode(5:8) = 'F+id'
        else
          mode(1:len(mode)) = ' '
          mode = 'dia-R12'
        end if

        call set_prc4op(rule%labels(1),mode,
     &       rule%labels(2:),rule%n_labels-1,
     &       op_info,str_info,orb_info)

      case(EVAL)

        if (form_test) exit loop
        call evaluate(rule%labels(1),
     &       op_info,form_info,str_info,strmap_info,orb_info)

      case(INVERT)

        if (rule%n_labels.ne.2)
     &     call quit(1,'process_me_lists','two labels expected for '
     &       //trim(INVERT))
        call form_parameters(+1,rule%parameters,
     &       rule%n_parameter_strings,title,imode,mode)

        call inv_op(rule%labels(2),rule%labels(1),mode,
     &       op_info,orb_info,str_info)

      case(ADD)

        call add_parameters(+1,rule%parameters,nfac,fac,maxfac)

        call add_op(rule%labels(1),fac,rule%labels(2),nfac,
     &       op_info,orb_info,str_info)

      case(SCALE)

        call scale_parameters(+1,rule%parameters,imode,
     &       nblk,idxblk,fac,maxfac)

        call scale_op(rule%labels(1),
     &       imode,idxblk,fac,rule%labels(2:),nblk,
     &       op_info,orb_info,str_info)

      case(EVALPROP)

        call evalprop_parameters(+1,rule%parameters,ndens,rank,env_type)

        call prop_evaluate(ndens,rank,rule%labels,
     &       env_type,op_info,str_info,orb_info)

      case(SOLVENLEQ)

        if (form_test) exit loop
        call solve_parameters(+1,rule%parameters,
     &       rule%n_parameter_strings,
     &       nopt,nroots,mode)

        if (rule%n_labels.lt.3*nopt+2)
     &       call quit(1,'process_me_lists',
     &       'too few labels to be passed for '//
     &       trim(SOLVENLEQ))

        nspecial = rule%n_labels-3*nopt-2
        ioff = 1
        if (nspecial.gt.0) ioff=2

        call solve_nleq(mode,nopt,
     &       rule%labels(1:nopt),               ! to be opt.
     &       rule%labels(nopt+1:nopt+nopt),     ! residual
     &       rule%labels(2*nopt+1:2*nopt+nopt), ! precond.
     &       rule%labels(3*nopt+1),             ! energy
     &       rule%labels(3*nopt+2),             ! formula
     &       rule%labels(3*nopt+ioff+1:
     &                   3*nopt+ioff+nspecial),
     &          nspecial,                       ! specials
     &       op_info,form_info,str_info,strmap_info,orb_info)

      case(SOLVELEQ)
        if (form_test) exit loop
        call solve_parameters(+1,rule%parameters,
     &       rule%n_parameter_strings,
     &       nopt,nroots,mode)

        if (rule%n_labels.lt.5*nopt+1)
     &       call quit(1,'process_me_lists',
     &       'incorrect number of labels to be passed for '//
     &       trim(SOLVELEQ))

        nspecial = rule%n_labels-(5*nopt+1)
        ioff = 0
        if (nspecial.gt.0) ioff=1

        call solve_leq(mode,nopt,nroots,
     &       rule%labels(1:nopt),               ! to be opt.
     &       rule%labels(  nopt+1:  nopt+nopt), ! precond.
     &       rule%labels(2*nopt+1:2*nopt+nopt), ! mvp-labels
     &       rule%labels(3*nopt+1:3*nopt+nopt), ! metric-labels
     &       rule%labels(4*nopt+1:4*nopt+nopt), ! rhs-labels
     &       rule%labels(5*nopt+1),             ! formula
     &       rule%labels(5*nopt+ioff+1:
     &                   5*nopt+ioff+nspecial),nspecial,
     &       op_info,form_info,str_info,strmap_info,orb_info)

      case(SOLVEEVP)
        if (form_test) exit loop
        call solve_parameters(+1,rule%parameters,
     &       rule%n_parameter_strings,
     &       nopt,nroots,mode)

        if (rule%n_labels.lt.4*nopt+1)
     &       call quit(1,'process_me_lists',
     &       'incorrect number of labels to be passed for '//
     &       trim(SOLVEEVP))
        
        nspecial = rule%n_labels-(4*nopt+1)
        ioff = 0
        if (nspecial.gt.0) ioff=1

        call solve_evp(mode,nopt,nroots,
     &       rule%labels(1:nopt),               ! to be opt.
     &       rule%labels(  nopt+1:  nopt+nopt), ! precond.
     &       rule%labels(2*nopt+1:2*nopt+nopt), ! mvp-labels
     &       rule%labels(3*nopt+1:3*nopt+nopt), ! metric-labels
     &       rule%labels(4*nopt+1),             ! formula
     &       rule%labels(4*nopt+ioff+1:
     &                   4*nopt+ioff+nspecial),nspecial,
     &       op_info,form_info,str_info,strmap_info,orb_info)

      case(SET_FREQ)

        idx = idx_mel_list(rule%labels(1),op_info)
        if(idx.lt.0)
     &       call quit(1,'process_me_lists','Label not on list: "'//
     &       trim(rule%labels(1))//'"')
        mel_pnt => op_info%mel_arr(idx)%mel
        call freq_parameters(+1,rule%parameters,freq)
        call set_frequency(mel_pnt,freq)

      case(PRINT_RES)
        if (form_test) exit loop
        idx = idx_mel_list(rule%labels(1),op_info)
        if(idx.lt.0)
     &       call quit(1,'process_me_lists','Label not on list: "'//
     &       trim(rule%labels(1))//'"')
        mel_pnt => op_info%mel_arr(idx)%mel
        allocate(ifreq_temp(maximum_order))
        call ord_parameters(+1,rule%parameters,
     &                      order,dummy,ifreq_temp)
        allocate(ifreq(order))
        ifreq = ifreq_temp(1:order)
        deallocate(ifreq_temp)
        call print_result(order,ifreq,mel_pnt,.false.,orb_info)
        deallocate(ifreq)

      case(PRINT_MEL)
        if (form_test) exit loop
        idx = idx_mel_list(rule%labels(1),op_info)
        mel_pnt => op_info%mel_arr(idx)%mel
        call form_parameters(+1,rule%parameters,
     &       rule%n_parameter_strings,title,imode,mode)
        call print_list(title,mel_pnt,mode,orb_info,str_info)

      case(SET_MEL)
        if (form_test) exit loop
        idx = idx_mel_list(rule%labels(1),op_info)
        if(idx.lt.0)
     &       call quit(1,'process_me_lists','Label not on list: "'//
     &       trim(rule%labels(1))//'"')

        mel_pnt => op_info%mel_arr(idx)%mel

        call scale_parameters(+1,rule%parameters,imode,
     &       nblk,idxblk,fac,maxfac)
        call zeroop(mel_pnt)
c        call diag_guess(mel_pnt,fac,idxblk,nblk,1,imode,
c     &       op_info,str_info,strmap_info,orb_info)
        call set_list(mel_pnt,idxblk,fac,nblk)

      case(EXTRACT_DIAG)
        if (form_test) exit loop
        if (rule%n_labels.ne.2)
     &     call quit(1,'process_me_lists',
     &       'two labels expected for '
     &       //trim(EXTRACT_DIAG))

        call dia_from_op(rule%labels(1),rule%labels(2),
     &       op_info,str_info,orb_info)

      case(REORDER_MEL)
        if (form_test) exit loop
        if (rule%n_labels.ne.2)
     &     call quit(1,'process_me_lists',
     &       'two labels expected for '
     &       //trim(REORDER_MEL))
        call form_parameters(+1,rule%parameters,
     &       rule%n_parameter_strings,title,imode,mode)

        call reo_mel(rule%labels(1),rule%labels(2),
     &       op_info,str_info,strmap_info,orb_info,imode)

      case default
        call quit(1,'process_me_lists','unknown command: '//
     &       trim(rule%command))
      end select

      end do loop

      return
      end
