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
     &     maxfac = 20
      real(8) ::
     &     fac(maxfac)
      integer ::
     &     idxblk(maxfac),
     &     idx, jdx, ioff, nfac, nblk, nspecial, imode,
     &     absym,casym,gamma,s2,ms,nopt,nroots,ndens,rank
      logical ::
     &     ms_fix
      type(me_list), pointer ::
     &     mel_pnt
      character(len_command_par) ::
     &     title, env_type, list_type, mode
      character(len_command_par), allocatable ::
     &     label_met(:)

      integer ::
     &     idx_formlist, order, dummy
      integer, external ::
     &     idx_mel_list

      if (rule%type.ne.ttype_opme)
     &     call quit(1,'process_me_lists',
     &               'called for wrong target type')
      if (rule%n_labels.lt.1)
     &     call quit(1,'process_me_lists',
     &               'at least one label expected')

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

        call zeroop(mel_pnt)

        call add_unity(1d0,mel_pnt,1,orb_info)
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

      case(DELETE_ME_LIST)
        
        if (rule%n_labels.ne.1)
     &     call quit(1,'process_me_lists','one label expected for '
     &       //trim(DELETE_ME_LIST))

        call del_me_list(rule%labels(1),op_info)

      case(IMPORT)

        call import_parameters(+1,rule%parameters,
     &       list_type,env_type)

        call import_op_el(rule%labels(1),
     &       list_type,env_type,
     &       op_info,str_info,orb_info)

      case(PRECONDITIONER)

        if (rule%n_labels.lt.2)
     &     call quit(1,'process_me_lists',
     &       'at least two labels expected for '
     &       //trim(PRECONDITIONER))

        ! quick fix:
        if (rule%n_labels.eq.2) then
          mode(1:len(mode)) = ' '
          mode = 'dia-F'
        else
          mode(1:len(mode)) = ' '
          mode = 'dia-R12'
        end if

        call set_prc4op(rule%labels(1),mode,
     &       rule%labels(2:),rule%n_labels-1,
     &       op_info,str_info,orb_info)

      case(EVAL)

        call evaluate(rule%labels(1),
     &       op_info,form_info,str_info,strmap_info,orb_info)

      case(INVERT)

        if (rule%n_labels.ne.2)
     &     call quit(1,'process_me_lists','two labels expected for '
     &       //trim(INVERT))

        call inv_op(rule%labels(2),rule%labels(1),
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
        call solve_parameters(+1,rule%parameters,
     &       rule%n_parameter_strings,
     &       nopt,nroots,mode)

        if (rule%n_labels.ne.4*nopt+1)
     &       call quit(1,'process_me_lists',
     &       'incorrect number of labels to be passed for '//
     &       trim(SOLVELEQ))

        call solve_leq(mode,nopt,nroots,
     &       rule%labels(1:nopt),               ! to be opt.
     &       rule%labels(  nopt+1:  nopt+nopt), ! precond.
     &       rule%labels(2*nopt+1:2*nopt+nopt), ! mvp-labels
     &       rule%labels(3*nopt+1:3*nopt+nopt), ! rhs-labels
     &       rule%labels(4*nopt+1),             ! formula
     &       op_info,form_info,str_info,strmap_info,orb_info)

      case(SOLVEEVP)
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

        call opt_parameters(+1,rule%parameters,
     &                      order,dummy)
        call set_frequency(mel_pnt,order)

      case default
        call quit(1,'process_me_lists','unknown command: '//
     &       trim(rule%command))
      end select

      return
      end
