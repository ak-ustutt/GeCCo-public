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

      integer ::
     &     idx, jdx, ioff,
     &     absym,casym,gamma,s2,ms,nopt,nroots,ndens,rank
      type(me_list), pointer ::
     &     mel_pnt
      character(len_command_par) ::
     &     title, env_type

      integer ::
     &     idx_formlist

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
     &       absym,casym,gamma,s2,ms)

        call define_me_list(rule%labels(1),rule%labels(2),
     &       absym,casym,gamma,s2,ms,
     &       -1,-1,
     &       op_info,orb_info,str_info,strmap_info)

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
     &       env_type)

        call import_op_el(rule%labels(1),
     &       op_info,
     &       env_type,str_info,orb_info)

      case(PRECONDITIONER)

        call set_prc4op(rule%labels(1),rule%labels(2),
     &       op_info,str_info,orb_info)

      case(EVAL)

        call evaluate(rule%labels(1),
     &       op_info,form_info,str_info,strmap_info,orb_info)

      case(EVALPROP)

c dbg
        print *,'lebenszeichen ',trim(rule%parameters(1))
c dbg
        call evalprop_parameters(+1,rule%parameters,ndens,rank,env_type)

c dbg
        print *,':>',ndens,rank,trim(env_type)
c dbg
        call prop_evaluate(ndens,rank,rule%labels,
     &       env_type,op_info,str_info,orb_info)

      case(SOLVENLEQ)

        call solve_parameters(+1,rule%parameters,
     &       nopt,nroots)

        if (rule%n_labels.ne.3*nopt+2)
     &       call quit(1,'process_me_lists',
     &       'incorrect number of labels to be passed for '//
     &       trim(SOLVENLEQ))

        call solve_nleq(nopt,
     &       rule%labels(1:nopt),               ! to be opt.
     &       rule%labels(nopt+1:nopt+nopt),     ! residual
     &       rule%labels(2*nopt+1:2*nopt+nopt), ! precond.
     &       rule%labels(3*nopt+1),             ! energy
     &       rule%labels(3*nopt+2),             ! formula
     &       op_info,form_info,str_info,strmap_info,orb_info)

      case(SOLVELEQ)
        call solve_parameters(+1,rule%parameters,
     &       nopt,nroots)

        if (rule%n_labels.ne.4*nopt+1)
     &       call quit(1,'process_me_lists',
     &       'incorrect number of labels to be passed for '//
     &       trim(SOLVELEQ))

        call solve_leq(nopt,nroots,
     &       rule%labels(1:nopt),               ! to be opt.
     &       rule%labels(  nopt+1:  nopt+nopt), ! precond.
     &       rule%labels(2*nopt+1:2*nopt+nopt), ! mvp-labels
     &       rule%labels(3*nopt+1:3*nopt+nopt), ! rhs-labels
     &       rule%labels(4*nopt+1),             ! formula
     &       op_info,form_info,str_info,strmap_info,orb_info)

      case(SOLVEEVP)
        call solve_parameters(+1,rule%parameters,
     &       nopt,nroots)

        if (rule%n_labels.ne.3*nopt+1)
     &       call quit(1,'process_me_lists',
     &       'incorrect number of labels to be passed for '//
     &       trim(SOLVEEVP))

        call solve_evp(nopt,nroots,
     &       rule%labels(1:nopt),               ! to be opt.
     &       rule%labels(  nopt+1:  nopt+nopt), ! precond.
     &       rule%labels(2*nopt+1:2*nopt+nopt), ! mvp-labels
     &       rule%labels(3*nopt+1),             ! formula
     &       op_info,form_info,str_info,strmap_info,orb_info)

      case default
        call quit(1,'process_me_lists','unknown command: '//
     &       trim(rule%command))
      end select

      return
      end
