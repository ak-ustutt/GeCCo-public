*----------------------------------------------------------------------*
      subroutine process_operators(rule,op_info,orb_info)
*----------------------------------------------------------------------*
*     process the operator rules
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_target.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'par_actions.h'

      type(action), intent(in) ::
     &     rule
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(operator), pointer ::
     &     op_pnt
      integer ::
     &     idx, idx_t,
     &     ncadiff, min_rank, max_rank, iformal
      character*(len_opname) ::
     &     name_template
      logical ::
     &     dagger

      integer, external ::
     &     idx_oplist2

      if (rule%type.ne.ttype_op)
     &     call quit(1,'process_operators',
     &     'called for wrong target type')
      if (rule%n_labels.ne.1)
     &     call quit(1,'process_operators','exactly one label expected')

      ! allocate a new entry
      call add_operator(trim(rule%labels(1)),op_info)
      idx = idx_oplist2(trim(rule%labels(1)),op_info)
      op_pnt => op_info%op_arr(idx)%op

      select case(trim(rule%command))
      case(DEF_GENERAL_OPERATOR)
        call quit(1,'process_operators',
     &       'not yet ready: '//DEF_GENERAL_OPERATOR)
      case(DEF_SCALAR)
        call set_hop(op_pnt,trim(rule%labels(1)),.false.,
     &       0,0,1,orb_info)
      case(DEF_HAMILTONIAN)
        if (rule%n_parameter_strings.lt.1)
     &       call quit(1,'process_operatorss',
     &       'no parameters provided for '//DEF_HAMILTONIAN)
        call hop_parameters(+1,rule%parameters,
     &                      min_rank,max_rank,iformal)
        call set_hop(op_pnt,trim(rule%labels(1)),.false.,
     &       min_rank,max_rank,iformal,orb_info)
      case(DEF_EXCITATION)
        if (rule%n_parameter_strings.lt.1)
     &       call quit(1,'process_operators',
     &       'no parameters provided for '//DEF_EXCITATION)
        call xop_parameters(+1,rule%parameters,
     &                      dagger,min_rank,max_rank,ncadiff,iformal)
        call set_xop(op_pnt,trim(rule%labels(1)),dagger,
     &       min_rank,max_rank,ncadiff,iformal,orb_info)
      case(DEF_DENSITY)
        if (rule%n_parameter_strings.lt.1)
     &       call quit(1,'process_operators',
     &       'no parameters provided for '//DEF_DENSITY)
        call dens_parameters(+1,rule%parameters,
     &                      min_rank,max_rank,iformal)
        call set_dens(op_pnt,trim(rule%labels(1)),.false.,
     &       min_rank,max_rank,iformal,orb_info)
      case(CLONE_OP)
        if (rule%n_parameter_strings.lt.1)
     &       call quit(1,'process_operators',
     &       'no parameters provided for '//CLONE_OP)
        call cloneop_parameters(+1,rule%parameters,name_template,dagger)
        idx_t = idx_oplist2(name_template,op_info)
        if (idx_t.le.0) then
          call quit(0,CLONE_OP,
     &         'template not defined: '//name_template)
        end if
        call clone_operator(op_pnt,op_info%op_arr(idx_t)%op,orb_info)
        op_pnt%dagger = op_pnt%dagger.xor.dagger
      case default
        call quit(1,'process_operators','unknown command: '//
     &       trim(rule%command))
      end select

      return
      end
