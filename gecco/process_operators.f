*----------------------------------------------------------------------*
      subroutine process_operators(rule,op_info,orb_info)
*----------------------------------------------------------------------*
*     process the operator rules
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
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
      integer, parameter ::
     &     ndef_max = 20
      integer ::
     &     idx, idx_t,
     &     ncadiff, min_rank, max_rank, iformal, ansatz,
     &     ndim1, ndim2, ndef, njoined,
     &     hpvx_constr(2,ngastp,2), gas_constr(2,orb_info%ngas,2,2),
     &     occ_def(2,ngastp,ndef_max)
      character*(len_opname) ::
     &     name_template
      logical ::
     &     dagger, explicit, set_p

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
        call genop_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       dagger,min_rank,max_rank,ncadiff,iformal,explicit,
     &       hpvx_constr,ndim1,gas_constr,ndim2)
        call set_genop(op_pnt,trim(rule%labels(1)),optyp_operator,
     &       dagger,
     &       min_rank,max_rank,ncadiff,hpvx_constr,gas_constr,
     &       iformal,orb_info)
      case(DEF_OP_FROM_OCC)
        call op_from_occ_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       occ_def,ndef,njoined,ndef_max)
        call set_uop2(op_pnt,trim(rule%labels(1)),
     &       occ_def,ndef,njoined,orb_info)
      case(DEF_SCALAR)
        call set_hop(op_pnt,trim(rule%labels(1)),.false.,
     &       0,0,1,.false.,orb_info)
      case(DEF_HAMILTONIAN)
        if (rule%n_parameter_strings.lt.1)
     &       call quit(1,'process_operatorss',
     &       'no parameters provided for '//DEF_HAMILTONIAN)
        call hop_parameters(+1,rule%parameters,
     &                      min_rank,max_rank,iformal,explicit)
        call set_hop(op_pnt,trim(rule%labels(1)),.false.,
     &       min_rank,max_rank,iformal,explicit,orb_info)
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
      case(DEF_CC_HBAR_OP)
        if (rule%n_parameter_strings.lt.1)
     &       call quit(1,'process_operators',
     &       'no parameters provided for '//DEF_DENSITY)
        call xop_parameters(+1,rule%parameters,
     &                      dagger,min_rank,max_rank,ncadiff,iformal)
        call set_cc_hbar(op_pnt,trim(rule%labels(1)),
     &       max_rank,orb_info)
      case(DEF_R12GEMINAL)
        if (rule%n_parameter_strings.lt.1)
     &       call quit(1,'process_operators',
     &       'no parameters provided for '//DEF_R12GEMINAL)
        call r12gem_parameters(+1,rule%parameters,
     &                      dagger,min_rank,ansatz)
        max_rank = 2
        call set_r12gem(op_pnt,trim(rule%labels(1)),dagger,
     &       min_rank,max_rank,ansatz,orb_info)        
      case(DEF_R12COEFF)
        if (rule%n_parameter_strings.lt.1)
     &       call quit(1,'process_operators',
     &       'no parameters provided for '//DEF_R12COEFF)
        call r12int_parameters(+1,rule%parameters,
     &                      set_p,min_rank,max_rank,ncadiff,iformal)
        call set_r12c(op_pnt,trim(rule%labels(1)),set_p,
     &       min_rank,max_rank,ncadiff,iformal,orb_info)        
      case(DEF_R12INT)
        if (rule%n_parameter_strings.lt.1)
     &       call quit(1,'process_operators',
     &       'no parameters provided for '//DEF_R12INT)
        call xop_parameters(+1,rule%parameters,
     &                      dagger,min_rank,max_rank,ncadiff,iformal)
        call set_r12i(op_pnt,trim(rule%labels(1)),dagger,
     &       min_rank,max_rank,ncadiff,iformal,orb_info)        
      case(DEF_R12INTERM)
        if (rule%n_parameter_strings.lt.1)
     &       call quit(1,'process_operators',
     &       'no parameters provided for '//DEF_R12INTERM)
        call xop_parameters(+1,rule%parameters,
     &                      dagger,min_rank,max_rank,ncadiff,iformal)
        call set_r12intm(op_pnt,trim(rule%labels(1)),dagger,
     &       min_rank,max_rank,ncadiff,iformal,op_info,orb_info)        
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
        call clone_operator(op_pnt,op_info%op_arr(idx_t)%op,
     &       dagger,orb_info)
c        op_pnt%dagger = op_pnt%dagger.xor.dagger
      case default
        call quit(1,'process_operators','unknown command: '//
     &       trim(rule%command))
      end select
c dbg
      if (op_pnt%dagger)
     &     call quit(1,'process_operators','op%dagger is obsolete!')
c dbg

      return
      end
