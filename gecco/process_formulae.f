*----------------------------------------------------------------------*
      subroutine process_formulae(rule,
     &                            form_info,op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     process formula rules
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_target.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'par_actions.h'

      type(action), intent(in) ::
     &     rule
      type(formula_info), intent(inout) ::
     &     form_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer, parameter ::
     &     maxterms = 20
      integer ::
     &     idxterms(maxterms), idx_sv(maxterms),
     &     idx, jdx, ioff, ncat, nint, nrename, nop,
     &     ansatz, ipos, idum, level, nterms
      type(formula), pointer ::
     &     form_pnt, form0_pnt
      character(len_command_par) ::
     &     title, strdum, approx, typ_str

      integer ::
     &     idx_formlist

      if (rule%type.ne.ttype_frm)
     &     call quit(1,'process_formulae',
     &     'called for wrong target type')
      if (rule%n_labels.lt.1)
     &     call quit(1,'process_formulae','at least one label expected')

      ! allocate a new entry, if necessary
      do idx = 1, rule%n_update
        jdx = idx_formlist(trim(rule%labels(idx)),form_info)
        if (jdx.gt.0) cycle
        call add_formula(form_info,trim(rule%labels(idx)))
      end do

      ! point to formula referenced by first label:
      idx = idx_formlist(trim(rule%labels(1)),form_info)
      form_pnt => form_info%form_arr(idx)%form

      select case(trim(rule%command))
      case(CHECK_FORMGEN)
        call check_formula_generators(form_pnt,op_info,orb_info)
      case(DEF_CC_LAGRANGIAN)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,strdum)
        ioff = rule%n_update
        call set_cc_lagrangian2(form_pnt,
     &       title,rule%labels(ioff+1),rule%labels(ioff+2),
     &             rule%labels(ioff+3),rule%labels(ioff+4),
     &       op_info)

      case(DEF_CC_HBAR)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,strdum)
        ioff = rule%n_update
        call set_cc_hbar_formula(form_pnt,
     &       title,rule%labels(ioff+1),
     &             rule%labels(ioff+2),rule%labels(ioff+3),
     &       op_info)
      case(DEF_HHAT)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,strdum)
        ioff = rule%n_update
        call set_hhat2(form_pnt,
     &       title,rule%labels(ioff+1),
     &             rule%labels(ioff+2),rule%labels(ioff+3),
     &       op_info)
      case(DEF_R12INTM_FORMAL)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,typ_str)
        ioff = rule%n_update
        call set_r12intm_formal3(form_pnt,
     &         title,rule%labels(ioff+1),rule%labels(ioff+2),
     &         rule%n_labels-ioff-1,typ_str,
     &         op_info,orb_info)

      case(DEF_R12INTM_CABS)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,ansatz,approx)
        ioff = rule%n_update
        call set_r12intm_cabs3(form_pnt,
     &       title,rule%labels(ioff+1),rule%n_labels-ioff,
     &       approx(1:2),ansatz,approx(3:),
     &       op_info,orb_info)

      case(DEF_MPR12_LAGRANGIAN)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,level,strdum)
c prelim        
        if (level.ne.2) call quit(1,'process_formulae',
     &       'MP: only level==2 implemented')
        ioff = rule%n_update
        call set_mp2_r12_lagrangian(form_pnt,
     &       title,rule%labels(ioff+1),rule%n_labels-ioff,
     &       op_info,orb_info)
c prelim
      case(DEF_CCR12_LAGRANGIAN)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,ansatz,strdum)
        ioff = rule%n_update
        call set_r12_lagrangian(form_pnt,
     &       title,rule%labels(ioff+1),rule%n_labels-ioff,ansatz,
     &       op_info,orb_info)
      case(EXPAND_OP_PRODUCT)
        call form_parameters2(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nop,idx_sv)
        ioff = rule%n_update
c dbg
        print *,'ioff,nop = ',ioff,nop
c dbg        
        call form_expand_op_product(form_pnt,
     &       title,rule%labels(ioff+1),rule%labels(ioff+2),nop,idx_sv,
     &       op_info,orb_info)
      case(FACTOR_OUT)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nint,strdum)
        ioff = rule%n_update
        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_factor_out(form_pnt,form0_pnt,
     &       title,
     &       nint,rule%labels(ioff+2),
     &       op_info,form_info
     &       )
      case(REPLACE)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nrename,strdum)
        ioff = rule%n_update
        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_op_replace_drv(form_pnt,form0_pnt,
     &       title,
     &       nrename,rule%labels(ioff+2),
     &       op_info
     &       )
      case(INVARIANT)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,strdum)
        ioff = rule%n_update
        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_invariant(form_pnt,form0_pnt,
     &       title,rule%labels(ioff+2),
     &       rule%n_labels-ioff-2,rule%labels(ioff+3),
     &       op_info
     &       )
      case(DERIVATIVE)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,strdum)
        ioff = rule%n_update
        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_deriv2(form_pnt,form0_pnt,
     &       title,1,
     &             rule%labels(ioff+2),
     &             rule%labels(ioff+3),
     &             rule%labels(ioff+4),
     &       op_info)
      case(LEQ_SPLIT)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,strdum)
        call leq_post(rule%labels(1),rule%labels(2),rule%labels(3),
     &       rule%labels(4),
     &       rule%labels(5),
     &       rule%labels(6),1,
     &       title,title,
     &       op_info,form_info)
      case(OPTIMIZE)
        ioff = 1
        call opt_parameters(+1,
     &       rule%parameters,ncat,nint)

        call form_opt(form_pnt,
     &       ncat,rule%labels(ioff+1),
     &       nint,rule%labels(ioff+ncat+1),
     &       form_info,op_info,str_info,orb_info)
c      case(CONTRACT)
c        call form_parameters(+1,
c     &       rule%parameters,rule%n_parameter_strings,
c     &       title,idum,strdum)
c        ioff = rule%n_update+1
c        jdx = rule%n_labels-ioff
c        call test_contract(form_pnt,title,
c     &       rule%labels(2),
c     &       rule%labels(ioff+1:ioff+jdx),jdx,
c     &       op_info)
c      case(EXTRACT_TERM)
      case(DEL_TERMS)
        call modify_parameters(+1,
     &       rule%parameters,nterms,idxterms,maxterms)
        ioff = rule%n_update        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_del_terms(form_pnt,form0_pnt,
     &       nterms,idxterms,-1,
     &       op_info)
      case(KEEP_TERMS)
        call modify_parameters(+1,
     &       rule%parameters,nterms,idxterms,maxterms)
        ioff = rule%n_update        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_del_terms(form_pnt,form0_pnt,
     &       nterms,idxterms,+1,
     &       op_info)
      case(MODIFY_FACTORIZATION)
        call modify_parameters(+1,
     &       rule%parameters,nterms,idxterms,maxterms)
        ioff = rule%n_update        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_mod_factorization(form_pnt,form0_pnt,
     &       nterms,idxterms,
     &       op_info)
      case default
        call quit(1,'process_formulae','unknown command: '//
     &       trim(rule%command))
      end select

      return
      end
