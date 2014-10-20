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
     &     maxterms = 20, maximum_order = 10
      integer ::
     &     idxterms(maxterms), idx_sv(maxterms),
     &     iblkmin(maxterms), iblkmax(maxterms),
     &     connect(maxterms*2), avoid(maxterms*2),
     &     inproj(maxterms*2), nconnect, navoid, ninproj,
     &     idx, jdx, ioff, ncat, nint, nrename, nop,
     &     ansatz, ipos, idum, level, order, nterms, mode, nint2,
     &     ninclude, ninclude_or, nexclude,
     &     iblk_include(maxterms), iblk_include_or(maxterms),
     &     iblk_exclude(maxterms)
      type(formula), pointer ::
     &     form_pnt, form0_pnt
      character(len_command_par) ::
     &     title, strdum, approx, typ_str, mode_str
      character(len=512) ::
     &     form_str

      integer ::
     &     idx_formlist
      integer, allocatable ::
     &     idxfreqdum(:), idxfreq(:), pop_idxdum(:), pop_idx(:)

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

      case(DEF_ECC_LAGRANGIAN)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,strdum)
        ioff = rule%n_update
        call set_ecc_lagrangian(form_pnt,
     &       title,rule%labels(ioff+1),rule%n_labels-1,ansatz,typ_str,
     &       op_info,orb_info)

      case(DEF_CCPT_LAGRANGIAN)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,ansatz,typ_str)
        ioff = rule%n_update
        call set_ccpt_lagrangian(form_pnt,
     &       title,rule%labels(ioff+1),rule%n_labels-1,ansatz,typ_str,
     &       op_info,orb_info)

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
     &       title,nint,strdum)
        ioff = rule%n_update
        if (rule%n_parameter_strings.le.1)  nint = 1
        call set_hhat2(form_pnt,
     &       title,rule%labels(ioff+1),
     &             rule%labels(ioff+2),rule%labels(ioff+3),
     &       nint,op_info)
      case(DEF_R12INTM_FORMAL)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,typ_str)
        ioff = rule%n_update
        call set_r12intm_formal4(form_pnt,
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
     &       op_info,orb_info,form_info)
      case(DEF_CCR12_METRIC)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,ansatz,strdum)
        ioff = rule%n_update
        call set_r12_metric(form_pnt,
     &       title,rule%labels(ioff+1),rule%n_labels-ioff,ansatz,
     &       op_info,orb_info)
      case(SPLIT_R12EXC_FORMULA)  
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,mode,strdum)
        jdx = idx_formlist(trim(rule%labels(2)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_r12exc_split(form_pnt,form0_pnt,
     &       rule%labels(3),
     &       mode,
     &       rule%labels(4:),rule%n_labels-3,
     &       op_info)
      case(DEF_EXP_FORMULA)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,ansatz,approx)
        ioff = rule%n_update
        call set_experimental_formula(form_pnt,
     &       title,rule%labels(ioff+1),rule%n_labels-ioff,
     &       ansatz,approx,
     &       op_info,orb_info)
      case(DEF_FORMULA)
c        call quit(1,DEF_FORMULA,'not yet')
        call def_form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       form_str,title)
        ioff = rule%n_update
        call set_formula(form_pnt,
     &       form_str,title,
     &       op_info)
      case(EXPAND_OP_PRODUCT)
        call expand_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nop,idx_sv,iblkmin,iblkmax,
     &       connect,nconnect,
     &       avoid,navoid,
     &       inproj,ninproj)
        ioff = rule%n_update
        call form_expand_op_product(.true.,form_pnt,1d0,
     &       title,rule%labels(ioff+1),rule%labels(ioff+2),nop,
     &       idx_sv,iblkmin,iblkmax,
     &       connect,nconnect,
     &       avoid,navoid,
     &       inproj,ninproj,
     &       strdum,0,
     &       .false.,op_info,orb_info)
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
      case(SUM_HERMIT)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nint,strdum)
        ioff = rule%n_update
        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_sum_hermite(form_pnt,form0_pnt,
     &       title,rule%labels(ioff+2),
     &       op_info)
      case(EXPAND)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nint,strdum)
        ioff = rule%n_update
        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_expand_subexpr(form_pnt,form0_pnt,
     &       title,0,
     &       nint,rule%labels(ioff+2),
     &       op_info,form_info
     &       )
      case(REPLACE)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nrename,strdum)
        ioff = rule%n_update
        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        if (jdx.le.0)
     &       call quit(1,'process_formulae',
     &       'label not found: '//trim(rule%labels(ioff+1)))
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
     &       .false.,op_info
     &       )
      case(DERIVATIVE)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nop,strdum)
        ioff = rule%n_update

        ! fix:
        if (rule%n_parameter_strings.eq.1) nop = 1
        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_deriv2(form_pnt,form0_pnt,
     &       title,nop,
     &       rule%labels(ioff+2),
     &       rule%labels(ioff+3),
     &       rule%labels(ioff+3+nop),
     &       op_info)
      case(LEQ_SPLIT)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nop,strdum)
        if (rule%n_parameter_strings.eq.1) nop = 1
        call leq_post(rule%labels(1),rule%labels(2),rule%labels(3),
     &       rule%labels(4),
     &       rule%labels(5),
     &       rule%labels(6:),nop,
     &       title,title,
     &       op_info,form_info)
      case(OPTIMIZE)
        ioff = 1
        call opt_parameters(+1,
     &       rule%parameters,ncat,nint)
c dbg fix by mh
        if (ioff+ncat+1.le.rule%n_labels) then
c dbg original
        call form_opt(form_pnt,
     &       ncat,rule%labels(ioff+1),
     &       nint,rule%labels(ioff+ncat+1),
     &       form_info,op_info,str_info,orb_info)
c dbg resume fix
        else
        call form_opt(form_pnt,
     &       ncat,rule%labels(ioff+1),
     &       nint,rule%labels(ioff+ncat),
     &       form_info,op_info,str_info,orb_info)
        end if
c dbg end fix
      case(PRINT_FORMULA)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,strdum)
        call print_formula_drv(form_pnt,title,'LONG',op_info)
      case(TEX_FORMULA)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,strdum)
        call tex_formula_drv(form_pnt,title,op_info)
      case(SELECT_TERMS)
        call select_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       ninclude,ninclude_or,nexclude,
     &       iblk_include,iblk_include_or,iblk_exclude)
        ioff = rule%n_update        
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)        
        form0_pnt => form_info%form_arr(jdx)%form
        call form_select_terms(form_pnt,form0_pnt,
     &       rule%labels(ioff+2:),
     &       ninclude,rule%labels(ioff+3:),iblk_include,
     &       ninclude_or,rule%labels(ioff+ninclude+3:),iblk_include_or,
     &       nexclude,rule%labels(ioff+ninclude+ninclude_or+3:),
     &                                              iblk_exclude,
     &       op_info)
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
      case(EXTRACT_ORDER)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nint,strdum)
        ioff = rule%n_update
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)
        form0_pnt => form_info%form_arr(jdx)%form
        call form_extract_order(form_pnt,form0_pnt,
     &       title, rule%labels(3), nint, op_info)
      case(EXTRACT_FREQ)
        allocate(idxfreqdum(maximum_order),pop_idxdum(100))
        call form_parameters2(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nint,idxfreqdum,nint2,pop_idxdum)
        allocate(idxfreq(nint),pop_idx(nint2))
        idxfreq = idxfreqdum(1:nint)
        pop_idx = pop_idxdum(1:nint2)
        deallocate(idxfreqdum,pop_idxdum)
        ioff = rule%n_update
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)
        form0_pnt => form_info%form_arr(jdx)%form
        call form_extract_freq(form_pnt,form0_pnt,
     &       title, rule%labels(3), nint, idxfreq, nint2, pop_idx,
     &       op_info)
        deallocate(idxfreq,pop_idx)
      case(CLASS_FORMULA)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idum,strdum)
        call class_formula_drv(form_pnt,title,op_info,idum)
      case(SELECT_HERMIT)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,nint,strdum)
        ioff = rule%n_update
        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)
        form0_pnt => form_info%form_arr(jdx)%form
        call form_select_hermitian(form_pnt,form0_pnt,
     &       title,rule%labels(ioff+2),
     &       op_info)
      case(SELECT_LINE)
        call form_parameters(+1,
     &       rule%parameters,rule%n_parameter_strings,
     &       title,idx,mode_str)
        ioff = rule%n_update

        jdx = idx_formlist(trim(rule%labels(ioff+1)),form_info)
        form0_pnt => form_info%form_arr(jdx)%form
        call form_select_line(form_pnt,form0_pnt,
     &       title,rule%labels(ioff+2),
     &       rule%n_labels-ioff-2,rule%labels(ioff+3),
     &       idx,mode_str,
     &       op_info
     &       )
      case default
        call quit(1,'process_formulae','unknown command: '//
     &       trim(rule%command))
      end select
      return
      end
