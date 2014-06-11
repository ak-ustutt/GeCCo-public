*----------------------------------------------------------------------*
      subroutine process_rule(rule,tgt_info,
     &     form_info,op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     new driver routine for processing rules
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_strmapinf.h'
      include 'par_actions.h'     
      include 'ifc_input.h'
      include 'ifc_targets.h'

      type(action), intent(in) ::
     &     rule
      type(target_info), intent(inout) ::
     &     tgt_info
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
     &     NEW = 0, OLD = 1, ANY = 2

      integer, parameter ::
     &     maxfac = 40, max_occ = 250, max_nj = 20,
     &     max_label = 200, maxterms = 50, max_pops = 100
      real(8) ::
     &     fac(maxfac), freq, xdum
      integer ::
     &     nblk, njoined, min_rank, max_rank, min_xrank, max_xrank,
     &     ncadiff, iformal, n_ap, ansatz, hermitian, iorder, spec,
     &     ninclude, ninclude_or, nexclude, norb, icase, icaseF,
     &     minblk, maxblk, idx, jdx, ioff, nfac, nspecial, imode,
     &     nop, nop2, nint, ncat, level, nconnect, navoid, ninproj,
     &     absym,casym,gamma,s2,ms,nopt,nroots,ndens,rank,nterms,ncmp,
     &     dgam, dms, nspcfrm, ndescr, ntmp, targ_root
      integer ::
     &     idxblk(maxfac), idxterms(maxterms), idx_sv(maxterms),
     &     iblkmin(maxterms), iblkmax(maxterms),
     &     connect(maxterms*2), avoid(maxterms*2),
     &     inproj(maxterms*2), 
     &     iblk_include(maxterms), iblk_include_or(maxterms),
     &     iblk_exclude(maxterms), iRdef(maxterms)
      logical ::
     &     dagger, explicit, ms_fix, form_test, init, arg_there, reo
      integer, pointer ::
     &     occ_def(:,:,:), nact(:), hpvx_constr(:), hpvxca_constr(:),
     &     gas_constr(:,:,:,:,:,:)
      type(operator), pointer ::
     &     op_pnt, op_pnt2
      type(formula), pointer ::
     &     form_pnt, form0_pnt
      type(me_list), pointer ::
     &     mel_pnt
      character(len=512) ::
     &     title, title2, form_str, mode, strscr
      character(len_command_par) ::
     &     env_type, list_type
      character(len_command_par) ::
     &     label, label2, label_list(max_label), descr(max_label)

      integer, allocatable ::
     &     ifreq(:), pop_idx(:) 

      integer, external ::
     &     idx_formlist, idx_mel_list, idx_oplist2

      ! form_test = true skips time consuming steps -> dry run
      call get_argument_value('general','form_test',lval=form_test)

*----------------------------------------------------------------------*
*     branch according to command
*      subsections (for search):
*        OPERATORS, FORMULAE, ME-LISTS, EVALUATE
*----------------------------------------------------------------------*
      select case(trim(rule%command))
*----------------------------------------------------------------------*
*     subsection OPERATORS
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      case(DEF_GENERAL_OPERATOR)
*----------------------------------------------------------------------*
        allocate(nact(2),
     &           hpvx_constr(2*ngastp),hpvxca_constr(2*2*ngastp),
     &           gas_constr(2,orb_info%ngas,2,2,orb_info%nspin,1))
        call get_arg('LABEL', rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),NEW)
        call get_arg('MIN_RANK',rule,tgt_info,val_int=min_rank)
        call get_arg('MAX_RANK',rule,tgt_info,val_int=max_rank)
        call get_arg('MIN_XRANK',rule,tgt_info,val_int=min_xrank)
        call get_arg('MAX_XRANK',rule,tgt_info,val_int=max_xrank)
        call get_arg('CHARGE',rule,tgt_info,val_int=ncadiff)
        call get_arg('FORMAL',rule,tgt_info,val_int=iformal)
        call get_arg('CONSTR_HPVX',rule,tgt_info,
     &       val_int_list=hpvx_constr)
        call get_arg('CONSTR_HPVX_CA',rule,tgt_info,
     &       val_int_list=hpvxca_constr)
        call get_arg('CONSTR_GAS',rule,tgt_info,val_restr=gas_constr)
        call get_arg('CORE',  rule,tgt_info,val_int_list=nact)
        call set_genop2(op_pnt,trim(label),optyp_operator,
     &       min_rank,max_rank,ncadiff,
     &       min_xrank,max_xrank,
     &       hpvx_constr,hpvxca_constr,
     &       gas_constr,iformal,nact,
     &       orb_info)
        deallocate(nact,hpvx_constr,hpvxca_constr,gas_constr)
*----------------------------------------------------------------------*
      case(DEF_OP_FROM_OCC)
*----------------------------------------------------------------------*
        allocate(occ_def(ngastp,2,max_occ),nact(2*max_nj))
        call get_arg('LABEL', rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),NEW)
        call get_arg('JOIN',  rule,tgt_info,val_int=njoined)
        call get_arg('CORE',  rule,tgt_info,val_int_list=nact)
        call get_arg('FORMAL',rule,tgt_info,val_int=iformal)
        call get_arg('DESCR', rule,tgt_info,val_str=strscr,
     &                             success=arg_there)
        if (.not.arg_there) then 
          call get_arg('BLOCKS',rule,tgt_info,val_int=nblk)
          call get_arg('OCC',   rule,tgt_info,val_occ=occ_def)
          call set_uop2(op_pnt,trim(label),
     &       occ_def,nblk,njoined,nact,iformal,orb_info)
        else
          call set_uop3(op_pnt,trim(label),
     &       strscr,njoined,nact,iformal,orb_info)
        end if
        deallocate(occ_def,nact)
*----------------------------------------------------------------------*
      case(DEF_SCALAR)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),NEW)
        call set_hop(op_pnt,trim(label),.false.,
     &       0,0,1,.false.,IEXTR,1,orb_info)
*----------------------------------------------------------------------*
      case(DEF_HAMILTONIAN)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),NEW)
        call get_arg('MIN_RANK',rule,tgt_info,val_int=min_rank)
        call get_arg('MAX_RANK',rule,tgt_info,val_int=max_rank)
        call get_arg('FORMAL',rule,tgt_info,val_int=iformal)
        call get_arg('SET_X',rule,tgt_info,val_log=explicit)
        call get_arg('X_SPCS',rule,tgt_info,
     &       val_int_list=iblk_exclude,ndim=nexclude)
        call set_hop(op_pnt,trim(label),.false.,
     &       min_rank,max_rank,iformal,explicit,
     &       iblk_exclude,nexclude,orb_info)        
*----------------------------------------------------------------------*
      case(DEF_EXCITATION)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),NEW)
        call get_arg('MIN_RANK',rule,tgt_info,val_int=min_rank)
        call get_arg('MAX_RANK',rule,tgt_info,val_int=max_rank)
        call get_arg('CHARGE',rule,tgt_info,val_int=ncadiff)
        call get_arg('ADJOINT',rule,tgt_info,val_log=dagger)
        call get_arg('FORMAL',rule,tgt_info,val_int=iformal)
        call set_xop(op_pnt,trim(label),dagger,
     &       min_rank,max_rank,0,ncadiff,iformal,orb_info)
C*----------------------------------------------------------------------*
C      case(DEF_DENSITY)
C*----------------------------------------------------------------------*
C*----------------------------------------------------------------------*
C      case(DEF_CC_HBAR_OP)
C*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      case(DEF_R12GEMINAL)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),NEW)
        call get_arg('MIN_RANK',rule,tgt_info,val_int=min_rank)
        call get_arg('MAX_RANK',rule,tgt_info,val_int=max_rank)
        call get_arg('ANSATZ',rule,tgt_info,val_int=ansatz)
        call get_arg('N_PART',rule,tgt_info,val_int=n_ap)
        call set_r12gem(op_pnt,trim(label),n_ap,
     &       min_rank,max_rank,ansatz,orb_info)        
*----------------------------------------------------------------------*
      case(DEF_R12COEFF)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),NEW)
        call get_arg('MIN_RANK',rule,tgt_info,val_int=min_rank)
        call get_arg('MAX_RANK',rule,tgt_info,val_int=max_rank)
        call get_arg('CHARGE',rule,tgt_info,val_int=ncadiff)
        call get_arg('ADJOINT',rule,tgt_info,val_log=dagger)
        call get_arg('FORMAL',rule,tgt_info,val_int=iformal)
        call set_r12c(op_pnt,trim(label),dagger,
     &       min_rank,max_rank,ncadiff,iformal,orb_info)        
*----------------------------------------------------------------------*
      case(DEF_R12INT)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),NEW)
        call get_arg('MIN_RANK',rule,tgt_info,val_int=min_rank)
        call get_arg('MAX_RANK',rule,tgt_info,val_int=max_rank)
        call get_arg('CHARGE',rule,tgt_info,val_int=ncadiff)
        call get_arg('N_PART',rule,tgt_info,val_int=n_ap)
        call get_arg('FORMAL',rule,tgt_info,val_int=iformal)
        call set_r12i(op_pnt,trim(label),n_ap,
     &       min_rank,max_rank,ncadiff,iformal,orb_info)        
*----------------------------------------------------------------------*
      case(DEF_R12INTERM)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),NEW)
        call get_arg('MIN_RANK',rule,tgt_info,val_int=min_rank)
        call get_arg('MAX_RANK',rule,tgt_info,val_int=max_rank)
        call get_arg('CHARGE',rule,tgt_info,val_int=ncadiff)
        call get_arg('N_PART',rule,tgt_info,val_int=n_ap)
        call get_arg('ADJOINT',rule,tgt_info,val_log=dagger)
        call get_arg('FORMAL',rule,tgt_info,val_int=iformal)
        call set_r12intm(op_pnt,trim(label),dagger,
     &       min_rank,max_rank,ncadiff,iformal,op_info,orb_info)        
*----------------------------------------------------------------------*
      case(CLONE_OP)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),NEW)
        call get_arg('TEMPLATE',rule,tgt_info,val_label=label2)
        call get_op(op_pnt2,trim(label2),OLD)
        call get_arg('ADJOINT',rule,tgt_info,val_log=dagger)
        call clone_operator(op_pnt,op_pnt2,dagger,orb_info)        
*----------------------------------------------------------------------*
      case(SET_ORDER)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),OLD)
        call get_arg('ORDER',rule,tgt_info,val_int=iorder)
        allocate(ifreq(iorder))
        call get_arg('IDX_FREQ',rule,tgt_info,val_int_list=ifreq)
        call get_arg('SPECIES',rule,tgt_info,val_int=spec)
        call set_pert_order(op_pnt,iorder,spec,ifreq)
        deallocate(ifreq)
*----------------------------------------------------------------------*
      case(SET_HERMIT)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_op(op_pnt,trim(label),OLD)
        call get_arg('CA_SYMMETRY',rule,tgt_info,val_int=hermitian)
        call set_hermitian(op_pnt,hermitian)

*----------------------------------------------------------------------*
*     subsection FORMULAE
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      case(CHECK_FORMGEN)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call check_formula_generators(form_pnt,op_info,orb_info)
*----------------------------------------------------------------------*
      case(DEF_CC_LAGRANGIAN)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OP_RES',rule,tgt_info,val_label=label_list(1))
        call get_arg('OP_H',rule,tgt_info,val_label=label_list(2))
        call get_arg('OP_L',rule,tgt_info,val_label=label_list(3))
        call get_arg('OP_T',rule,tgt_info,val_label=label_list(4))
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_cc_lagrangian2(form_pnt,
     &       title,label_list(1),label_list(2),
     &             label_list(3),label_list(4),
     &       op_info)
*----------------------------------------------------------------------*
      case(DEF_ECC_LAGRANGIAN)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OP_RES',rule,tgt_info,val_label=label_list(1))
        call get_arg('OP_H',rule,tgt_info,val_label=label_list(2))
        call get_arg('OP_L',rule,tgt_info,val_label=label_list(3))
        call get_arg('OP_T',rule,tgt_info,val_label=label_list(4))
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_ecc_lagrangian(form_pnt,
     &       title,label_list,4,0,mode,
     &       op_info,orb_info)
*----------------------------------------------------------------------*
      case(DEF_CCPT_LAGRANGIAN)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list,ndim=nop)
        call get_arg('ANSATZ',rule,tgt_info,val_int=ansatz)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_ccpt_lagrangian(form_pnt,
     &       title,label_list,nop,ansatz,mode,
     &       op_info,orb_info)
*----------------------------------------------------------------------*
      case(DEF_MRCC_LAGRANGIAN)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OP_RES',rule,tgt_info,
     &               val_label=label_list(1))
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list(2:),ndim=nop)
        call get_arg('MAXCOM_RES',rule,tgt_info,val_int=ansatz)
        call get_arg('MAXCOM_EN',rule,tgt_info,val_int=nint)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_mrcc_lagrangian(form_pnt,
     &       title,label_list,nop+1,
     &       ansatz,nint,mode,
     &       op_info,orb_info)
*----------------------------------------------------------------------*
      case(DEF_CC_HBAR)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OP_HBAR',rule,tgt_info,val_label=label_list(1))
        call get_arg('OP_H',rule,tgt_info,val_label=label_list(2))
        call get_arg('OP_T',rule,tgt_info,val_label=label_list(3))
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_cc_hbar_formula(form_pnt,
     &       title,label_list(1),
     &             label_list(2),
     &             label_list(3),
     &       op_info)
*----------------------------------------------------------------------*
      case(DEF_HHAT)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OP_HHAT',rule,tgt_info,val_label=label_list(1))
        call get_arg('OP_H',rule,tgt_info,val_label=label_list(2))
        call get_arg('OP_T',rule,tgt_info,val_label=label_list(3))
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_hhat2(form_pnt,
     &       title,label_list(1),
     &             label_list(2),label_list(3),
     &       nint,op_info)
*----------------------------------------------------------------------*
      case(DEF_R12INTM_FORMAL)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('INTERM',rule,tgt_info,
     &               val_label=label)
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list,ndim=nop)
        call get_arg('ANSATZ',rule,tgt_info,val_int=ansatz)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_r12intm_formal4(form_pnt,
     &         title,label,label_list,
     &         nop,mode,
     &         op_info,orb_info)
*----------------------------------------------------------------------*
      case(DEF_R12INTM_CABS)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('INTERM',rule,tgt_info,
     &               val_label=label_list(1))
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list(2:),ndim=nop)
        call get_arg('ANSATZ',rule,tgt_info,val_int=ansatz)
        call get_arg('APPROX',rule,tgt_info,val_str=mode)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_r12intm_cabs3(form_pnt,
     &       title,label_list,nop+1,
     &       mode(1:2),ansatz,mode(3:),
     &       op_info,orb_info)
*----------------------------------------------------------------------*
      case(DEF_MPR12_LAGRANGIAN)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list,ndim=nop)
        call get_arg('ANSATZ',rule,tgt_info,val_int=ansatz)
        call get_arg('LEVEL',rule,tgt_info,val_int=level)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
c prelim        
        if (level.ne.2) call quit(1,'process_formulae',
     &       'MP: only level==2 implemented')
        call set_mp2_r12_lagrangian(form_pnt,
     &       title,label_list,nop,
     &       op_info,orb_info)
c prelim
*----------------------------------------------------------------------*
      case(DEF_CCR12_LAGRANGIAN)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list,ndim=nop)
        call get_arg('ANSATZ',rule,tgt_info,val_int=ansatz)
c        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_r12_lagrangian(form_pnt,
     &       title,label_list,nop,ansatz,
     &       op_info,orb_info,form_info)
*----------------------------------------------------------------------*
      case(DEF_CCR12_METRIC)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list,ndim=nop)
        call get_arg('ANSATZ',rule,tgt_info,val_int=ansatz)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_r12_metric(form_pnt,
     &       title,label_list,nop,ansatz,
     &       op_info,orb_info)
*----------------------------------------------------------------------*
      case(SPLIT_R12EXC_FORMULA)  
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OP_RES',rule,tgt_info,
     &               val_label=label)
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list,ndim=nop)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call form_r12exc_split(form_pnt,form0_pnt,
     &       label,
     &       mode,
     &       label_list,nop,
     &       op_info)
*----------------------------------------------------------------------*
      case(DEF_EXP_FORMULA)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OP_RES',rule,tgt_info,
     &               val_label=label_list(1))
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list(2:),ndim=nop)
        call get_arg('SWITCH',rule,tgt_info,val_int=ansatz)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_experimental_formula(form_pnt,
     &       title,label_list,nop+1,
     &       ansatz,mode,
     &       op_info,orb_info)
*----------------------------------------------------------------------*
      case(DEF_FORMULA)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('FORMULA',rule,tgt_info,val_str=form_str)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_formula(form_pnt,
     &       form_str,title,
     &       op_info)
*----------------------------------------------------------------------*
      case(EXPAND_OP_PRODUCT)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_arg('NEW',rule,tgt_info,val_log=init)
        if (init) then
          call get_form(form_pnt,trim(label),NEW)
        else
          call get_form(form_pnt,trim(label),OLD)
         end if
        call get_arg('OP_RES',rule,tgt_info,val_label=label)
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list,ndim=nop)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call get_arg('IDX_SV',rule,tgt_info,val_int_list=idx_sv)
        call get_arg('BLK_MIN',rule,tgt_info,val_int_list=iblkmin)
        call get_arg('BLK_MAX',rule,tgt_info,val_int_list=iblkmax)
        call get_arg('CONNECT',rule,tgt_info,val_int_list=connect,
     &               ndim=nconnect)
        nconnect = nconnect/2
        !call get_arg('N_CONNECT',rule,tgt_info,val_int=nconnect)
        call get_arg('AVOID',rule,tgt_info,val_int_list=avoid,
     &               ndim=navoid)
        navoid = navoid/2
        !call get_arg('N_AVOID',rule,tgt_info,val_int=navoid)
        call get_arg('INPROJ',rule,tgt_info,val_int_list=inproj,
     &               ndim=ninproj)
        ninproj = ninproj/4
        !call get_arg('N_INPROJ',rule,tgt_info,val_int=ninproj)
        call get_arg('DESCR',rule,tgt_info,val_label_list=descr,
     &               ndim=ndescr)
        !call get_arg('N_DESCR',rule,tgt_info,val_int=ndescr)
        call get_arg('FAC',rule,tgt_info,val_rl8_list=fac)
        call get_arg('FIX_VTX',rule,tgt_info,val_log=ms_fix)
        call form_expand_op_product(init,form_pnt,fac,
     &       title,label,label_list,nop,
     &       idx_sv,iblkmin,iblkmax,
     &       connect,nconnect,
     &       avoid,navoid,
     &       inproj,ninproj,
     &       descr,ndescr,
     &       ms_fix,op_info,orb_info)
*----------------------------------------------------------------------*
      case(FACTOR_OUT)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('INTERM',rule,tgt_info,
     &       val_label_list=label_list,ndim=nint)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call form_factor_out(form_pnt,form0_pnt,
     &       title,
     &       nint,label_list,
     &       op_info,form_info
     &       )
*----------------------------------------------------------------------*
      case(SUM_HERMIT)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OP_RES',rule,tgt_info,val_label=label)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call form_sum_hermite(form_pnt,form0_pnt,
     &       title,label,
     &       op_info)
*----------------------------------------------------------------------*
      case(EXPAND)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('INTERM',rule,tgt_info,
     &       val_label_list=label_list,ndim=nint)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call get_arg('IMODE',rule,tgt_info,val_int=imode)
        call form_expand_subexpr(form_pnt,form0_pnt,
     &       title,imode,
     &       nint,label_list,
     &       op_info,form_info
     &       )
*----------------------------------------------------------------------*
      case(REPLACE)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OP_LIST',rule,tgt_info,
     &               val_label_list=label_list,ndim=nop)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call form_op_replace_drv(form_pnt,form0_pnt,
     &       title,
     &       nop/2,label_list,
     &       op_info
     &       )
*----------------------------------------------------------------------*
      case(INVARIANT)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call get_arg('OP_RES',rule,tgt_info,
     &       val_label=label)
        call get_arg('OPERATORS',rule,tgt_info,
     &       val_label_list=label_list,ndim=nop)
        call get_arg('REORDER',rule,tgt_info,val_log=reo)
        call form_invariant(form_pnt,form0_pnt,
     &       title,label,
     &       nop,label_list,
     &       reo,op_info
     &       )
*----------------------------------------------------------------------*
      case(REORDER_FORMULA)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call form_reorder(form_pnt,form0_pnt,
     &       title,op_info
     &       )
*----------------------------------------------------------------------*
      case(SUMTERMS)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call get_arg('THRESH',rule,tgt_info,val_rl8=fac(1))
        call form_sum_terms(form_pnt,form0_pnt,title,fac(1),op_info
     &       )
*----------------------------------------------------------------------*
      case(DERIVATIVE)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OP_RES',rule,tgt_info,
     &       val_label=label)
        call get_arg('OP_DERIV',rule,tgt_info,
     &       val_label_list=label_list,ndim=nop)
        call get_arg('OP_MULT',rule,tgt_info,
     &       val_label_list=label_list(nop+1:),ndim=nop2)
        call form_deriv2(form_pnt,form0_pnt,
     &       title,nop,
     &       label,
     &       label_list(1:),
     &       label_list(nop+1:),
     &       op_info)
*----------------------------------------------------------------------*
      case(LEQ_SPLIT)
*----------------------------------------------------------------------*
        call get_arg('LABEL_TRF',rule,tgt_info,val_label=label_list(1))
        call get_form(form_pnt,trim(label_list(1)),ANY) ! pointer not used here
        call get_arg('LABEL_RHS',rule,tgt_info,val_label=label_list(2))
        call get_form(form_pnt,trim(label_list(2)),ANY)
        call get_arg('LABEL_RAW',rule,tgt_info,val_label=label_list(3))
        call get_form(form_pnt,trim(label_list(3)),OLD)
        call get_arg('OP_TRF',rule,tgt_info,val_label=label_list(4))
        call get_arg('OP_RHS',rule,tgt_info,val_label=label_list(5))
        call get_arg('OP_X',rule,tgt_info,
     &       val_label_list=label_list(6:),ndim=nop)
        call get_arg('TITLE_TRF',rule,tgt_info,val_str=title)
        call get_arg('TITLE_RHS',rule,tgt_info,val_str=title2)
        call leq_post(label_list(1),label_list(2),label_list(3),
     &       label_list(4),label_list(5),label_list(6:),nop,
     &       title,title2,
     &       op_info,form_info)
*----------------------------------------------------------------------*
      case(OPTIMIZE)
*----------------------------------------------------------------------*
        call get_arg('LABEL_OPT',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('LABELS_IN',rule,tgt_info,
     &       val_label_list=label_list,ndim=ncat)
        call get_arg('INTERM',rule,tgt_info,
     &       val_label_list=label_list(ncat+1:),ndim=nint)
        call form_opt(form_pnt,
     &       ncat,label_list(1:),
     &       nint,label_list(ncat+1:),
     &       form_info,op_info,str_info,orb_info)

        ! just in case that additional graphs were added:
        call update_strmap(str_info,strmap_info)
*----------------------------------------------------------------------*
      case(PRINT_FORMULA)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),OLD)
        call get_arg('OUTPUT',rule,tgt_info,val_str=title)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call print_formula_drv(form_pnt,title,mode,op_info)
*----------------------------------------------------------------------*
      case(TEX_FORMULA)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),OLD)
        call get_arg('OUTPUT',rule,tgt_info,val_str=title)
        call tex_formula_drv(form_pnt,title,op_info)
*----------------------------------------------------------------------*
      case(SELECT_TERMS)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OP_RES',rule,tgt_info,
     &       val_label=label)
        call get_arg('OP_INCL',rule,tgt_info,
     &       val_label_list=label_list,ndim=ninclude)
        call get_arg('BLK_INCL',rule,tgt_info,
     &       val_int_list=iblk_include,ndim=ntmp)
        if (ninclude.ne.ntmp) call quit(1,'process_rule',
     &       'use same number of elements for OP_INCL, BLK_INCL')
        call get_arg('OP_INCL_OR',rule,tgt_info,
     &       val_label_list=label_list(ninclude+1:),ndim=ninclude_or)
        call get_arg('BLK_INCL_OR',rule,tgt_info,
     &       val_int_list=iblk_include_or,ndim=ntmp)
        if (ninclude_or.ne.ntmp) call quit(1,'process_rule',
     &       'use same number of elements for OP_INCL_OR, BLK_INCL_OR')
        call get_arg('OP_EXCL',rule,tgt_info,
     &       val_label_list=label_list(ninclude+ninclude_or+1:),
     &                                                 ndim=nexclude)
        call get_arg('BLK_EXCL',rule,tgt_info,
     &       val_int_list=iblk_exclude,ndim=ntmp)
        if (nexclude.ne.ntmp) call quit(1,'process_rule',
     &       'use same number of elements for OP_EXCL, BLK_EXCL')
        call form_select_terms(form_pnt,form0_pnt,
     &       label,
     &       ninclude,label_list,iblk_include,
     &       ninclude_or,label_list(ninclude+1:),iblk_include_or,
     &       nexclude,label_list(ninclude+ninclude_or+1:),iblk_exclude,
     &       op_info)
*----------------------------------------------------------------------*
      case(SELECT_SPECIAL)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OPERATORS',rule,tgt_info,
     &       val_label_list=label_list,ndim=nop)
        call get_arg('TYPE',rule,tgt_info,val_str=form_str)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call form_select_special(form_pnt,form0_pnt,
     &       label_list,nop,
     &       form_str,mode,
     &       op_info)
*----------------------------------------------------------------------*
      case(DEL_TERMS)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('TERMS',rule,tgt_info,
     &       val_int_list=idxterms,ndim=nterms)
        call form_del_terms(form_pnt,form0_pnt,
     &       nterms,idxterms,-1,
     &       op_info)
*----------------------------------------------------------------------*
      case(KEEP_TERMS)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('TERMS',rule,tgt_info,
     &       val_int_list=idxterms,ndim=nterms)
        call form_del_terms(form_pnt,form0_pnt,
     &       nterms,idxterms,+1,
     &       op_info)
*----------------------------------------------------------------------*
      case(MODIFY_FACTORIZATION)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('MODIFY',rule,tgt_info,
     &       val_int_list=idxterms,ndim=nterms)
        call form_mod_factorization(form_pnt,form0_pnt,
     &       nterms,idxterms,
     &       op_info)
*----------------------------------------------------------------------*
      case(EXTRACT_ORDER)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OP_RES',rule,tgt_info,val_label=label)
        call get_arg('ORDER',rule,tgt_info,val_int=iorder)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call form_extract_order(form_pnt,form0_pnt,
     &       title, label, iorder, op_info)
*----------------------------------------------------------------------*
      case(EXTRACT_FREQ)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OP_RES',rule,tgt_info,val_label=label)
        call get_arg('ORDER',rule,tgt_info,val_int=iorder)
        allocate(ifreq(iorder),pop_idx(max_pops))
        call get_arg('IDX_FREQ',rule,tgt_info,val_int_list=ifreq)
        call get_arg('NCOMPONENTS',rule,tgt_info,val_int=ncmp)
        call get_arg('IDX_POPS',rule,tgt_info,val_int_list=pop_idx)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call form_extract_freq(form_pnt,form0_pnt,
     &       title, label, iorder, ifreq, ncmp, pop_idx,
     &       op_info)
        deallocate(ifreq,pop_idx)
*----------------------------------------------------------------------*
      case(CLASS_FORMULA)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),OLD)
        call get_arg('OUTPUT',rule,tgt_info,val_str=title)
        call class_formula_drv(form_pnt,title,op_info,1)
*----------------------------------------------------------------------*
      case(SELECT_HERMIT)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OP_RES',rule,tgt_info,val_label=label)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call form_select_hermitian(form_pnt,form0_pnt,
     &       title,label,
     &       op_info)
*----------------------------------------------------------------------*
      case(SELECT_LINE)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OP_RES',rule,tgt_info,val_label=label)
        call get_arg('OP_INCL',rule,tgt_info,
     &       val_label_list=label_list,ndim=ninclude)
        call get_arg('IGAST',rule,tgt_info,val_int=idx)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call form_select_line(form_pnt,form0_pnt,
     &       title,label,
     &       ninclude,label_list,
     &       idx,mode,
     &       op_info
     &       )
*----------------------------------------------------------------------*
      case(DEF_CUMULANTS)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('OP_RES',rule,tgt_info,
     &               val_label=label_list(1))
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list(2:),ndim=nop)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('LEVEL',rule,tgt_info,val_int=level)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_cumulants(form_pnt,
     &       title,label_list,nop+1,mode,level,op_info)
*----------------------------------------------------------------------*
      case(INSERT)
*----------------------------------------------------------------------*
        call get_arg('LABEL_RES',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),ANY)
        call get_arg('LABEL_IN',rule,tgt_info,val_label=label)
        call get_form(form0_pnt,trim(label),OLD)
        call get_arg('OP_RES',rule,tgt_info,val_label=label)
        call get_arg('OP_INS',rule,tgt_info,val_label=label2)
        call get_arg('OP_INCL',rule,tgt_info,
     &       val_label_list=label_list,ndim=ninclude)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call form_insert_op(form_pnt,form0_pnt,
     &       title,label,label2,
     &       ninclude,label_list,
     &       op_info)
*----------------------------------------------------------------------*
      case(DEF_MRCC_INTM)
*----------------------------------------------------------------------*
        call get_arg('LABEL',rule,tgt_info,val_label=label)
        call get_form(form_pnt,trim(label),NEW)
        call get_arg('INTERM',rule,tgt_info,
     &               val_label=label)
        call get_arg('OPERATORS',rule,tgt_info,
     &               val_label_list=label_list,ndim=nop)
        call get_arg('MAXCOM',rule,tgt_info,val_int=ansatz)
        call get_arg('FAC',rule,tgt_info,val_rl8_list=fac)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('TITLE',rule,tgt_info,val_str=title)
        call set_mrcc_intermediates(form_pnt,
     &         title,label,label_list,
     &         nop,ansatz,fac,mode,op_info)

*----------------------------------------------------------------------*
*     subsection ME-LISTS
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      case(DEF_ME_LIST)
*----------------------------------------------------------------------*
        call get_arg('LIST',rule,tgt_info,val_label=label)
        call get_arg('OPERATOR',rule,tgt_info,val_label=label2)
        call get_arg('AB_SYM',rule,tgt_info,val_int=absym)
        call get_arg('CA_SYM',rule,tgt_info,val_int=casym)
        call get_arg('IRREP',rule,tgt_info,val_int=gamma)
        call get_arg('S2',rule,tgt_info,val_int=s2)
        call get_arg('2MS',rule,tgt_info,val_int=ms)
        call get_arg('MS_FIX',rule,tgt_info,val_log=ms_fix)
        call get_arg('DIAG_TYPE',rule,tgt_info,val_int=imode)
        call get_arg('DIAG_IRREP',rule,tgt_info,val_int=dgam)
        call get_arg('DIAG_MS',rule,tgt_info,val_int=dms)
        call get_arg('MIN_REC',rule,tgt_info,val_int=min_rank)
        call get_arg('MAX_REC',rule,tgt_info,val_int=max_rank)
        call get_arg('REC',rule,tgt_info,val_int=idx)
        call define_me_list(label,label2,
     &       absym,casym,gamma,s2,ms,ms_fix,
     &       idx,min_rank,max_rank,imode,dgam,dms,
     &       op_info,orb_info,str_info,strmap_info)
*----------------------------------------------------------------------*
      case(UNITY)
*----------------------------------------------------------------------*
        call get_arg('LIST',rule,tgt_info,val_label=label)
        call get_arg('FAC',rule,tgt_info,val_rl8=fac(1))
        call get_arg('INIT',rule,tgt_info,val_log=init)
        call get_arg('MIN_BLK',rule,tgt_info,val_int=minblk)
        call get_arg('MAX_BLK',rule,tgt_info,val_int=maxblk)
        call get_arg('MS_SYM_SIGN',rule,tgt_info,val_int=imode)

        call add_unity_drv(label,fac(1),imode,init,minblk,maxblk,
     &       op_info,orb_info,str_info)

*----------------------------------------------------------------------*
      case(ASSIGN_ME2OP)
*----------------------------------------------------------------------*

        call get_arg('LIST',rule,tgt_info,val_label=label)
        call get_arg('OPERATOR',rule,tgt_info,val_label=label2)        

        call assign_me_list(label,label2,op_info)

*----------------------------------------------------------------------*
      case(RES_ME_LIST)
*----------------------------------------------------------------------*
        
        call get_arg('LIST',rule,tgt_info,val_label=label)

        call reset_me_list(label,op_info)

*----------------------------------------------------------------------*
      case(DELETE_ME_LIST)
*----------------------------------------------------------------------*
        
        call get_arg('LIST',rule,tgt_info,val_label=label)

        call del_me_list(label,op_info)

*----------------------------------------------------------------------*
      case(IMPORT)
*----------------------------------------------------------------------*

        call get_arg('LIST',rule,tgt_info,val_label=label)
        call get_arg('TYPE',rule,tgt_info,val_str=list_type)
        call get_arg('ENV',rule,tgt_info,val_str=env_type)

        if (form_test) return

        call import_op_el(label,
     &       list_type,env_type,
     &       op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
      case(GETEST)
*----------------------------------------------------------------------*

        call get_arg('LIST',rule,tgt_info,val_label=label)
        call get_arg('R-SYS',rule,tgt_info,val_int_list=iRdef,
     &               ndim=norb)
c dbg
        print *,'norb = ',norb
        print *,'iRdef = ',iRdef(1:norb)
c dbg
        ! trap, if we get so far ...
        if (norb.gt.maxterms) 
     &       call quit(1,'process_rule','norb.gt.maxterms')
        call get_arg('CASE',rule,tgt_info,val_int=icase)
        call get_arg('SPLIT-FOCK',rule,tgt_info,val_int=icaseF)

        if (form_test) return

        call mod_op_for_ge_test(label,
     &       iRdef,norb,icase,icaseF,
     &       op_info,str_info,strmap_info,orb_info)

*----------------------------------------------------------------------*
      case(SET_FREQ)
*----------------------------------------------------------------------*

        call get_arg('LIST',rule,tgt_info,val_label=label)
        call get_arg('FREQ',rule,tgt_info,val_rl8=freq)
        call get_mel(mel_pnt,label,OLD)
        call set_frequency(mel_pnt,freq)

*----------------------------------------------------------------------*
      case(PRINT_RES)
*----------------------------------------------------------------------*

        call get_arg('LIST',rule,tgt_info,val_label=label)
        call get_arg('ORDER',rule,tgt_info,val_int=iorder)
        allocate(ifreq(iorder))
        call get_arg('IDX_FREQ',rule,tgt_info,val_int_list=ifreq)
        call get_mel(mel_pnt,label,OLD)
        
        if (form_test) return

        call print_result(iorder,ifreq,mel_pnt,.false.,orb_info)
        deallocate(ifreq)

*----------------------------------------------------------------------*
      case(PRINT_MEL)
*----------------------------------------------------------------------*

        call get_arg('LIST',rule,tgt_info,val_label=label)
        call get_arg('COMMENT',rule,tgt_info,val_str=title)
        call get_arg('FORMAT',rule,tgt_info,val_str=mode)
        call get_arg('CHECK_THRESH',rule,tgt_info,val_rl8=fac(1))
        call get_arg('EXPECTED',rule,tgt_info,val_rl8=fac(2))

        if (form_test) return

        call get_mel(mel_pnt,label,OLD)
        call print_list(title,mel_pnt,mode,fac(1),fac(2),
     &                  orb_info,str_info)

*----------------------------------------------------------------------*
      case(PRINT_MEL_INFO_)
*----------------------------------------------------------------------*

        call get_arg('LIST',rule,tgt_info,val_label=label)

        if (form_test) return

        call get_mel(mel_pnt,label,OLD)
        call print_mel_info(luout,mel_pnt)

*----------------------------------------------------------------------*
      case(PRINT_)
*----------------------------------------------------------------------*

        call get_arg('STRING',rule,tgt_info,val_str=strscr)

        write(luout,*) trim(strscr)

*----------------------------------------------------------------------*
      case(SET_MEL)
*----------------------------------------------------------------------*

        call get_arg('LIST',rule,tgt_info,val_label=label)
        call get_arg('VAL_LIST',rule,tgt_info,val_rl8_list=fac)
        call get_arg('IDX_LIST',rule,tgt_info,
     &       val_int_list=idxblk,ndim=nblk)

        if (form_test) return

        call get_mel(mel_pnt,label,OLD)

        call zeroop(mel_pnt)
        call set_list(mel_pnt,idxblk,fac,nblk)

*----------------------------------------------------------------------*
      case(EXTRACT_DIAG)
*----------------------------------------------------------------------*
        call get_arg('LIST_RES',rule,tgt_info,val_label=label)
        call get_arg('LIST_IN',rule,tgt_info,val_label=label2)
        call get_arg('MODE',rule,tgt_info,val_str=mode)

        if (form_test) return

        call dia_from_op(label,label2,mode,
     &       op_info,str_info,orb_info)

*----------------------------------------------------------------------*
      case(REORDER_MEL)
*----------------------------------------------------------------------*
        call get_arg('LIST_RES',rule,tgt_info,val_label=label)
        call get_arg('LIST_IN',rule,tgt_info,val_label=label2)
        call get_arg('FROMTO',rule,tgt_info,val_int=idx)
        call get_arg('ADJOINT',rule,tgt_info,val_log=dagger)
        call get_arg('SEARCH',rule,tgt_info,val_log=init)

        if (form_test) return

        call reo_mel(label,label2,init,
     &       op_info,str_info,strmap_info,orb_info,idx,dagger)

*----------------------------------------------------------------------*
      case(SPIN_PROJECT)
*----------------------------------------------------------------------*
        call get_arg('LIST',rule,tgt_info,val_label=label)
        call get_arg('S2',rule,tgt_info,val_int=s2)

        if (form_test) return

        call spin_prj_list_drv(label,s2,
     &       op_info,str_info,strmap_info,orb_info)

*----------------------------------------------------------------------*
      case(ORB_FLIP)
*----------------------------------------------------------------------*

        call get_arg('LIST',rule,tgt_info,val_label=label)

        if (form_test) return

        call get_mel(mel_pnt,label,OLD)
        call orb_flip_mel(mel_pnt,str_info,orb_info)

*----------------------------------------------------------------------*
*     subsection EVALUATE
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      case(PRECONDITIONER)
*----------------------------------------------------------------------*

        call get_arg('LIST_PRC',rule,tgt_info,val_label=label)
        call get_arg('LIST_INP',rule,tgt_info,val_label_list=label_list,
     &       ndim=nop)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('SHIFT',rule,tgt_info,val_rl8=fac(1))
        call get_arg('THRES',rule,tgt_info,val_rl8=fac(2))

        if (form_test) return

C set these directly now:
c          mode = 'dia-F'
c          if (rule%n_parameter_strings.eq.2) mode(5:5) = 'H'
c          if (rule%n_parameter_strings.eq.3) mode(5:8) = 'F+id'
c          mode = 'dia-R12'

        call set_prc4op(label,mode,fac(1),
     &       label_list,nop,fac(2),
     &       op_info,str_info,orb_info)

*----------------------------------------------------------------------*
      case(EVAL)
*----------------------------------------------------------------------*

        call get_arg('FORM',rule,tgt_info,val_label=label)
        call get_arg('INIT',rule,tgt_info,val_log=init)

        if (form_test) return
        call evaluate(label,init,
     &       op_info,form_info,str_info,strmap_info,orb_info)

*----------------------------------------------------------------------*
      case(INVERT)
*----------------------------------------------------------------------*
        call get_arg('LIST_INV',rule,tgt_info,
     &               val_label_list=label_list,ndim=nop)
        call get_arg('LIST',rule,tgt_info,
     &               val_label_list=label_list(nop+1:),ndim=nop2)
        call get_arg('MODE',rule,tgt_info,val_str=mode)

        if (form_test) return
        call inv_op(nop,label_list(1:),nop2,label_list(nop+1:),mode,
     &       op_info,orb_info,str_info,strmap_info)

*----------------------------------------------------------------------*
      case(ADD)
*----------------------------------------------------------------------*

        call get_arg('LIST_SUM',rule,tgt_info,val_label=label)
        call get_arg('LISTS',rule,tgt_info,
     &       val_label_list=label_list,ndim=nfac)
        call get_arg('FAC',rule,tgt_info,val_rl8_list=fac,ndim=nfac)
        call get_arg('REPLACE',rule,tgt_info,val_log=init)
        
        call add_op(label,fac,label_list,nfac,
     &       op_info,orb_info,str_info,init)

*----------------------------------------------------------------------*
      case(SCALE)
*----------------------------------------------------------------------*

        call get_arg('LIST_RES',rule,tgt_info,val_label=label)
        call get_arg('LIST_INP',rule,tgt_info,val_label=label_list(1))
        call get_arg('LIST_SCAL',rule,tgt_info,val_label=label_list(2))
        call get_arg('FAC',rule,tgt_info,val_rl8_list=fac,ndim=nfac)
        !call get_arg('NFAC',rule,tgt_info,val_int=nfac)
        call get_arg('IDX_LIST',rule,tgt_info,val_int_list=idxblk)

        imode = 2
        if (len_trim(label_list(2)).eq.0.or.label_list(2)(1:1).eq.'-')
     &       imode = 1
        call scale_op(label,
     &       imode,idxblk,fac,label_list,nfac,
     &       op_info,orb_info,str_info)

*----------------------------------------------------------------------*
      case(SCALE_COPY)
*----------------------------------------------------------------------*

        call get_arg('LIST_RES',rule,tgt_info,val_label=label)
        call get_arg('LIST_INP',rule,tgt_info,val_label=label_list(1))
        call get_arg('LIST_SHAPE',rule,tgt_info,
     &               val_label_list=label_list(2:),ndim=nspcfrm)
        call get_arg('FAC',rule,tgt_info,val_rl8_list=fac,ndim=nfac)
        call get_arg('MODE',rule,tgt_info,val_str=mode)

        if (form_test) return

        call scale_copy_op(label,label_list,fac,nfac,mode,nspcfrm,
     &       op_info,orb_info,str_info)

*----------------------------------------------------------------------*
      case(EVALPROP)
*----------------------------------------------------------------------*

        call get_arg('DENS',rule,tgt_info,
     &       val_label_list=label_list,ndim=ndens)
        call get_arg('ENV',rule,tgt_info,val_str=env_type)
        call get_arg('RANK',rule,tgt_info,val_int=rank)

        call prop_evaluate(ndens,rank,label_list,
     &       env_type,op_info,str_info,orb_info)

*----------------------------------------------------------------------*
      case(SOLVENLEQ)
*----------------------------------------------------------------------*

        call get_arg('LIST_OPT',rule,tgt_info,
     &       val_label_list=label_list(1:),ndim=nopt)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('LIST_RESID',rule,tgt_info,
     &       val_label_list=label_list(nopt+1:))
        call get_arg('LIST_PRC',rule,tgt_info,
     &       val_label_list=label_list(2*nopt+1:))
        call get_arg('LIST_SPC',rule,tgt_info,
     &       val_label_list=label_list(3*nopt+1:),ndim=nspecial)
        call get_arg('LIST_E',rule,tgt_info,
     &       val_label=label)
        call get_arg('FORM',rule,tgt_info,
     &       val_label=label2)
        call get_arg('FORM_SPC',rule,tgt_info,
     &       val_label_list=label_list(3*nopt+nspecial+1:),ndim=nspcfrm)

        if (form_test) return

        call solve_nleq(mode,nopt,
     &       label_list(1:nopt),               ! to be opt.
     &       label_list(nopt+1:nopt+nopt),     ! residual
     &       label_list(2*nopt+1:2*nopt+nopt), ! precond.
     &       label,                            ! energy
     &       label2,                           ! formula
     &       label_list(3*nopt+1:
     &                  3*nopt+nspecial),
     &          nspecial,
     &       label_list(3*nopt+nspecial+1:
     &                  3*nopt+nspecial+nspcfrm),
     &          nspcfrm,                       ! specials
     &       op_info,form_info,str_info,strmap_info,orb_info)

*----------------------------------------------------------------------*
      case(SOLVELEQ)
*----------------------------------------------------------------------*

        call get_arg('LIST_OPT',rule,tgt_info,
     &       val_label_list=label_list(1:),ndim=nopt)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('N_ROOTS',rule,tgt_info,val_int=nroots)
        call get_arg('LIST_PRC',rule,tgt_info,
     &       val_label_list=label_list(1*nopt+1:))
        call get_arg('OP_MVP',rule,tgt_info,
     &       val_label_list=label_list(2*nopt+1:))
        call get_arg('OP_SVP',rule,tgt_info,
     &       val_label_list=label_list(3*nopt+1:))
        call get_arg('OP_RHS',rule,tgt_info,
     &       val_label_list=label_list(4*nopt+1:))
        call get_arg('LIST_SPC',rule,tgt_info,
     &       val_label_list=label_list(5*nopt+1:),ndim=nspecial)
        call get_arg('FORM_SPC',rule,tgt_info,
     &       val_label_list=label_list(5*nopt+nspecial+1:),ndim=nspcfrm)
        call get_arg('FORM',rule,tgt_info,
     &       val_label=label)

        if (form_test) return

        call solve_leq(mode,nopt,nroots,
     &       label_list(1:nopt),               ! to be opt.
     &       label_list(  nopt+1:  nopt+nopt), ! precond.
     &       label_list(2*nopt+1:2*nopt+nopt), ! mvp-labels
     &       label_list(3*nopt+1:3*nopt+nopt), ! metric-labels
     &       label_list(4*nopt+1:4*nopt+nopt), ! rhs-labels
     &       xdum,                             ! dummy
     &       label,                            ! formula
     &       label_list(5*nopt+1:
     &                   5*nopt+nspecial),nspecial,
     &       label_list(5*nopt+nspecial+1:
     &                  5*nopt+nspecial+nspcfrm),nspcfrm,0d0,
     &       op_info,form_info,str_info,strmap_info,orb_info)

*----------------------------------------------------------------------*
      case(SOLVEEVP)
*----------------------------------------------------------------------*

        call get_arg('LIST_OPT',rule,tgt_info,
     &       val_label_list=label_list(1:),ndim=nopt)
        call get_arg('MODE',rule,tgt_info,val_str=mode)
        call get_arg('N_ROOTS',rule,tgt_info,val_int=nroots)
        call get_arg('TARG_ROOT',rule,tgt_info,val_int=targ_root)
        if (targ_root.le.0) targ_root=nroots
        call get_arg('LIST_PRC',rule,tgt_info,
     &       val_label_list=label_list(1*nopt+1:))
        call get_arg('OP_MVP',rule,tgt_info,
     &       val_label_list=label_list(2*nopt+1:))
        call get_arg('OP_SVP',rule,tgt_info,
     &       val_label_list=label_list(3*nopt+1:))
        call get_arg('LIST_SPC',rule,tgt_info,
     &       val_label_list=label_list(4*nopt+1:),ndim=nspecial)
        call get_arg('FORM_SPC',rule,tgt_info,
     &       val_label_list=label_list(4*nopt+nspecial+1:),ndim=nspcfrm)
        call get_arg('FORM',rule,tgt_info,
     &       val_label=label)

        if (form_test) return

        call solve_evp(mode,nopt,nroots,targ_root,
     &       label_list(1:nopt),               ! to be opt.
     &       label_list(  nopt+1:  nopt+nopt), ! precond.
     &       label_list(2*nopt+1:2*nopt+nopt), ! mvp-labels
     &       label_list(3*nopt+1:3*nopt+nopt), ! metric-labels
     &       label,                            ! formula
     &       label_list(4*nopt+1:
     &                   4*nopt+nspecial),nspecial,
     &       label_list(4*nopt+nspecial+1:     ! spec. form.
     &                  4*nopt+nspecial+nspcfrm),
     &          nspcfrm,0d0,
     &       op_info,form_info,str_info,strmap_info,orb_info)


*----------------------------------------------------------------------*
*     unknown:
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      case default
*----------------------------------------------------------------------*
        call quit(1,'process_rule','unknown command: '//
     &       trim(rule%command))
      end select

      ! do not put anything here as we might have returned
      ! from this subroutine previously (see above)

*----------------------------------------------------------------------*
*     end of subroutine
*----------------------------------------------------------------------*

      return

*----------------------------------------------------------------------*
*     embedded subroutines follow
*----------------------------------------------------------------------*
      contains

      subroutine get_op(op_pnt,label,status)
      ! return pointer to operator
      ! status NEW: label must not exist yet, generate new entry
      ! status OLD: label must exist already
      ! status ANY: return existing operator, or generate new entry

      type(operator), intent(out), pointer ::
     &     op_pnt
      character(len=*), intent(in) ::
     &     label
      integer, intent(in) ::
     &     status

      integer ::
     &     idx

      idx = idx_oplist2(trim(label),op_info)

      select case(status)
      case(NEW)
        if (idx.gt.0)
     &       call quit(0,'process_rule',
     &       'operator does already exist: '//trim(label))
        call add_operator(trim(label),op_info)
        idx = idx_oplist2(trim(label),op_info)
      case(OLD)
        if (idx.le.0)
     &       call quit(0,'process_rule',
     &       'operator does not exist: '//trim(label))
      case(ANY)
        if (idx.le.0) then
          call add_operator(trim(label),op_info)
          idx = idx_oplist2(trim(label),op_info)
        end if
      end select

      op_pnt => op_info%op_arr(idx)%op

      return

      end subroutine

      subroutine get_mel(mel_pnt,label,status)
      ! return pointer to operator
      ! status NEW: label must not exist yet, generate new entry
      ! status OLD: label must exist already
      ! status ANY: return existing operator, or generate new entry

      type(me_list), intent(out), pointer ::
     &     mel_pnt
      character(len=*), intent(in) ::
     &     label
      integer, intent(in) ::
     &     status

      integer ::
     &     idx

      idx = idx_mel_list(trim(label),op_info)

      select case(status)
      case(NEW)
        call quit(0,'process_rule',
     &       'get_mel(NEW) '//trim(label)//': use DEF_ME_LIST')
      case(OLD)
        if (idx.le.0)
     &       call quit(0,'process_rule',
     &       'ME list does not exist: '//trim(label))
      case(ANY)
        call quit(0,'process_rule',
     &       'get_mel(ANY) '//trim(label)//': use DEF_ME_LIST')
      end select

      mel_pnt => op_info%mel_arr(idx)%mel

      return

      end subroutine

      subroutine get_form(form_pnt,label,status)
      ! return pointer to formula
      ! status NEW: label must not exist yet, generate new entry
      ! status OLD: label must exist already
      ! status ANY: return existing formula, or generate new entry

      type(formula), intent(out), pointer ::
     &     form_pnt
      character(len=*), intent(in) ::
     &     label
      integer, intent(in) ::
     &     status

      integer ::
     &     idx

      idx = idx_formlist(trim(label),form_info)

      select case(status)
      case(NEW)
        if (idx.gt.0)
     &       call quit(0,'process_rule',
     &       'formula does already exist: '//trim(label))
        call add_formula(form_info,trim(label))
        idx = idx_formlist(trim(label),form_info)
      case(OLD)
        if (idx.le.0)
     &       call quit(0,'process_rule',
     &       'formula does not exist: '//trim(label))
      case(ANY)
        if (idx.le.0) then
          call add_formula(form_info,trim(label))
          idx = idx_formlist(trim(label),form_info)
        end if
      end select

      form_pnt => form_info%form_arr(idx)%form

      return

      end subroutine

      end
