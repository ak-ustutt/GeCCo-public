*----------------------------------------------------------------------*
      subroutine set_r12_general_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
*     set targets needed in more or less all kinds of R12 calculations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     hpvx_constr(2,ngastp,2),
     &     gas_constr(2,orb_info%ngas,2,2)

      integer ::
     &     min_rank, max_rank, ansatz,
     &     isim, ncat, nint, icnt, nlab,
     &     isym, ms, msc, sym_arr(8)
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters(3)
      character(12) ::
     &     approx

      character(*), intent(in) ::
     &     env_type

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting general targets for R12 ...'

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! the formal R12 geminal: P12 r12|0>
      call add_target(op_r12,ttype_op,.false.,tgt_info)
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      min_rank = 2  ! 1 is a possibility 
      call r12gem_parameters(-1,parameters,
     &                   .false.,min_rank,ansatz)
      call set_rule(op_r12,ttype_op,DEF_R12GEMINAL,
     &              op_r12,1,1,
     &              parameters,1,tgt_info)

      ! the adjoint (should be obsolete soon)
      call add_target(op_rba,ttype_op,.false.,tgt_info)
      call set_dependency(op_rba,op_r12,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_r12,.true.) ! <- dagger=.true.
      call set_rule(op_rba,ttype_op,CLONE_OP,
     &              op_rba,1,1,
     &              parameters,1,tgt_info)

      ! the coefficients
      call add_target(op_c12,ttype_op,.false.,tgt_info)
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call xop_parameters(-1,parameters,
     &     .false.,min_rank,max_rank,0,max_rank+1)
      call set_rule(op_c12,ttype_op,DEF_R12COEFF,
     &              op_c12,1,1,
     &              parameters,1,tgt_info)

      ! Lagrange multipliers associated with coefficients
      call add_target(op_cba,ttype_op,.false.,tgt_info)
      call set_dependency(op_cba,op_c12,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_c12,.true.) ! <- dagger=.true.
      call set_rule(op_cba,ttype_op,CLONE_OP,
     &              op_cba,1,1,
     &              parameters,1,tgt_info)

      ! Preconditioner
      call add_target(op_diar12,ttype_op,.false.,tgt_info)
      call set_dependency(op_diar12,op_c12,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_c12,.false.) ! <- dagger=.false.
      call set_rule(op_diar12,ttype_op,CLONE_OP,
     &              op_diar12,1,1,
     &              parameters,1,tgt_info)

      ! Now: the operators associated with the actual R12 integrals:
      !  <pq'|r12|ij> 
      call add_target(op_rint,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,min_rank,2,0,2)
      call set_rule(op_rint,ttype_op,DEF_R12INT,
     &              op_rint,1,1,
     &              parameters,1,tgt_info)
      
      ! soon obsolete: adjoint of the above ints ...
      call add_target(op_rinba,ttype_op,.false.,tgt_info)
      call set_dependency(op_rinba,op_rint,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rint,.true.) ! <- dagger=.true.
      call set_rule(op_rinba,ttype_op,CLONE_OP,
     &              op_rinba,1,1,
     &              parameters,1,tgt_info)
      
      ! commutator integrals <ab|[T1+T2,r12]|cd>
      call add_target(op_ttr,ttype_op,.false.,tgt_info)
      call set_dependency(op_ttr,op_rint,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rint,.false.) ! <- dagger=.false.
      call set_rule(op_ttr,ttype_op,CLONE_OP,
     &              op_ttr,1,1,
     &              parameters,1,tgt_info)
      
c      ! the adjoint
c      call add_target(op_ttr_bar,ttype_op,.false.,tgt_info)
c      call set_dependency(op_ttr_bar,op_ttr,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_ttr,.true.) ! <- dagger=.true.
c      call set_rule(op_ttr_bar,ttype_op,CLONE_OP,
c     &              op_ttr_bar,1,1,
c     &              parameters,1,tgt_info)
      
      ! V^{ij}_{kl}
      call add_target(op_v_inter,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,2,2,0,2)
      call set_rule(op_v_inter,ttype_op,DEF_R12INTERM,
     &              op_v_inter,1,1,
     &              parameters,1,tgt_info)
      
      ! the adjoint
      call add_target(op_vbar_inter,ttype_op,.false.,tgt_info)
      call set_dependency(op_vbar_inter,op_v_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_v_inter,.false.) ! <- dagger=.false.
                                   ! we really need the transposed list
      call set_rule(op_vbar_inter,ttype_op,CLONE_OP,
     &              op_vbar_inter,1,1,
     &              parameters,1,tgt_info)
      
      ! R12^{2} integrals
      call add_target(op_f2,ttype_op,.false.,tgt_info)
      call set_dependency(op_f2,op_v_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_v_inter,.false.) ! <- dagger=.false.
      call set_rule(op_f2,ttype_op,CLONE_OP,
     &              op_f2,1,1,
     &              parameters,1,tgt_info)

      ! B intermediate
      call add_target(op_b_inter,ttype_op,.false.,tgt_info)
      call set_dependency(op_b_inter,op_v_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_v_inter,.false.) ! <- dagger=.false.
      call set_rule(op_b_inter,ttype_op,CLONE_OP,
     &              op_b_inter,1,1,
     &              parameters,1,tgt_info)

      ! X intermediate
      call add_target(op_x_inter,ttype_op,.false.,tgt_info)
      call set_dependency(op_x_inter,op_v_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_v_inter,.false.) ! <- dagger=.false.
      call set_rule(op_x_inter,ttype_op,CLONE_OP,
     &              op_x_inter,1,1,
     &              parameters,1,tgt_info)

      ! inverse of B
      call add_target(op_b_inv,ttype_op,.false.,tgt_info)
      call set_dependency(op_b_inv,op_v_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_v_inter,.false.) ! <- dagger=.false.
      call set_rule(op_b_inv,ttype_op,CLONE_OP,
     &              op_b_inv,1,1,
     &              parameters,1,tgt_info)

      ! inverse of X
      call add_target(op_x_inv,ttype_op,.false.,tgt_info)
      call set_dependency(op_x_inv,op_v_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_v_inter,.false.) ! <- dagger=.false.
      call set_rule(op_x_inv,ttype_op,CLONE_OP,
     &              op_x_inv,1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*
      ! set approx string
      approx(1:12) = ' '
      call get_argument_value('method.R12','approx',str=approx)

      ! formal definition of V
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_vint
      labels(2) = op_v_inter
      labels(3) = op_ham
      labels(4) = op_r12
      call add_target(form_r12_vint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_vint,op_v_inter,tgt_info)
      call set_dependency(form_r12_vint,op_ham,tgt_info)
      call set_dependency(form_r12_vint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vint,0,'gxr')
      call set_rule(form_r12_vint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to V
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_vcabs
      labels(2) = op_v_inter
      labels(3) = op_ham
      labels(4) = op_rint
      labels(5) = op_unity
      call add_target(form_r12_vcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_vcabs,op_v_inter,tgt_info)
      call set_dependency(form_r12_vcabs,op_unity,tgt_info)
      call set_dependency(form_r12_vcabs,op_ham,tgt_info)
      call set_dependency(form_r12_vcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vcabs,ansatz,'V '//approx)
      call set_rule(form_r12_vcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! to be removed soon: transpose of V
      ! formal definition of Vbar (viz. V^dagger)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_vbint
      labels(2) = op_vbar_inter
      labels(3) = op_rba
      labels(4) = op_ham
      call add_target(form_r12_vbint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_vbint,op_vbar_inter,tgt_info)
      call set_dependency(form_r12_vbint,op_ham,tgt_info)
      call set_dependency(form_r12_vbint,op_rba,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vbint,0,'rxg')
      call set_rule(form_r12_vbint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to Vbar
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_vbcabs
      labels(2) = op_vbar_inter
      labels(3) = op_rinba
      labels(4) = op_ham
      labels(5) = op_unity
      call add_target(form_r12_vbcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_vbcabs,op_vbar_inter,tgt_info)
      call set_dependency(form_r12_vbcabs,op_unity,tgt_info)
      call set_dependency(form_r12_vbcabs,op_ham,tgt_info)
      call set_dependency(form_r12_vbcabs,op_rinba,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vbcabs,ansatz,'V+'//approx)
      call set_rule(form_r12_vbcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! formal definition of X
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xint
      labels(2) = op_x_inter
      labels(3) = op_rba
      labels(4) = op_r12
      call add_target(form_r12_xint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xint,op_x_inter,tgt_info)
      call set_dependency(form_r12_xint,op_rba,tgt_info)
      call set_dependency(form_r12_xint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xint,0,'rxr')
      call set_rule(form_r12_xint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to X
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xcabs
      labels(2) = op_x_inter
      labels(3) = op_rinba
      labels(4) = op_rint
      labels(5) = op_f2
      call add_target(form_r12_xcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xcabs,op_x_inter,tgt_info)
      call set_dependency(form_r12_xcabs,op_f2,tgt_info)
      call set_dependency(form_r12_xcabs,op_rinba,tgt_info)
      call set_dependency(form_r12_xcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xcabs,ansatz,'X '//approx)
      call set_rule(form_r12_xcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! formal definition of B
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_bint
      labels(2) = op_b_inter
      labels(3) = op_rba
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_r12_bint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_bint,op_b_inter,tgt_info)
      call set_dependency(form_r12_bint,op_rba,tgt_info)
      call set_dependency(form_r12_bint,op_ham,tgt_info)
      call set_dependency(form_r12_bint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_bint,0,'rfxr')
      call set_rule(form_r12_bint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)


      ! CABS approximation to B
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_bcabs
      labels(2) = op_b_inter
      labels(3) = op_rinba
      labels(4) = op_ttr
      labels(5) = op_unity
      nlab = 5
      call add_target(form_r12_bcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_bcabs,op_b_inter,tgt_info)
      call set_dependency(form_r12_bcabs,op_unity,tgt_info)
      call set_dependency(form_r12_bcabs,op_rinba,tgt_info)
      call set_dependency(form_r12_bcabs,op_rint,tgt_info)
      if (approx(1:2).ne.'A ') then
        call set_dependency(form_r12_bcabs,op_x_inter,tgt_info)
        call set_dependency(form_r12_bcabs,form_r12_xcabs,tgt_info)
        call set_dependency(form_r12_bcabs,op_ham,tgt_info)
        labels(6) = op_x_inter
        labels(7) = op_ham
        nlab = 7
      end if
      approx(12:12) = 'S' ! set symmetrization flag
      call form_parameters(-1,
     &     parameters,2,title_r12_bcabs,ansatz,'B '//approx)
      approx(12:12) = ' ' ! unset flag
      call set_rule(form_r12_bcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,nlab,1,
     &              parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*
      ! set V
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_vcabs
      labels(2) = form_r12_vcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_vcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_vcabs,form_r12_vcabs,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_v_def,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_ham,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_vcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! set V+
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_vbcabs
      labels(2) = form_r12_vbcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_vbcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_vbcabs,form_r12_vbcabs,tgt_info)
      call set_dependency(fopt_r12_vbcabs,mel_vbar_def,tgt_info)
      call set_dependency(fopt_r12_vbcabs,mel_ham,tgt_info)
      call set_dependency(fopt_r12_vbcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_vbcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! set X
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_xcabs
      labels(2) = form_r12_xcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_xcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_xcabs,form_r12_xcabs,tgt_info)
      call set_dependency(fopt_r12_xcabs,mel_x_def,tgt_info)
      call set_dependency(fopt_r12_xcabs,mel_f2,tgt_info)
      call set_dependency(fopt_r12_xcabs,mel_rinba,tgt_info)
      call set_dependency(fopt_r12_xcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_xcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! set B
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_bcabs
      labels(2) = form_r12_bcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_bcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_bcabs,form_r12_bcabs,tgt_info)
      call set_dependency(fopt_r12_bcabs,mel_b_def,tgt_info)
      call set_dependency(fopt_r12_bcabs,mel_ham,tgt_info)
      call set_dependency(fopt_r12_bcabs,mel_rinba,tgt_info)      
      call set_dependency(fopt_r12_bcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_bcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*
      ! ----------------------------
      ! A) integrals to be imported:
      ! ----------------------------
      ! R12integrals
      call add_target(mel_rint,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_rint,op_rint,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rint
      labels(2) = op_rint
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_rint,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rint
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_rint,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! to be changed soon:
      ! adjoint of R12integrals
      call add_target(mel_rinba,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_rinba,op_rinba,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rinba
      labels(2) = op_rinba
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_rinba,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rinba
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_rinba,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! [T1+T2,R12] integrals
      call add_target(mel_ttr,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ttr,op_ttr,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ttr
      labels(2) = op_ttr
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_ttr,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ttr
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_ttr,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! R12^2 integrals
      call add_target(mel_f2,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_f2,op_f2,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_f2
      labels(2) = op_f2
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_f2,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_f2
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_f2,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! ---------------------------------------
      ! B) definition of list for intermediates
      ! ---------------------------------------
      ! V-list
      call add_target(mel_v_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_v_def,op_v_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_v_inter
      labels(2) = op_v_inter
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_v_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! V+-list
      call add_target(mel_vbar_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_vbar_def,op_vbar_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_vbar_inter
      labels(2) = op_vbar_inter
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_vbar_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)


      ! X-list
      call add_target(mel_x_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_x_def,op_x_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_x_inter
      labels(2) = op_x_inter
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_x_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! B-list
      call add_target(mel_b_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_b_def,op_b_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_b_inter
      labels(2) = op_b_inter
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_b_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! B^-1 for "diagonal"
      call add_target(mel_b_inv,ttype_opme,.true.,tgt_info)
      call set_dependency(mel_b_inv,op_diar12,tgt_info)
      call set_dependency(mel_b_inv,eval_r12_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_b_inv
      labels(2) = op_diar12
c      labels(2) = op_b_inter ! actually, B^-1 should have the contravariant shape
c                             ! but as long as we do not formally calculate with
c                             ! this entity this does not matter
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_b_inv,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1) = mel_b_inv   ! output
      labels(2) = mel_b_inter ! input
      call set_rule(mel_b_inv,ttype_opme,INVERT,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! diagonal of B(ij) for testing
      call add_target(mel_b_dia,ttype_opme,.true.,tgt_info)
      call set_dependency(mel_b_dia,op_diar12,tgt_info)
      call set_dependency(mel_b_dia,eval_r12_inter,tgt_info)
      call set_dependency(mel_b_dia,mel_ham,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_b_dia
      labels(2) = op_diar12
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_b_dia,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1) = mel_b_dia   ! output
      labels(2) = mel_ham     ! input
      labels(3) = mel_b_inter ! input
      labels(4) = mel_x_inter ! input
      call set_rule(mel_b_dia,ttype_opme,PRECONDITIONER,
     &              labels,4,1,
     &              parameters,1,tgt_info)
      
      ! X^-1 for testing
      call add_target(mel_x_inv,ttype_opme,.true.,tgt_info)
      call set_dependency(mel_x_inv,op_diar12,tgt_info)
      call set_dependency(mel_x_inv,eval_r12_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_x_inv
      labels(2) = op_diar12
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_x_inv,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1) = mel_x_inv   ! output
      labels(2) = mel_x_inter ! input
      call set_rule(mel_x_inv,ttype_opme,INVERT,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      

      
*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      ! test
      call add_target(eval_r12_inter,ttype_gen,.false.,tgt_info)
      call set_dependency(eval_r12_inter,mel_ham,tgt_info)
      call set_dependency(eval_r12_inter,mel_rint,tgt_info)
      call set_dependency(eval_r12_inter,mel_rinba,tgt_info)
      call set_dependency(eval_r12_inter,mel_ttr,tgt_info)
      call set_dependency(eval_r12_inter,mel_f2,tgt_info)
      call set_dependency(eval_r12_inter,mel_v_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_vbar_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_x_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_b_def,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_vcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_vbcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_xcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_bcabs,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_vcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = fopt_r12_vbcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = fopt_r12_xcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = fopt_r12_bcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      return

      contains

      subroutine set_fvirt_constr(hpvx_constr,gas_constr)

      implicit none

      integer, intent(out) ::
     &     hpvx_constr(2,ngastp,2),
     &     gas_constr(2,orb_info%ngas,2,2)
      
      integer ::
     &     idx

      hpvx_constr(1:2,1:ngastp,1:2) = 0
      gas_constr(1:2,1:orb_info%ngas,1:2,1:2)  = 0

      print *,'1: ',hpvx_constr

      do idx = 1, ngastp
        if (idx.eq.IHOLE.or.idx.eq.IVALE) cycle
        hpvx_constr(2,idx,1:2) = 1
      end do

      do idx = 1, orb_info%ngas
        if (orb_info%ihpvgas(idx).eq.IHOLE.or.
     &      orb_info%ihpvgas(idx).eq.IVALE) cycle
        gas_constr(2,idx,1:2,1) = 1
      end do

      end subroutine

      end
