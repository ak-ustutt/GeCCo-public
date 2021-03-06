*----------------------------------------------------------------------*
      subroutine set_ccmp_general_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed in more or less all kinds of MP/CC calculations
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
      include 'ifc_targets.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     min_rank, max_rank, t1ext_mode, max_rank_guess,
     &     isim, ncat, nint, icnt, ansatz,
     &     isym, ms, msc, sym_arr(8),
     &     occ_def(ngastp,2,60), ndef
      logical ::
     &     l_exist, explicit
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters(2)

      if (iprlvl.gt.0)
     &     write(lulog,*) 'setting general targets for MP/CC ...'

      msc = 1 ! presently assuming closed shell

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! T1 transformed Hamiltonian
      ansatz=0
      explicit = is_keyword_set('method.R12').gt.0
      if (explicit.and.orb_info%caborb.gt.0)
     &     call get_argument_value('method.R12','ansatz',ival=ansatz)

      call add_target(op_hhat,ttype_op,.false.,tgt_info)
C      if (ansatz.gt.1) then
C        call hop_parameters(-1,parameters,
C     &                   0,2,2,.true.) ! 1-external at most
C        call set_rule(op_hhat,ttype_op,DEF_HAMILTONIAN,
C     &              op_hhat,1,1,
C     &              parameters,1,tgt_info)
Cc        call set_dependency(op_hhat,op_ham,tgt_info)
Cc        call cloneop_parameters(-1,parameters,
Cc     &       op_ham,.false.)
Cc        call set_rule(op_hhat,ttype_op,CLONE_OP,
Cc     &       op_hhat,1,1,
Cc     &       parameters,1,tgt_info)
C      else
        call hop_parameters(-1,parameters,
     &                   0,2,1,.false.)  ! avoid any X blocks
        call set_rule(op_hhat,ttype_op,DEF_HAMILTONIAN,
     &              op_hhat,1,1,
     &              parameters,1,tgt_info)
C      end if

      call get_argument_value('method.CC','minexc',ival=min_rank)
      call get_argument_value('method.CC','maxexc',ival=max_rank)
      if (is_keyword_set('method.ECC').gt.0) then
        call get_argument_value('method.ECC','minexc',ival=min_rank)
        call get_argument_value('method.ECC','maxexc',ival=max_rank)
      end if
      call get_argument_value('method.CC','T1ext',ival=t1ext_mode)
      if (t1ext_mode.eq.0.and.is_keyword_set('method.R12').gt.0) then
        call get_argument_value('method.R12','T1ext',ival=t1ext_mode)
      end if
      if (t1ext_mode.eq.0.and.is_keyword_set('method.ECC').gt.0) then
        call get_argument_value('method.ECC','T1ext',ival=t1ext_mode)
      end if
      call get_argument_value('method.CC','maxexc_guess',
     &                        ival=max_rank_guess)
      call get_argument_value('calculate.solve.non_linear','restart',
     &     lval=l_exist)
      if (max_rank_guess.gt.0.and..not.l_exist)
     &   call quit(1,'set_ccmp_general_targets',
     &             'Use of initial guess for T requires restart=T')

      ! T operator
      call add_target(op_top,ttype_op,.false.,tgt_info)
      if (t1ext_mode.eq.0) then
        call xop_parameters(-1,parameters,
     &                   .false.,min_rank,max_rank,0,1)
        call set_rule(op_top,ttype_op,DEF_EXCITATION,
     &              op_top,1,1,
     &              parameters,1,tgt_info)
      else
        if (min_rank.ne.1.or.max_rank.gt.2)
     &       call quit(1,'set_ccmp_general_targets',
     &       'experimental T1ext option is for CCSD only')
        occ_def = 0
        ! 1 -- T1
        occ_def(IPART,1,1) = 1
        occ_def(IHOLE,2,1) = 1
        ! 2 -- T1'
        occ_def(IEXTR,1,2) = 1
        occ_def(IHOLE,2,2) = 1
        ndef = 2
        if (max_rank.eq.2) then
          ! 3 -- T2
          occ_def(IPART,1,3) = 2
          occ_def(IHOLE,2,3) = 2
          ndef = 3
        end if
        call op_from_occ_parameters(-1,parameters,2,
     &                       occ_def,ndef,1,(/     0,     0/),ndef)
        call set_rule(op_top,ttype_op,DEF_OP_FROM_OCC,
     &                op_top,1,1,
     &                parameters,2,tgt_info)

      end if

      ! Hbar intermediate
      call add_target(op_hbar,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &                    .false.,min_rank,max_rank,0,1)
      call set_rule(op_hbar,ttype_op,DEF_CC_HBAR_OP,
     &              op_hbar,1,1,
     &              parameters,1,tgt_info)
      
      ! Tbar
      call add_target(op_tbar,ttype_op,.false.,tgt_info)
      call set_dependency(op_tbar,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.true.)
      call set_rule(op_tbar,ttype_op,CLONE_OP,
     &              op_tbar,1,1,
     &              parameters,1,tgt_info)

      ! Residual
      call add_target(op_omg,ttype_op,.false.,tgt_info)
      call set_dependency(op_omg,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.false.)
      call set_rule(op_omg,ttype_op,CLONE_OP,
     &              op_omg,1,1,
     &              parameters,1,tgt_info)

      ! Diagonal
      call add_target(op_dia,ttype_op,.false.,tgt_info)
      call set_dependency(op_dia,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.false.)
      call set_rule(op_dia,ttype_op,CLONE_OP,
     &              op_dia,1,1,
     &              parameters,1,tgt_info)

      if (max_rank_guess.gt.0) then
        ! operator for initial guess
        call add_target('Tguess',ttype_op,.false.,tgt_info)
        call xop_parameters(-1,parameters,
     &                   .false.,min_rank,max_rank_guess,0,1)
        call set_rule('Tguess',ttype_op,DEF_EXCITATION,
     &                'Tguess',1,1,
     &                parameters,1,tgt_info)
      end if

*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cchhat
      labels(2) = op_hhat
      labels(3) = op_ham
      labels(4) = op_top
      call add_target(form_cchhat,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cchhat,op_hhat,tgt_info)
      call set_dependency(form_cchhat,op_ham,tgt_info)
      call set_dependency(form_cchhat,op_top,tgt_info)
      call set_rule(form_cchhat,ttype_frm,DEF_HHAT,
     &              labels,4,1,
     &              title_cchhat,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cchbar
      labels(2) = op_hbar
      labels(3) = op_ham
      labels(4) = op_top
      call add_target(form_cchbar,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cchbar,op_hbar,tgt_info)
      call set_dependency(form_cchbar,op_ham,tgt_info)
      call set_dependency(form_cchbar,op_top,tgt_info)
      call set_rule(form_cchbar,ttype_frm,DEF_CC_HBAR,
     &              labels,4,1,
     &              title_cchbar,1,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*
      ! T-list definition
      call add_target(mel_topdef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_topdef,op_top,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_top
      labels(2) = op_top
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_topdef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! OMG-list definition
      call add_target(mel_omgdef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_omgdef,op_omg,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_omg
      labels(2) = op_omg
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_omgdef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! TBAR-list definition
      call add_target(mel_tbardef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_tbardef,op_tbar,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_tbar
      labels(2) = op_tbar
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_tbardef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! HHat definition
      call add_target(mel_hhatdef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_hhatdef,op_hhat,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_hhat
      labels(2) = op_hhat
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_hhatdef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! Hbar definition
      call add_target(meldef_hbar,ttype_opme,.false.,tgt_info)
      call set_dependency(meldef_hbar,op_hbar,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_hbar
      labels(2) = op_hbar
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(meldef_hbar,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! DIA-list for all symmetries:
      do isym = 1,  orb_info%nsym
        call me_list_label(me_label,mel_dia,isym,0,0,0,.false.)
        call add_target(me_label,ttype_opme,.false.,tgt_info)
        call set_dependency(me_label,mel_ham,tgt_info)
        call set_dependency(me_label,op_dia,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = me_label
        labels(2) = op_dia
        call me_list_parameters(-1,parameters,
     &       0,0,isym,0,0,.false.)
        call set_rule(me_label,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
        labels(1) = me_label
        labels(2) = mel_ham
        call set_rule(me_label,ttype_opme,PRECONDITIONER,
     &              labels,2,1,
     &              parameters,0,tgt_info)
      end do
      
      if (max_rank_guess.gt.0) then
        ! read in initial guess
        call add_target('DEF_ME_Tguess',ttype_opme,.true.,tgt_info)
        call set_dependency('DEF_ME_Tguess',mel_topdef,tgt_info)
        call set_dependency('DEF_ME_Tguess','Tguess',tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = 'ME_Tguess'
        labels(2) = 'Tguess'
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule('DEF_ME_Tguess',ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        inquire(file='ME_Tguess_list.da',exist=l_exist)
        if (.not.l_exist) call quit(1,'set_ccmp_general_targets',
     &           'Amplitude file for initial guess not found!')
        ! first initialize T op with zeros
        call set_rule2('DEF_ME_Tguess',UNITY,tgt_info)
        call set_arg('DEF_ME_Tguess',UNITY,'LIST',1,tgt_info,
     &               val_label=(/mel_top/))
        call set_arg('DEF_ME_Tguess',UNITY,'FAC',1,tgt_info,
     &               val_rl8=(/0d0/)) ! not needed since no diag. blks
        call set_arg('DEF_ME_Tguess',UNITY,'INIT',1,tgt_info,
     &               val_log=(/.true./))
        ! now copy the initial guess
        call set_rule2('DEF_ME_Tguess',SCALE_COPY,tgt_info)
        call set_arg('DEF_ME_Tguess',SCALE_COPY,'LIST_RES',1,tgt_info,
     &               val_label=(/mel_top/))
        call set_arg('DEF_ME_Tguess',SCALE_COPY,'LIST_INP',1,tgt_info,
     &               val_label=(/'ME_Tguess'/))
        call set_arg('DEF_ME_Tguess',SCALE_COPY,'FAC',1,tgt_info,
     &               val_rl8=(/1d0/))
      end if

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*

      return
      end
