*----------------------------------------------------------------------*
      subroutine set_ccmp_general_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed in more or less all kinds of MP/CC calculations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
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
     &     min_rank, max_rank,
     &     isim, ncat, nint, icnt, ansatz,
     &     isym, ms, msc, sym_arr(8)
      logical ::
     &     needed, explicit
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting general targets for MP/CC ...'

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! T1 transformed Hamiltonian
      ansatz=0
      explicit = is_keyword_set('method.R12').gt.0
      if (explicit.and.orb_info%caborb.gt.0)
     &     call get_argument_value('method.R12','ansatz',ival=ansatz)

      call add_target(op_hhat,ttype_op,.false.,tgt_info)
      if (ansatz.gt.1) then
        call hop_parameters(-1,parameters,
     &                   0,2,2,.true.) ! 1-external at most
        call set_rule(op_hhat,ttype_op,DEF_HAMILTONIAN,
     &              op_hhat,1,1,
     &              parameters,1,tgt_info)
c        call set_dependency(op_hhat,op_ham,tgt_info)
c        call cloneop_parameters(-1,parameters,
c     &       op_ham,.false.)
c        call set_rule(op_hhat,ttype_op,CLONE_OP,
c     &       op_hhat,1,1,
c     &       parameters,1,tgt_info)
      else
        call hop_parameters(-1,parameters,
     &                   0,2,1,.false.)  ! avoid any X blocks
        call set_rule(op_hhat,ttype_op,DEF_HAMILTONIAN,
     &              op_hhat,1,1,
     &              parameters,1,tgt_info)
      end if

      ! T operator
      call add_target(op_top,ttype_op,.false.,tgt_info)
      call get_argument_value('method.CC','minexc',ival=min_rank)
      call get_argument_value('method.CC','maxexc',ival=max_rank)
      if (is_keyword_set('method.ECC').gt.0) then
        call get_argument_value('method.ECC','minexc',ival=min_rank)
        call get_argument_value('method.ECC','maxexc',ival=max_rank)
      end if

      ! Hbar intermediate
      call add_target(op_hbar,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &                    .false.,min_rank,max_rank,0,1)
      call set_rule(op_hbar,ttype_op,DEF_CC_HBAR_OP,
     &              op_hbar,1,1,
     &              parameters,1,tgt_info)
      
      call xop_parameters(-1,parameters,
     &                   .false.,min_rank,max_rank,0,1)
      call set_rule(op_top,ttype_op,DEF_EXCITATION,
     &              op_top,1,1,
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
     &     0,0,1,0,0)
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
     &     0,0,1,0,0)
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
     &     0,0,1,0,0)
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
     &     0,0,1,0,0)
      call set_rule(mel_hhatdef,ttype_opme,DEF_ME_LIST,
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
     &       0,0,isym,0,0)
        call set_rule(me_label,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
        labels(1) = me_label
        labels(2) = mel_ham
        call set_rule(me_label,ttype_opme,PRECONDITIONER,
     &              labels,2,1,
     &              parameters,0,tgt_info)
      end do
      
*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*

      return
      end
