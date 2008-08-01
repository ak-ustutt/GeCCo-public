*----------------------------------------------------------------------*
      subroutine set_general_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
*     set targets needed in more or less all kinds of CC calculations
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
      character(*), intent(in) ::
     &     env_type

      integer ::
     &     min_rank, max_rank,
     &     isim, ncat, nint, icnt, iformal,
     &     isym, ms, msc, sym_arr(8),trunc_type
      logical ::
     &     needed, explicit, truncate
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting general targets ...'

c      call get_argument_value('method.R12','truncate',lval=truncate)
      truncate = is_keyword_set('method.truncate').gt.0
      call get_argument_value('method.truncate','trunc_type',
     &     ival=trunc_type)

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! Unit operator (purely formal)
      ! we use DEF_HAMILTONIAN with rank 0 to 0 and iformal=0
      call add_target(op_unity,ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &                   0,0,0,.false.)
      call set_rule(op_unity,ttype_op,DEF_HAMILTONIAN,
     &              op_unity,1,1,
     &              parameters,1,tgt_info)

      ! Hamiltonian
      iformal = 1
      explicit = is_keyword_set('method.R12').gt.0
      if (explicit.and.orb_info%caborb.gt.0.and.(.not.truncate
     &     .or.(truncate.and.trunc_type.gt.0)))
     &     iformal = 4
      if (explicit.and.orb_info%caborb.gt.0.and.truncate
     &     .and.trunc_type.eq.0)
     &     iformal = 3
      call add_target(op_ham,ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &                   0,2,iformal,explicit)
      call set_rule(op_ham,ttype_op,DEF_HAMILTONIAN,
     &              op_ham,1,1,
     &              parameters,1,tgt_info)

      ! Fock-Operator
      call add_target(op_fock,ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &                   1,1,1,explicit)
      call set_rule(op_fock,ttype_op,DEF_HAMILTONIAN,
     &              op_fock,1,1,
     &              parameters,1,tgt_info)


*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*
      if (is_keyword_set('calculate.check.formulae')) then
        call add_target(form_test,ttype_frm,.true.,tgt_info)
        call set_rule(form_test,ttype_frm,CHECK_FORMGEN,
     &                form_test,1,1,
     &                parameters,0,tgt_info)
      end if
      
*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! Hamilton list:
      call add_target(mel_ham,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ham,op_ham,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ham
      labels(2) = op_ham
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(mel_ham,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ham
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_ham,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)
      
*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*

      return
      end
