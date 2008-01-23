*----------------------------------------------------------------------*
      subroutine set_ccr12_special_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed specifically in CC-R12 calculations
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
     &     min_rank, max_rank, ansatz,
     &     isim, ncat, nint, icnt,
     &     isym, ms, msc, sym_arr(8)
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters(2)

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting special targets for MP-R12 ...'

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! Lagrange functional
      call add_target(op_ccr12lg,ttype_op,.false.,tgt_info)
      call set_rule(op_ccr12lg,ttype_op,DEF_SCALAR,
     &              op_ccr12lg,1,1,
     &              parameters,0,tgt_info)
      
      ! Energy
      call add_target(op_ccr12en,ttype_op,.false.,tgt_info)
      call set_rule(op_ccr12en,ttype_op,DEF_SCALAR,
     &              op_ccr12en,1,1,
     &              parameters,0,tgt_info)

      ! residual
      call add_target(op_omgr12,ttype_op,.false.,tgt_info)
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call xop_parameters(-1,parameters,
     &     .false.,min_rank,max_rank,0,2)
      call set_rule(op_omgr12,ttype_op,DEF_R12INTERM,
     &              op_omgr12,1,1,
     &              parameters,1,tgt_info)

      ! diagonal
      call add_target(op_diar12,ttype_op,.false.,tgt_info)
      call set_dependency(op_diar12,op_omgr12,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_omgr12,.false.) ! <- dagger=.false.
      call set_rule(op_diar12,ttype_op,CLONE_OP,
     &              op_diar12,1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*
      call get_argument_value('method.R12','ansatz',ival=ansatz)

      call add_target(form_ccr12lg0,ttype_frm,.false.,tgt_info)
      ! (a) set formal Lagrangian (in 'complete' basis)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12lg0
      labels(2) = op_ccr12lg
      labels(3) = op_ham
      labels(4) = op_r12
      labels(5) = op_rba
      labels(6) = op_tbar
      labels(7) = op_cba
      labels(8) = op_top
      labels(9) = op_c12
      call set_dependency(form_ccr12lg0,op_ccr12lg,tgt_info)
      call set_dependency(form_ccr12lg0,op_ham,tgt_info)
      call set_dependency(form_ccr12lg0,op_r12,tgt_info)
      call set_dependency(form_ccr12lg0,op_rba,tgt_info)
      call set_dependency(form_ccr12lg0,op_tbar,tgt_info)
      call set_dependency(form_ccr12lg0,op_top,tgt_info)
      call set_dependency(form_ccr12lg0,op_cba,tgt_info)
      call set_dependency(form_ccr12lg0,op_c12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_ccr12lg0,ansatz,'---')
      call set_rule(form_ccr12lg0,ttype_frm,DEF_CCR12_LAGRANGIAN,
     &              labels,9,1,
     &              parameters,2,tgt_info)
      ! (b) Factor out the R12 intermediates 
      ! (effectively removing all reference to the complete basis)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12lg0 ! output formula (itself)
      labels(2) = form_ccr12lg0 ! input formula
      labels(3) = form_r12_vint    ! the intermediates to be factored
      labels(4) = form_r12_vbint
      labels(5) = form_r12_xint
      labels(6) = form_r12_bint
      call set_dependency(form_ccr12lg0,form_r12_vint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_vbint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_xint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_bint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_ccr12lg0,4,'---')
      call set_rule(form_ccr12lg0,ttype_frm,FACTOR_OUT,
     &              labels,6,1,
     &              parameters,2,tgt_info)
      ! (c) post-processing: remove terms which do not contribute for
      !     the given R12-approximation
      ! .... to come

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12en0
      labels(2) = form_ccr12lg0
      labels(3) = op_ccr12en
      labels(4) = op_tbar
      labels(5) = op_cba
      call add_target(form_ccr12en0,ttype_frm,.true.,tgt_info)
      call set_dependency(form_ccr12en0,form_ccr12lg0,tgt_info)
      call set_dependency(form_ccr12en0,op_ccr12en,tgt_info)
      call set_rule(form_ccr12en0,ttype_frm,INVARIANT,
     &              labels,5,1,
     &              title_ccr12en0,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12rs_t
      labels(2) = form_ccr12lg0
      labels(3) = op_omg
      labels(4) = op_tbar
      labels(5) = ' '
      call add_target(form_ccr12rs_t,ttype_frm,.true.,tgt_info)
      call set_dependency(form_ccr12rs_t,form_ccr12lg0,tgt_info)
      call set_dependency(form_ccr12rs_t,op_omg,tgt_info)
      call set_rule(form_ccr12rs_t,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_ccr12rs_t,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12rs_c
      labels(2) = form_ccr12lg0
      labels(3) = op_omgr12
      labels(4) = op_cba
      labels(5) = ' '
      call add_target(form_ccr12rs_c,ttype_frm,.true.,tgt_info)
      call set_dependency(form_ccr12rs_c,form_ccr12lg0,tgt_info)
      call set_dependency(form_ccr12rs_c,op_omgr12,tgt_info)
      call set_rule(form_ccr12rs_c,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_ccr12rs_c,1,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*
      
*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*

      return
      end
