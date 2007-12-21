*----------------------------------------------------------------------*
      subroutine set_cc_operator_targets(tgt_info)
*----------------------------------------------------------------------*
*     set all operator targets needed in CC calculations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info

      integer ::
     &     min_rank, max_rank
      character ::
     &     parameters*(len_command_par)

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting operator targets for CC ...'

      ! Hamiltonian
      call add_target(op_ham,ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &                   0,2,1)
      call set_rule(op_ham,ttype_op,DEF_HAMILTONIAN,
     &              op_ham,1,1,
     &              parameters,1,tgt_info)

      ! T1 transformed Hamiltonian
      call add_target(op_hhat,ttype_op,.false.,tgt_info)
      call set_dependency(op_hhat,op_ham,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_ham,.false.)
      call set_rule(op_hhat,ttype_op,CLONE_OP,
     &              op_hhat,1,1,
     &              parameters,1,tgt_info)

      ! T operator
      call add_target(op_top,ttype_op,.false.,tgt_info)
      call get_argument_value('method.CC','minexc',ival=min_rank)
      call get_argument_value('method.CC','maxexc',ival=max_rank)
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

      ! Lagrange functional
      call add_target(op_cclg,ttype_op,.false.,tgt_info)
      call set_rule(op_cclg,ttype_op,DEF_SCALAR,
     &              op_cclg,1,1,
     &              parameters,0,tgt_info)
      
      ! Energy
      call add_target(op_ccen,ttype_op,.false.,tgt_info)
      call set_rule(op_ccen,ttype_op,DEF_SCALAR,
     &              op_ccen,1,1,
     &              parameters,0,tgt_info)

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


      ! Entries needed for left-hand equation:
      ! RHS vector eta
      call add_target(op_eta,ttype_op,.false.,tgt_info)
      call set_dependency(op_eta,op_tbar,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tbar,.false.)
      call set_rule(op_eta,ttype_op,CLONE_OP,
     &              op_eta,1,1,
     &              parameters,1,tgt_info)
      
      ! MV-product Tbar.A
      call add_target(op_tbar_a,ttype_op,.false.,tgt_info)
      call set_dependency(op_tbar_a,op_tbar,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tbar,.false.)
      call set_rule(op_tbar_a,ttype_op,CLONE_OP,
     &              op_tbar_a,1,1,
     &              parameters,1,tgt_info)
      
      ! one-particle density
      call add_target(op_1dens,ttype_op,.false.,tgt_info)
      call dens_parameters(-1,parameters,
     &                     1,1,1)
      call set_rule(op_1dens,ttype_op,DEF_DENSITY,
     &              op_1dens,1,1,
     &              parameters,1,tgt_info)

      ! particle conserving response
      ! right-response vector
      call add_target(op_r,ttype_op,.false.,tgt_info)
      call set_dependency(op_r,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.false.)
      call set_rule(op_r,ttype_op,CLONE_OP,
     &              op_r,1,1,
     &              parameters,1,tgt_info)
      
      ! right-response vector times Jacobian
      call add_target(op_a_r,ttype_op,.false.,tgt_info)
      call set_dependency(op_a_r,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.false.)
      call set_rule(op_a_r,ttype_op,CLONE_OP,
     &              op_a_r,1,1,
     &              parameters,1,tgt_info)

      ! left-response vector
      call add_target(op_l,ttype_op,.false.,tgt_info)
      call set_dependency(op_l,op_tbar,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tbar,.false.)
      call set_rule(op_l,ttype_op,CLONE_OP,
     &              op_l,1,1,
     &              parameters,1,tgt_info)

      ! left-response vector times Jacobian
      call add_target(op_l_a,ttype_op,.false.,tgt_info)
      call set_dependency(op_l_a,op_tbar,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tbar,.false.)
      call set_rule(op_l_a,ttype_op,CLONE_OP,
     &              op_l_a,1,1,
     &              parameters,1,tgt_info)

      ! particle removing response
      ! right-response vector
      call add_target(op_rip,ttype_op,.false.,tgt_info)
      call get_argument_value('method.CC','minexc',ival=min_rank)
      call get_argument_value('method.CC','maxexc',ival=max_rank)
      call xop_parameters(-1,parameters,
     &                   .false.,min_rank,max_rank,-1,1)
      call set_rule(op_rip,ttype_op,DEF_EXCITATION,
     &              op_rip,1,1,
     &              parameters,1,tgt_info)

      ! the appropriate diagonal:
      call add_target(op_dia_ip,ttype_op,.false.,tgt_info)
      call set_dependency(op_dia_ip,op_rip,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rip,.false.)
      call set_rule(op_dia_ip,ttype_op,CLONE_OP,
     &              op_dia_ip,1,1,
     &              parameters,1,tgt_info)
      
      ! right-response vector times Jacobian
      call add_target(op_a_rip,ttype_op,.false.,tgt_info)
      call set_dependency(op_a_rip,op_rip,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rip,.false.)
      call set_rule(op_a_rip,ttype_op,CLONE_OP,
     &              op_a_rip,1,1,
     &              parameters,1,tgt_info)


      return
      end
