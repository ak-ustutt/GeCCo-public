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
     &     isym, ms, msc, sym_arr(8),
     &     occ_def(ngastp,2,20)
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

*----------------------------------------------------------------------*
      if (iprlvl.gt.0)
     &     write(luout,*) 'setting general targets for R12 ...'

*----------------------------------------------------------------------*
*     read input
*----------------------------------------------------------------------*
      ! set approx string
      approx(1:12) = ' '
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! the formal R12 geminal: P12 r12|0>
      call add_target(op_r12,ttype_op,.false.,tgt_info)
      min_rank = 2  ! 1 is a possibility 
      call r12gem_parameters(-1,parameters,
     &                   .false.,min_rank,ansatz)
      call set_rule(op_r12,ttype_op,DEF_R12GEMINAL,
     &              op_r12,1,1,
     &              parameters,1,tgt_info)

c      ! the adjoint (should be obsolete soon)
c      call add_target(op_rba,ttype_op,.false.,tgt_info)
c      call set_dependency(op_rba,op_r12,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_r12,.true.) ! <- dagger=.true.
c      call set_rule(op_rba,ttype_op,CLONE_OP,
c     &              op_rba,1,1,
c     &              parameters,1,tgt_info)

      ! the coefficients
      call add_target(op_c12,ttype_op,.false.,tgt_info)
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
      call r12int_parameters(-1,parameters,
     &     .false.,min_rank,2,0,2)
      call set_rule(op_rint,ttype_op,DEF_R12INT,
     &              op_rint,1,1,
     &              parameters,1,tgt_info)
      
      ! (pq)_frozen/(pq)_ae block of 2e-Hamiltonian
      call add_target(op_g_x,ttype_op,.false.,tgt_info)
      min_rank = 2 
      call r12int_parameters(-1,parameters,
     &     .true.,min_rank,2,0,2)
      call set_rule(op_g_x,ttype_op,DEF_R12INT,
     &              op_g_x,1,1,
     &              parameters,1,tgt_info)

      ! soon obsolete: adjoint of the above ints ...
      call add_target(op_rinba,ttype_op,.false.,tgt_info)
      call set_dependency(op_rinba,op_rint,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rint,.true.) ! <- dagger=.true.
      call set_rule(op_rinba,ttype_op,CLONE_OP,
     &              op_rinba,1,1,
     &              parameters,1,tgt_info)

      ! extended list of R12 integrals
      call add_target(op_rintx,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! 1
      occ_def(IEXTR,1,1) = 2
      occ_def(IHOLE,2,1) = 2
      ! 2
      occ_def(IHOLE,1,2) = 2
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1
      ! 3
      occ_def(IHOLE,1,3) = 1
      occ_def(IPART,1,3) = 1
      occ_def(IHOLE,2,3) = 1
      occ_def(IPART,2,3) = 1
      ! 4
      occ_def(IHOLE,1,4) = 1
      occ_def(IEXTR,1,4) = 1
      occ_def(IHOLE,2,4) = 1
      occ_def(IPART,2,4) = 1
      ! 5
      occ_def(IPART,1,5) = 2
      occ_def(IHOLE,2,5) = 1
      occ_def(IPART,2,5) = 1
      ! 6
      occ_def(IPART,1,6) = 1
      occ_def(IEXTR,1,6) = 1
      occ_def(IHOLE,2,6) = 1
      occ_def(IPART,2,6) = 1
      ! 7
      occ_def(IHOLE,1,7) = 2
      occ_def(IHOLE,2,7) = 1
      occ_def(IEXTR,2,7) = 1
      ! 8
      occ_def(IHOLE,1,8) = 1
      occ_def(IPART,1,8) = 1
      occ_def(IHOLE,2,8) = 1
      occ_def(IEXTR,2,8) = 1
      ! 9
      occ_def(IHOLE,1,9) = 1
      occ_def(IEXTR,1,9) = 1
      occ_def(IHOLE,2,9) = 1
      occ_def(IEXTR,2,9) = 1
      ! 10
      occ_def(IPART,1,10) = 2
      occ_def(IHOLE,2,10) = 1
      occ_def(IEXTR,2,10) = 1
      ! 11
      occ_def(IPART,1,11) = 1
      occ_def(IEXTR,1,11) = 1
      occ_def(IHOLE,2,11) = 1
      occ_def(IEXTR,2,11) = 1
      call op_from_occ_parameters(-1,parameters,2,
     &       occ_def,11,1,11)
      call set_rule(op_rintx,ttype_op,DEF_OP_FROM_OCC,
     &              op_rintx,1,1,
     &              parameters,2,tgt_info)
      
      ! commutator integrals <ab|[T1+T2,r12]|ij>
      call add_target(op_ttr,ttype_op,.false.,tgt_info)
      call set_dependency(op_ttr,op_rint,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rint,.false.) ! <- dagger=.false.
      call set_rule(op_ttr,ttype_op,CLONE_OP,
     &              op_ttr,1,1,
     &              parameters,1,tgt_info)
      
      ! (f+k) modified integrals r12bar
      call add_target(op_rintbar,ttype_op,.false.,tgt_info)
      call set_dependency(op_rintbar,op_rint,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rint,.false.) 
      call set_rule(op_rintbar,ttype_op,CLONE_OP,
     &              op_rintbar,1,1,
     &              parameters,1,tgt_info)
      
      ! k modified integrals r12tilde
      call add_target(op_rinttilde,ttype_op,.false.,tgt_info)
      call set_dependency(op_rinttilde,op_rint,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rint,.false.) 
      call set_rule(op_rinttilde,ttype_op,CLONE_OP,
     &              op_rinttilde,1,1,
     &              parameters,1,tgt_info)
      
      ! commutator integrals <kl|r12[T1+T2,r12]|ij>
      call add_target(op_rttr,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,2,2,0,2)
      call set_rule(op_rttr,ttype_op,DEF_R12INTERM,
     &              op_rttr,1,1,
     &              parameters,1,tgt_info)


c      ! the adjoint
c      call add_target(op_ttr_bar,ttype_op,.false.,tgt_info)
c      call set_dependency(op_ttr_bar,op_ttr,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_ttr,.true.) ! <- dagger=.true.
c      call set_rule(op_ttr_bar,ttype_op,CLONE_OP,
c     &              op_ttr_bar,1,1,
c     &              parameters,1,tgt_info)

      ! (G.R)^{ij}_{pq}
      call add_target(op_gr,ttype_op,.false.,tgt_info)
      call set_dependency(op_gr,op_v_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_v_inter,.false.) ! <- dagger=.false.
      call set_rule(op_gr,ttype_op,CLONE_OP,
     &              op_gr,1,1,
     &              parameters,1,tgt_info)
            
      ! (R.G)^{pq}_{ij}
      call add_target(op_rg,ttype_op,.false.,tgt_info)
      call set_dependency(op_rg,op_v_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_v_inter,.true.) ! <- dagger=.true.
      call set_rule(op_rg,ttype_op,CLONE_OP,
     &              op_rg,1,1,
     &              parameters,1,tgt_info)
            
      ! V^{ij}_{pq}
      call add_target(op_v_inter,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,2,2,0,2)
      call set_rule(op_v_inter,ttype_op,DEF_R12INTERM,
     &              op_v_inter,1,1,
     &              parameters,1,tgt_info)
      
      ! the adjoint
      call add_target(op_vbar_inter,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,2,2,0,2)
      call set_rule(op_vbar_inter,ttype_op,DEF_R12INTERM,
     &              op_vbar_inter,1,1,
     &              parameters,1,tgt_info)
c      call set_dependency(op_vbar_inter,op_v_inter,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_v_inter,.true.) ! <- dagger=.false.
c                                   ! we really need the transposed list
c      call set_rule(op_vbar_inter,ttype_op,CLONE_OP,
c     &              op_vbar_inter,1,1,
c     &              parameters,1,tgt_info)
      
      ! B intermediate
      call add_target(op_b_inter,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,2,2,0,2)
      call set_rule(op_b_inter,ttype_op,DEF_R12INTERM,
     &              op_b_inter,1,1,
     &              parameters,1,tgt_info)

c      call add_target(op_b_inter,ttype_op,.false.,tgt_info)
c      call set_dependency(op_b_inter,op_v_inter,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_v_inter,.false.) ! <- dagger=.false.
c      call set_rule(op_b_inter,ttype_op,CLONE_OP,
c     &              op_b_inter,1,1,
c     &              parameters,1,tgt_info)

      ! R12^{2} integrals
      call add_target(op_f2,ttype_op,.false.,tgt_info)
      if (approx(1:1).eq.'A') then
        call set_dependency(op_f2,op_b_inter,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
        call set_rule(op_f2,ttype_op,CLONE_OP,
     &              op_f2,1,1,
     &              parameters,1,tgt_info)
      else
        occ_def = 0
        ! 1
        occ_def(IHOLE,1,1) = 2
        occ_def(IHOLE,2,2) = 2
        ! 2
        occ_def(IHOLE,1,3) = 1
        occ_def(IPART,1,3) = 1
        occ_def(IHOLE,2,4) = 2
        ! 3
        occ_def(IHOLE,1,5) = 1
        occ_def(IEXTR,1,5) = 1
        occ_def(IHOLE,2,6) = 2
        call op_from_occ_parameters(-1,parameters,2,
     &       occ_def,3,2,6)
        call set_rule(op_f2,ttype_op,DEF_OP_FROM_OCC,
     &                op_f2,1,1,
     &                parameters,2,tgt_info)
      end if

      ! X intermediate
      call add_target(op_x_inter,ttype_op,.false.,tgt_info)
      call set_dependency(op_x_inter,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
      call set_rule(op_x_inter,ttype_op,CLONE_OP,
     &              op_x_inter,1,1,
     &              parameters,1,tgt_info)

c      ! Xbar = X(f+k) intermediate
c      call add_target(op_xbar_inter,ttype_op,.false.,tgt_info)
c      call set_dependency(op_xbar_inter,op_b_inter,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_b_inter,.false.) ! <- dagger=.false.
c      call set_rule(op_xbar_inter,ttype_op,CLONE_OP,
c     &              op_xbar_inter,1,1,
c     &              parameters,1,tgt_info)

c      ! Xtilde = X k intermediate
c      call add_target(op_xtilde_inter,ttype_op,.false.,tgt_info)
c      call set_dependency(op_xtilde_inter,op_b_inter,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_b_inter,.false.) ! <- dagger=.false.
c      call set_rule(op_xtilde_inter,ttype_op,CLONE_OP,
c     &              op_xtilde_inter,1,1,
c     &              parameters,1,tgt_info)

      ! C intermediate
      call add_target(op_c_inter,ttype_op,.false.,tgt_info)
      occ_def = 0
      occ_def(IPART,1,1) = 2
      occ_def(IHOLE,2,1) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,1,1,1)
      call set_rule(op_c_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_c_inter,1,1,
     &              parameters,2,tgt_info)

      ! Cbar intermediate
      call add_target(op_cbar_inter,ttype_op,.false.,tgt_info)
      occ_def = 0
      occ_def(IHOLE,1,1) = 2
      occ_def(IPART,2,1) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,1,1,1)
      call set_rule(op_cbar_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_cbar_inter,1,1,
     &              parameters,2,tgt_info)

      ! P intermediate
      call add_target(op_p_inter,ttype_op,.false.,tgt_info)
      call set_dependency(op_p_inter,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
      call set_rule(op_p_inter,ttype_op,CLONE_OP,
     &              op_p_inter,1,1,
     &              parameters,1,tgt_info)

      ! Z intermediate
      call add_target(op_z_inter,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,1,2) = 1
      occ_def(IHOLE,2,2) = 1
      occ_def(IHOLE,2,3) = 2
      ! 2
      occ_def(IHOLE,1,4) = 2
      occ_def(IHOLE,1,5) = 1
      occ_def(IPART,2,5) = 1
      occ_def(IHOLE,2,6) = 2
      ! 3
      occ_def(IHOLE,1,7) = 2
      occ_def(IPART,1,8) = 1
      occ_def(IHOLE,2,8) = 1
      occ_def(IHOLE,2,9) = 2
      ! 4
      occ_def(IHOLE,1,10) = 2
      occ_def(IPART,1,11) = 1
      occ_def(IPART,2,11) = 1
      occ_def(IHOLE,2,12) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,4,3,12)
      call set_rule(op_z_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_z_inter,1,1,
     &              parameters,2,tgt_info)

      ! inverse of B
      call add_target(op_b_inv,ttype_op,.false.,tgt_info)
      call set_dependency(op_b_inv,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
      call set_rule(op_b_inv,ttype_op,CLONE_OP,
     &              op_b_inv,1,1,
     &              parameters,1,tgt_info)

      ! inverse of X
      call add_target(op_x_inv,ttype_op,.false.,tgt_info)
      call set_dependency(op_x_inv,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
      call set_rule(op_x_inv,ttype_op,CLONE_OP,
     &              op_x_inv,1,1,
     &              parameters,1,tgt_info)

      ! Exchange operator, K.
      call add_target(op_exchange,ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &     1,1,2,.true.)
      call set_rule(op_exchange,ttype_op,DEF_HAMILTONIAN,
     &              op_exchange,1,1,
     &              parameters,1,tgt_info)
      
      ! Hartree operator, F+K.
      call add_target(op_hartree,ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &     1,1,2,.true.)
      call set_rule(op_hartree,ttype_op,DEF_HAMILTONIAN,
     &              op_hartree,1,1,
     &              parameters,1,tgt_info)
      

*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*
      ! R-bar intermediate
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_r12bar
      labels(2) = op_rintbar
      labels(3) = op_rint
      labels(4) = op_rintx
      labels(5) = op_hartree
      ! for GBC: pass op_exchange instead
      call add_target(form_r12_r12bar,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_r12bar,op_rintbar,tgt_info)
      call set_dependency(form_r12_r12bar,op_rint,tgt_info)
      call set_dependency(form_r12_r12bar,op_rintx,tgt_info)
      call set_dependency(form_r12_r12bar,op_hartree,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_rbar,ansatz,'RB'//approx)
      call set_rule(form_r12_r12bar,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      
      ! R-tilde intermediate
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_r12tilde
      labels(2) = op_rinttilde
      labels(3) = op_rint
      labels(4) = op_rintx
      labels(5) = op_exchange
      call add_target(form_r12_r12tilde,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_r12tilde,op_rinttilde,tgt_info)
      call set_dependency(form_r12_r12tilde,op_rint,tgt_info)
      call set_dependency(form_r12_r12tilde,op_rintx,tgt_info)
      call set_dependency(form_r12_r12tilde,op_exchange,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_rtilde,ansatz,'RT'//approx)
      call set_rule(form_r12_r12tilde,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)

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
c     &     parameters,2,title_r12_vint,0,'gxxr')
      call set_rule(form_r12_vint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to V
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_vcabs
      labels(2) = op_v_inter
      labels(3) = op_g_x !op_ham
      labels(4) = op_rint
      labels(5) = op_unity
      ! F12: op_gr
      call add_target(form_r12_vcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_vcabs,op_v_inter,tgt_info)
      call set_dependency(form_r12_vcabs,op_unity,tgt_info)
c      call set_dependency(form_r12_vcabs,op_ham,tgt_info)
      call set_dependency(form_r12_vcabs,op_g_x,tgt_info)
      call set_dependency(form_r12_vcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vcabs,ansatz,'V '//approx)
      call set_rule(form_r12_vcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! to be removed soon: transpose of V
      ! formal definition of Vbar (viz. V^dagger)
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_r12_vbint
c      labels(2) = op_vbar_inter
c      labels(3) = op_r12
c      labels(4) = op_ham
c      call add_target(form_r12_vbint,ttype_frm,.false.,tgt_info)
c      call set_dependency(form_r12_vbint,op_vbar_inter,tgt_info)
c      call set_dependency(form_r12_vbint,op_ham,tgt_info)
c      call set_dependency(form_r12_vbint,op_rba,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_r12_vbint,0,'rxg')
cc     &     parameters,2,title_r12_vbint,0,'rxxg')
c      call set_rule(form_r12_vbint,ttype_frm,DEF_R12INTM_FORMAL,
c     &              labels,4,1,
c     &              parameters,2,tgt_info)

c      ! CABS approximation to Vbar
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_r12_vbcabs
c      labels(2) = op_vbar_inter
c      labels(3) = op_rinba
c      labels(4) = op_ham
cc      labels(5) = op_unity
c      call add_target(form_r12_vbcabs,ttype_frm,.false.,tgt_info)
c      call set_dependency(form_r12_vbcabs,op_vbar_inter,tgt_info)
c      call set_dependency(form_r12_vbcabs,op_unity,tgt_info)
c      call set_dependency(form_r12_vbcabs,op_ham,tgt_info)
c      call set_dependency(form_r12_vbcabs,op_rinba,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_r12_vbcabs,ansatz,'V+'//approx)
c      call set_rule(form_r12_vbcabs,ttype_frm,DEF_R12INTM_CABS,
c     &              labels,5,1,
c     &              parameters,2,tgt_info)

      ! formal definition of X
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xint
      labels(2) = op_x_inter
c      labels(3) = op_rba
      labels(3) = op_r12
      labels(4) = op_r12
      call add_target(form_r12_xint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xint,op_x_inter,tgt_info)
c      call set_dependency(form_r12_xint,op_rba,tgt_info)
      call set_dependency(form_r12_xint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xint,0,'rxr')
c     &     parameters,2,title_r12_xint,0,'rxxr')
      call set_rule(form_r12_xint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to X
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xcabs
      labels(2) = op_x_inter
      labels(3) = op_rint
      labels(4) = op_rint
      labels(5) = op_f2
      call add_target(form_r12_xcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xcabs,op_x_inter,tgt_info)
      call set_dependency(form_r12_xcabs,op_f2,tgt_info)
c      call set_dependency(form_r12_xcabs,op_rinba,tgt_info)
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
c      labels(3) = op_rba
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_r12_bint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_bint,op_b_inter,tgt_info)
c      call set_dependency(form_r12_bint,op_rba,tgt_info)
      call set_dependency(form_r12_bint,op_ham,tgt_info)
      call set_dependency(form_r12_bint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_bint,0,'rfxr')
c     &     parameters,2,title_r12_bint,0,'rfxxr')
      call set_rule(form_r12_bint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)


      ! CABS approximation to B
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_bcabs
      labels(2) = op_b_inter
      labels(3) = op_rint
      labels(4) = op_ttr
      labels(5) = op_unity
      nlab = 5
      call add_target(form_r12_bcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_bcabs,op_b_inter,tgt_info)
      call set_dependency(form_r12_bcabs,op_unity,tgt_info)
c      call set_dependency(form_r12_bcabs,op_rinba,tgt_info)
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


      ! formal definition of C
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_cint
      labels(2) = op_c_inter
      labels(3) = op_ham
      labels(4) = op_r12
      call add_target(form_r12_cint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_cint,op_c_inter,tgt_info)
      call set_dependency(form_r12_cint,op_r12,tgt_info)
      call set_dependency(form_r12_cint,op_ham,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_cint,0,'fr')
      call set_rule(form_r12_cint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! formal definition of C+
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_r12_cbint
c      labels(2) = op_cbar_inter
c      labels(3) = op_rba
c      labels(4) = op_ham
c      call add_target(form_r12_cbint,ttype_frm,.true.,tgt_info)
c      call set_dependency(form_r12_cbint,op_cbar_inter,tgt_info)
c      call set_dependency(form_r12_cbint,op_rba,tgt_info)
c      call set_dependency(form_r12_cbint,op_ham,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_r12_cbint,0,'rf')
c      call set_rule(form_r12_cbint,ttype_frm,DEF_R12INTM_FORMAL,
c     &              labels,4,1,
c     &              parameters,2,tgt_info)

      ! formal definition of P
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_pint
      labels(2) = op_p_inter
c      labels(3) = op_rba
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_r12_pint,ttype_frm,.true.,tgt_info)
      call set_dependency(form_r12_pint,op_p_inter,tgt_info)
c      call set_dependency(form_r12_pint,op_rba,tgt_info)
      call set_dependency(form_r12_pint,op_ham,tgt_info)
      call set_dependency(form_r12_pint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_pint,0,'rgxr')
      call set_rule(form_r12_pint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! formal definition of Z
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_zint
      labels(2) = op_z_inter
c      labels(3) = op_rba
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_r12_zint,ttype_frm,.true.,tgt_info)
      call set_dependency(form_r12_zint,op_z_inter,tgt_info)
c      call set_dependency(form_r12_zint,op_rba,tgt_info)
      call set_dependency(form_r12_zint,op_r12,tgt_info)
      call set_dependency(form_r12_zint,op_ham,tgt_info)
      call set_dependency(form_r12_zint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_zint,0,'rxgxr')
      call set_rule(form_r12_zint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
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
      call set_dependency(fopt_r12_vcabs,mel_gintx,tgt_info)
c      call set_dependency(fopt_r12_vcabs,mel_ham,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_vcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

c      ! set V+
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = fopt_r12_vbcabs
c      labels(2) = form_r12_vbcabs
c      ncat = 1
c      nint = 0
c      call add_target(fopt_r12_vbcabs,ttype_frm,.false.,tgt_info)
c      call set_dependency(fopt_r12_vbcabs,form_r12_vbcabs,tgt_info)
c      call set_dependency(fopt_r12_vbcabs,mel_vbar_def,tgt_info)
c      call set_dependency(fopt_r12_vbcabs,mel_ham,tgt_info)
c      call set_dependency(fopt_r12_vbcabs,mel_rint,tgt_info)      
c      call opt_parameters(-1,parameters,ncat,nint)
c      call set_rule(fopt_r12_vbcabs,ttype_frm,OPTIMIZE,
c     &              labels,ncat+nint+1,1,
c     &              parameters,1,tgt_info)

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
c      call set_dependency(fopt_r12_xcabs,mel_rinba,tgt_info)
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
c      call set_dependency(fopt_r12_bcabs,mel_rinba,tgt_info)      
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

      ! R12integrals (extended list)
      call add_target(mel_rintx,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_rintx,op_rintx,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rintx
      labels(2) = op_rintx
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_rintx,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rintx
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_rintx,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! special two-electron integral list
      call add_target(mel_gintx,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_gintx,op_g_x,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_gintx
      labels(2) = op_g_x
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_gintx,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_gintx
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_gintx,ttype_opme,IMPORT,
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
      call add_target(mel_ttr,ttype_opme,.true.,tgt_info)
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

      ! R12[T1+T2,R12] intergrals (for f(R12))
      call add_target(mel_rttr,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_rttr,op_rttr,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rttr
      labels(2) = op_rttr
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_rttr,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rttr
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_rttr,ttype_opme,IMPORT,
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

      ! R12.G integrals (for f(R12))
      call add_target(mel_rg,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_rg,op_rg,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rg
      labels(2) = op_rg
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_rg,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rg
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_rg,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! Exchange integrals, K.
      call add_target(mel_exchange,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_exchange,op_exchange,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_exchange
      labels(2) = op_exchange
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_exchange,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_exchange
      call set_rule(mel_exchange,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! Hartree integrals, F+K.
      call add_target(mel_hartree,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_hartree,op_hartree,tgt_info)
      call set_dependency(mel_hartree,mel_exchange,tgt_info)
      call set_dependency(mel_hartree,mel_ham,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_hartree
      labels(2) = op_hartree
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_hartree,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) calculate
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_hartree
      labels(2) = mel_ham
      labels(3) = mel_exchange
      call add_parameters(-1,parameters,2,(/1d0,1d0/),2)
      call set_rule(mel_hartree,ttype_opme,ADD,
     &              labels,3,1,
     &              parameters,1,tgt_info)

      ! ----------------------------------------
      ! B) definition of lists for intermediates
      ! ----------------------------------------
      ! RBAR-list
      call add_target(meldef_rintbar,ttype_opme,.false.,tgt_info)
      call set_dependency(meldef_rintbar,op_rintbar,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rintbar
      labels(2) = op_rintbar
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(meldef_rintbar,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! RTILDE-list
      call add_target(meldef_rinttilde,ttype_opme,.false.,tgt_info)
      call set_dependency(meldef_rinttilde,op_rinttilde,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rinttilde
      labels(2) = op_rinttilde
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(meldef_rinttilde,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

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

c      ! V+-list
c      call add_target(mel_vbar_def,ttype_opme,.false.,tgt_info)
c      call set_dependency(mel_vbar_def,op_vbar_inter,tgt_info)
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_vbar_inter
c      labels(2) = op_vbar_inter
c      call me_list_parameters(-1,parameters,
c     &     0,0,1,0,0)
c      call set_rule(mel_vbar_def,ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)


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
      call add_target(mel_b_inv,ttype_opme,.false.,tgt_info)
c      call set_dependency(mel_b_inv,op_diar12,tgt_info)
      call set_dependency(mel_b_inv,op_b_inv,tgt_info)
      call set_dependency(mel_b_inv,eval_r12_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_b_inv
c      labels(2) = op_diar12
      labels(2) = op_b_inv ! actually, B^-1 should have the 
c                             ! contravariant shape
c                             ! but as long as we do not formally 
c                             ! calculate with
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
      call add_target(mel_b_dia,ttype_opme,.false.,tgt_info)
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
      call add_target(mel_x_inv,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_x_inv,op_diar12,tgt_info)
      call set_dependency(mel_x_inv,eval_r12_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_x_inv
      labels(2) = op_x_inter
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
      call set_dependency(eval_r12_inter,mel_gintx,tgt_info)
c      call set_dependency(eval_r12_inter,mel_rinba,tgt_info)
      call set_dependency(eval_r12_inter,mel_ttr,tgt_info)
      call set_dependency(eval_r12_inter,mel_f2,tgt_info)
      call set_dependency(eval_r12_inter,mel_v_def,tgt_info)
c      call set_dependency(eval_r12_inter,mel_vbar_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_x_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_b_def,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_vcabs,tgt_info)
c      call set_dependency(eval_r12_inter,fopt_r12_vbcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_xcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_bcabs,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_vcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
c      labels(1) = fopt_r12_vbcabs
c      call set_rule(eval_r12_inter,ttype_opme,EVAL,
c     &     labels,1,0,
c     &     parameters,0,tgt_info)
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
        if (orb_info%ihpvgas(idx,1).eq.IHOLE.or.
     &      orb_info%ihpvgas(idx,1).eq.IVALE) cycle
        gas_constr(2,idx,1:2,1) = 1
      end do

      end subroutine

      end
