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
     &     needed, r12fix
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character(12) ::
     &     approx, F_appr, K_appr

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
      F_appr(1:12) = ' '
      K_appr(1:12) = ' '
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','F_appr',str=F_appr)
      call get_argument_value('method.R12','K_appr',str=K_appr)
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.R12','fixed',lval=r12fix)

      ! assemble approx string
      select case(trim(F_appr))
      case('none')
        write(luout,*) 'no approximations wrt. Fock made'
      case('no_Z')
        write(luout,*) 'Z matrix omitted'
        approx(4:6) = 'noZ'
      case('GBC','EBC')
        write(luout,*)
     &  'GBC/EBC are currently only possible be supplying the'
        write(luout,*)
     &  'suitable integrals. Make that sure and restart w/o'
        write(luout,*)
     &  'GBC/EBC flag'
        call quit(0,'set_r12_general_targets','GBC/EBC?')
      case default
        call quit(0,'set_r12_general_targets',
     &       'F_appr unknown: "'//trim(F_appr)//'"')
      end select

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

      ! Only need coefficients if optimising the R12 contribution.
c      if(.not.r12fix)then
        ! the coefficients
        call add_target(op_c12,ttype_op,.false.,tgt_info)
c        if(r12fix)then
c          call set_dependency(op_c12,op_unity,tgt_info)
c        endif
        call xop_parameters(-1,parameters,
     &       .false.,min_rank,max_rank,0,max_rank+1)
        call set_rule(op_c12,ttype_op,DEF_R12COEFF,
     &                op_c12,1,1,
     &                parameters,1,tgt_info)

        ! Lagrange multipliers associated with coefficients
        call add_target(op_cba,ttype_op,.false.,tgt_info)
        call set_dependency(op_cba,op_c12,tgt_info)
c        if(r12fix)then
c          call set_dependency(op_cba,op_unity,tgt_info)
c        endif
        call cloneop_parameters(-1,parameters,
     &                          op_c12,.true.) ! <- dagger=.true.
        call set_rule(op_cba,ttype_op,CLONE_OP,
     &                op_cba,1,1,
     &                parameters,1,tgt_info)

      if(.not.r12fix)then
        ! Preconditioner
        call add_target(op_diar12,ttype_op,.false.,tgt_info)
        call set_dependency(op_diar12,op_c12,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_c12,.false.) ! <- dagger=.false.
        call set_rule(op_diar12,ttype_op,CLONE_OP,
     &                op_diar12,1,1,
     &                parameters,1,tgt_info)
      endif

      ! Now: the operators associated with the actual R12 integrals:
      !  <pq'|r12|ij> 
      call add_target(op_rint,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
c     &     .false.,min_rank,2,0,2)
     &     .false.,min_rank,2,0,3) ! 3: two externals needed
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

      ! extended list of R12 integrals (UNUSED at present)
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
      call r12int_parameters(-1,parameters,
     &     .false.,min_rank,2,0,2)
      call set_rule(op_ttr,ttype_op,DEF_R12INT,
     &              op_ttr,1,1,
     &              parameters,1,tgt_info)
      
      ! A(f+k) modified integrals r12bar
      call add_target(op_rintbar,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
     &     .false.,min_rank,2,0,2)
      call set_rule(op_rintbar,ttype_op,DEF_R12INT,
     &              op_rintbar,1,1,
     &              parameters,1,tgt_info)
      
      ! C(f+k) modified integrals r12bar+
      call add_target(op_rdagbar,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
     &     .false.,min_rank,2,0,2)
      call set_rule(op_rdagbar,ttype_op,DEF_R12INT,
     &              op_rdagbar,1,1,
     &              parameters,1,tgt_info)
      
      ! C(f) modified integrals r12breve
      call add_target(op_rintbreve,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
     &     .false.,min_rank,2,0,2)
      call set_rule(op_rintbreve,ttype_op,DEF_R12INT,
     &              op_rintbreve,1,1,
     &              parameters,1,tgt_info)
      
      ! C k modified integrals r12tilde (cloning rint -> 2ext usually)
      call add_target(op_rinttilde,ttype_op,.false.,tgt_info)
      call set_dependency(op_rinttilde,op_rint,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rint,.false.) 
      call set_rule(op_rinttilde,ttype_op,CLONE_OP,
     &              op_rinttilde,1,1,
     &              parameters,1,tgt_info)
      
      ! C f modified integrals r12c (=C-intermediate of approx C)
c      call add_target(op_rintc,ttype_op,.false.,tgt_info)
c      call xop_parameters(-1,parameters,
c     &     .false.,2,2,0,1)
c      call set_rule(op_rintc,ttype_op,DEF_EXCITATION,
c     &              op_rintc,1,1,
c     &              parameters,1,tgt_info)
c      call set_dependency(op_rintc,op_c_inter,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_c_inter,.false.) 
c      call set_rule(op_rintc,ttype_op,CLONE_OP,
c     &              op_rintc,1,1,
c     &              parameters,1,tgt_info)
      
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

      ! V^{ij}_{pq}
      call add_target(op_v_inter,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,2,2,0,2)
      call set_rule(op_v_inter,ttype_op,DEF_R12INTERM,
     &              op_v_inter,1,1,
     &              parameters,1,tgt_info)
      
c      if(r12fix)then
c        call add_target(op_v0_inter,ttype_op,.false.,tgt_info)
c        occ_def = 0
c        call op_from_occ_parameters(-1,parameters,2,
c     &       occ_def,1,1,1)
c        call set_rule(op_v0_inter,ttype_op,DEF_OP_FROM_OCC,
c     &                op_v0_inter,1,1,
c     &                parameters,2,tgt_info)
c      endif

c      ! the adjoint
c      call add_target(op_vbar_inter,ttype_op,.false.,tgt_info)
c      call xop_parameters(-1,parameters,
c     &     .false.,2,2,0,2)
c      call set_rule(op_vbar_inter,ttype_op,DEF_R12INTERM,
c     &              op_vbar_inter,1,1,
c     &              parameters,1,tgt_info)
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

c      if(r12fix)then
c        call add_target(op_b0_inter,ttype_op,.false.,tgt_info)
c        occ_def = 0
c        call op_from_occ_parameters(-1,parameters,2,
c     &       occ_def,1,1,1)
c        call set_rule(op_b0_inter,ttype_op,DEF_OP_FROM_OCC,
c     &                op_b0_inter,1,1,
c     &                parameters,2,tgt_info)
c      endif

c      call add_target(op_b_inter,ttype_op,.false.,tgt_info)
c      call set_dependency(op_b_inter,op_v_inter,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_v_inter,.false.) ! <- dagger=.false.
c      call set_rule(op_b_inter,ttype_op,CLONE_OP,
c     &              op_b_inter,1,1,
c     &              parameters,1,tgt_info)

      ! R12^{2} integrals
      call add_target(op_ff,ttype_op,.false.,tgt_info)
c      if (approx(1:1).eq.'A') then
      call set_dependency(op_ff,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     op_b_inter,.false.)  ! <- dagger=.false.
      call set_rule(op_ff,ttype_op,CLONE_OP,
     &              op_ff,1,1,
     &              parameters,1,tgt_info)
c      else
c        occ_def = 0
c        ! 1
c        occ_def(IHOLE,1,1) = 2
c        occ_def(IHOLE,2,2) = 2
c        ! 2
c        occ_def(IHOLE,1,3) = 1
c        occ_def(IPART,1,3) = 1
c        occ_def(IHOLE,2,4) = 2
c        ! 3
c        occ_def(IHOLE,1,5) = 1
c        occ_def(IEXTR,1,5) = 1
c        occ_def(IHOLE,2,6) = 2
c        call op_from_occ_parameters(-1,parameters,2,
c     &       occ_def,3,2,6)
c        call set_rule(op_ff,ttype_op,DEF_OP_FROM_OCC,
c     &                op_ff,1,1,
c     &                parameters,2,tgt_info)
c      end if

      ! {R12^2}BAR integrals
      call add_target(op_ffbar,ttype_op,.false.,tgt_info)
      call set_dependency(op_ffbar,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
      call set_rule(op_ffbar,ttype_op,CLONE_OP,
     &              op_ffbar,1,1,
     &              parameters,1,tgt_info)

      ! X intermediate
      call add_target(op_x_inter,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,2,2,0,2)
      call set_rule(op_x_inter,ttype_op,DEF_R12INTERM,
     &              op_x_inter,1,1,
     &              parameters,1,tgt_info)

c      if(r12fix)then
c        call add_target(op_x1_inter,ttype_op,.false.,tgt_info)
c        occ_def = 0
c        occ_def(ihole,1,1) = 1
c        occ_def(ihole,2,2) = 1
c        call op_from_occ_parameters(-1,parameters,2,
c     &       occ_def,1,2,1)
c        call set_rule(op_x1_inter,ttype_op,DEF_OP_FROM_OCC,
c     &                op_x1_inter,1,1,
c     &                parameters,2,tgt_info)
c      endif


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
      ! 1
c      occ_def(IHOLE,1,1) = 1
c      occ_def(IHOLE,1,2) = 1
c      occ_def(IHOLE,2,2) = 2
      ! 2
c      occ_def(IHOLE,1,3) = 1
c      occ_def(IPART,1,4) = 1
c      occ_def(IHOLE,2,4) = 2
c      ! 3 - symmetrized: currently a problem
      occ_def(IPART,1,1) = 2
      occ_def(IHOLE,2,1) = 2
      ! 3 - not symmetrized:
c      occ_def(IPART,1,5) = 1
c      occ_def(IPART,1,6) = 1
c      occ_def(IHOLE,2,6) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,1,1,6)
      call set_rule(op_c_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_c_inter,1,1,
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
      call set_dependency(op_x_inv,op_x_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_x_inter,.false.) ! <- dagger=.false.
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
      

c dbg
      ! F12 integral contract with C12 for testing.
      if(r12fix)then
        call add_target(op_f12int_fix,ttype_op,.false.,tgt_info)
        call set_dependency(op_f12int_fix,op_rint,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_rint,.false.)
        call set_rule(op_f12int_fix,ttype_op,CLONE_OP,
     &                op_f12int_fix,1,1,
     &                parameters,1,tgt_info)

      ! [T12,F12] integral contract with C12 for testing.
        call add_target(op_t12f12int_fix,ttype_op,.false.,tgt_info)
        call set_dependency(op_t12f12int_fix,op_ttr,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_ttr,.false.)
        call set_rule(op_t12f12int_fix,ttype_op,CLONE_OP,
     &                op_t12f12int_fix,1,1,
     &                parameters,1,tgt_info)

      ! G12.F12 integral contract with C12 for testing.
        call add_target(op_g12f12int_fix,ttype_op,.false.,tgt_info)
        call set_dependency(op_g12f12int_fix,op_gr,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_gr,.false.)
        call set_rule(op_g12f12int_fix,ttype_op,CLONE_OP,
     &                op_g12f12int_fix,1,1,
     &                parameters,1,tgt_info)

      endif
c dbg

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
      call set_rule(form_r12_vint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

c      if(r12fix)then
c        ! formal definition of V0
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = form_r12_v0int
c        labels(2) = op_v0_inter
c        labels(3) = op_ham
c        labels(4) = op_r12
c        call add_target(form_r12_v0int,ttype_frm,.false.,tgt_info)
c        call set_dependency(form_r12_v0int,op_v0_inter,tgt_info)
c        call set_dependency(form_r12_v0int,op_ham,tgt_info)
c        call set_dependency(form_r12_v0int,op_r12,tgt_info)
c        call form_parameters(-1,
c     &       parameters,2,title_r12_v0int,0,'gr')
c        call set_rule(form_r12_v0int,ttype_frm,DEF_R12INTM_FORMAL,
c     &                labels,4,1,
c     &                parameters,2,tgt_info)
c      endif

      ! CABS approximation to V
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_vcabs
      labels(2) = op_v_inter
      labels(3) = op_g_x !op_ham
      labels(4) = op_rint
c      labels(5) = op_unity
      labels(5) = op_gr   
      ! F12: op_gr
      call add_target(form_r12_vcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_vcabs,op_v_inter,tgt_info)
c      call set_dependency(form_r12_vcabs,op_unity,tgt_info)
      call set_dependency(form_r12_vcabs,op_gr,tgt_info)
      call set_dependency(form_r12_vcabs,op_g_x,tgt_info)
      call set_dependency(form_r12_vcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vcabs,ansatz,'V '//approx)
      call set_rule(form_r12_vcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)

c      ! CABS approxiamtion to V0
c      if(r12fix)then
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = form_r12_v0cabs
c        labels(2) = op_v0_inter
cc        labels(3) = op_v_inter
cc dbg
c        labels(3) = op_g_x      !op_ham
c        labels(4) = op_rint
c        labels(5) = op_gr   
cc dbg
c        ! F12: op_gr
c        call add_target(form_r12_v0cabs,ttype_frm,.false.,tgt_info)
c        call set_dependency(form_r12_v0cabs,op_v0_inter,tgt_info)
c        call set_dependency(form_r12_vcabs,op_unity,tgt_info)
c        call set_dependency(form_r12_v0cabs,op_v_inter,tgt_info)
cc dbg
c        call set_dependency(form_r12_v0cabs,op_gr,tgt_info)
c        call set_dependency(form_r12_v0cabs,op_g_x,tgt_info)
c        call set_dependency(form_r12_v0cabs,op_rint,tgt_info)
cc dbg
c        call form_parameters(-1,
c     &       parameters,2,title_r12_v0cabs,ansatz,'V0'//approx)
c        call set_rule(form_r12_v0cabs,ttype_frm,DEF_R12INTM_CABS,
c     &                labels,5,1,
c     &                parameters,2,tgt_info)
c      endif

      ! formal definition of X
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xint
      labels(2) = op_x_inter
c      labels(3) = op_rba
      labels(3) = op_r12
      labels(4) = op_r12
      call add_target(form_r12_xint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xint,op_x_inter,tgt_info)
      call set_dependency(form_r12_xint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xint,0,'rxr')
      call set_rule(form_r12_xint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

c      if(r12fix)then
c        ! formal definition of X1
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = form_r12_x1int
c        labels(2) = op_x1_inter
cc        labels(3) = op_rba
c        labels(3) = op_r12
c        labels(4) = op_r12
c        call add_target(form_r12_x1int,ttype_frm,.false.,tgt_info)
c        call set_dependency(form_r12_x1int,op_x1_inter,tgt_info)
c        call set_dependency(form_r12_x1int,op_r12,tgt_info)
c        call form_parameters(-1,
c     &       parameters,2,title_r12_x1int,0,'rxr')
c        call set_rule(form_r12_x1int,ttype_frm,DEF_R12INTM_FORMAL,
c     &                labels,4,1,
c     &                parameters,2,tgt_info)
c      endif

      ! CABS approximation to X
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xcabs
      labels(2) = op_x_inter
      labels(3) = op_rint
      labels(4) = op_rint
      labels(5) = op_ff
      call add_target(form_r12_xcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xcabs,op_x_inter,tgt_info)
      call set_dependency(form_r12_xcabs,op_ff,tgt_info)
      call set_dependency(form_r12_xcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xcabs,ansatz,'X '//approx)
      call set_rule(form_r12_xcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)

c      if(r12fix)then
c        ! CABS approximation to X1
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = form_r12_x1cabs
c        labels(2) = op_x1_inter
cc        labels(3) = op_x_inter
cc dbg
c        labels(3) = op_rint
c        labels(4) = op_rint
c        labels(5) = op_ff
cc dbg
c        call add_target(form_r12_x1cabs,ttype_frm,.false.,tgt_info)
c        call set_dependency(form_r12_x1cabs,op_x1_inter,tgt_info)
cc        call set_dependency(form_r12_x1cabs,op_x_inter,tgt_info)
cc dbg
c        call set_dependency(form_r12_xcabs,op_ff,tgt_info)
c        call set_dependency(form_r12_xcabs,op_rint,tgt_info)
cc dbg
c        call form_parameters(-1,
c     &       parameters,2,title_r12_x1cabs,ansatz,'X1'//approx)
c        call set_rule(form_r12_x1cabs,ttype_frm,DEF_R12INTM_CABS,
c     &                labels,5,1,
c     &                parameters,2,tgt_info)
c      endif

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
      call set_rule(form_r12_bint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)

c      if(r12fix)then
c        ! formal definition of B0
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = form_r12_b0int
c        labels(2) = op_b0_inter
cc        labels(3) = op_rba
c        labels(3) = op_r12
c        labels(4) = op_ham
c        labels(5) = op_r12
c        call add_target(form_r12_b0int,ttype_frm,.false.,tgt_info)
c        call set_dependency(form_r12_b0int,op_b0_inter,tgt_info)
cc      call set_dependency(form_r12_bint,op_rba,tgt_info)
c        call set_dependency(form_r12_b0int,op_ham,tgt_info)
c        call set_dependency(form_r12_b0int,op_r12,tgt_info)
c        call form_parameters(-1,
c     &       parameters,2,title_r12_b0int,0,'rfr')
c        call set_rule(form_r12_b0int,ttype_frm,DEF_R12INTM_FORMAL,
c     &                labels,5,1,
c     &                parameters,2,tgt_info)
c      endif

      ! CABS approximation to B
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_r12_bcabs
      labels(2) = op_b_inter
      labels(3) = op_rint
      labels(4) = op_ttr
c      labels(5) = op_unity
      labels(5) = op_rttr
      nlab = 5
      call add_target(form_r12_bcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_bcabs,op_b_inter,tgt_info)
      call set_dependency(form_r12_bcabs,op_rttr,tgt_info)
c      call set_dependency(form_r12_bcabs,op_unity,tgt_info)
      call set_dependency(form_r12_bcabs,op_rint,tgt_info)
      if (approx(1:2).eq.'A''') then
        call set_dependency(form_r12_bcabs,op_x_inter,tgt_info)
        call set_dependency(form_r12_bcabs,form_r12_xcabs,tgt_info)
        call set_dependency(form_r12_bcabs,op_ham,tgt_info)
        labels(6) = op_x_inter
        labels(7) = op_ham
        nlab = 7        
      else if (approx(1:1).eq.'B') then
        call set_dependency(form_r12_bcabs,op_x_inter,tgt_info)
        call set_dependency(form_r12_bcabs,form_r12_xcabs,tgt_info)
        call set_dependency(form_r12_bcabs,op_rintbar,tgt_info)        
        call set_dependency(form_r12_bcabs,op_rinttilde,tgt_info)        
        call set_dependency(form_r12_bcabs,op_ffbar,tgt_info)        
        call set_dependency(form_r12_bcabs,op_rintbreve,tgt_info)        
        labels(6) = op_x_inter
        labels(7) = op_ham
        labels(8) = op_rintbar
        labels(9) = op_rinttilde
        labels(10) = op_ffbar
        labels(11) = '-'
        labels(12) = op_rintbreve
        nlab = 12
      else if (approx(1:1).eq.'C') then
        call set_dependency(form_r12_bcabs,op_rdagbar,tgt_info)        
        call set_dependency(form_r12_bcabs,op_rinttilde,tgt_info)        
        call set_dependency(form_r12_bcabs,op_ffbar,tgt_info)        
        call set_dependency(form_r12_bcabs,op_rintbreve,tgt_info)        
        labels(6) = '-'
        labels(7) = '-'
        labels(8) = op_rdagbar
        labels(9) = op_rinttilde
        labels(10) = op_ffbar
        labels(11) = '-'
        labels(12) = op_rintbreve
        nlab = 12
      end if
      if (ansatz.gt.1) then
        call set_dependency(form_r12_bcabs,op_c_inter,tgt_info)
        labels(13) = op_c_inter
        nlab = 13
      end if
      approx(12:12) = 'S' ! set symmetrization flag
c test
c      approx(8:10)  = 'HY1'
c
      call form_parameters(-1,
     &     parameters,2,title_r12_bcabs,ansatz,'B '//approx)
      approx(12:12) = ' ' ! unset flag
c test
      approx(8:10)  = '   '
c
      call set_rule(form_r12_bcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,nlab,1,
     &              parameters,2,tgt_info)

c      if(r12fix)then
c        ! CABS approximation to B0
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = form_r12_b0cabs
c        labels(2) = op_b0_inter
c        labels(3) = op_b_inter
c        nlab = 3
c        call add_target(form_r12_b0cabs,ttype_frm,.false.,tgt_info)
c        call set_dependency(form_r12_b0cabs,op_b0_inter,tgt_info)
c        call set_dependency(form_r12_b0cabs,op_b_inter,tgt_info)
c        call form_parameters(-1,
c     &       parameters,2,title_r12_b0cabs,ansatz,'B0'//approx)
c        call set_rule(form_r12_b0cabs,ttype_frm,DEF_R12INTM_CABS,
c     &                labels,nlab,1,
c     &                parameters,2,tgt_info)
c      endif

      ! formal definition of C intermediate
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
c     &     parameters,2,title_r12_cint,0,'fxr')
     &     parameters,2,title_r12_cint,0,'fr')
      call set_rule(form_r12_cint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to C intermediate
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_ccabs
      labels(2) = op_c_inter
      if (approx.eq.'B') then
        labels(3) = op_rint
        labels(4) = op_ham
        labels(5) = op_ttr
        labels(6) = op_rintbar
        labels(7) = op_rinttilde
        nlab = 7
      else ! use this one for 3A and 3C:
        labels(3) = op_rint
        labels(4) = op_ham
        nlab = 4
      end if
      call add_target(form_r12_ccabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_ccabs,op_c_inter,tgt_info)
      call set_dependency(form_r12_ccabs,op_rint,tgt_info)
      call set_dependency(form_r12_ccabs,op_ham,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_ccabs,ansatz,'C '//approx)
c     &     parameters,2,title_r12_ccabs,ansatz,'C '//
c     &     'C           ')
      call set_rule(form_r12_ccabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,nlab,1,
     &              parameters,2,tgt_info)

      ! formal definition of P
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_pint
      labels(2) = op_p_inter
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_r12_pint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_pint,op_p_inter,tgt_info)
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
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_r12_zint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_zint,op_z_inter,tgt_info)
      call set_dependency(form_r12_zint,op_r12,tgt_info)
      call set_dependency(form_r12_zint,op_ham,tgt_info)
      call set_dependency(form_r12_zint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_zint,0,'rxgxr')
      call set_rule(form_r12_zint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)

c dbg
      ! Test evaluation of integrals connected to fixed amplitudes.
      if(r12fix)then
        
        ! F12 integral connects once to t'.
        labels(1:)(1:len_target_name) = ' '
        labels(1) = form_f12int_fix
        labels(2) = op_f12int_fix
        labels(3) = op_c12
        labels(4) = op_rint
        call add_target(form_f12int_fix,ttype_frm,.false.,tgt_info)
        call set_dependency(form_f12int_fix,op_f12int_fix,tgt_info)
        call set_dependency(form_f12int_fix,op_rint,tgt_info)
        call set_dependency(form_f12int_fix,op_c12,tgt_info)
        call form_parameters(-1,
     &       parameters,2,title_f12int_fix,0,'F12.C12')
        call set_rule(form_f12int_fix,ttype_frm,CONTRACT,
     &       labels,4,1,
     &       parameters,2,tgt_info)

        ! [T12,F12] integral connects once to t'.
        labels(1:)(1:len_target_name) = ' '
        labels(1) = form_t12f12int_fix
        labels(2) = op_t12f12int_fix
        labels(3) = op_c12
        labels(4) = op_ttr
        call add_target(form_t12f12int_fix,ttype_frm,.false.,tgt_info)
        call set_dependency(form_t12f12int_fix,
     &       op_t12f12int_fix,tgt_info)
        call set_dependency(form_t12f12int_fix,op_ttr,tgt_info)
        call set_dependency(form_t12f12int_fix,op_c12,tgt_info)
        call form_parameters(-1,
     &       parameters,2,title_t12f12int_fix,0,'T12F12.C12')
        call set_rule(form_t12f12int_fix,ttype_frm,CONTRACT,
     &       labels,4,1,
     &       parameters,2,tgt_info)

        ! g12.F12 integral connects once to t'.
        labels(1:)(1:len_target_name) = ' '
        labels(1) = form_g12f12int_fix
        labels(2) = op_g12f12int_fix
        labels(3) = op_c12
        labels(4) = op_gr
        call add_target(form_g12f12int_fix,ttype_frm,.false.,tgt_info)
        call set_dependency(form_g12f12int_fix,
     &       op_g12f12int_fix,tgt_info)
        call set_dependency(form_g12f12int_fix,op_gr,tgt_info)
        call set_dependency(form_g12f12int_fix,op_c12,tgt_info)
        call form_parameters(-1,
     &       parameters,2,title_g12f12int_fix,0,'G12.F12.C12')
        call set_rule(form_g12f12int_fix,ttype_frm,CONTRACT,
     &       labels,4,1,
     &       parameters,2,tgt_info)

      endif
c dbg

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
      call set_dependency(fopt_r12_vcabs,mel_gr,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_v_def,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_gintx,tgt_info)
c      call set_dependency(fopt_r12_vcabs,mel_ham,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_vcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

c      if(r12fix)then
c        ! set V0
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = fopt_r12_v0cabs
c        labels(2) = form_r12_v0cabs
c        ncat = 1
c        nint = 0
c        call add_target(fopt_r12_v0cabs,ttype_frm,.false.,tgt_info)
c        call set_dependency(fopt_r12_v0cabs,form_r12_v0cabs,tgt_info)
c        call set_dependency(fopt_r12_v0cabs,mel_v0_def,tgt_info)
cc        call set_dependency(fopt_r12_v0cabs,mel_v_def,tgt_info)
cc dbg
c        call set_dependency(fopt_r12_v0cabs,mel_gr,tgt_info)
c        call set_dependency(fopt_r12_v0cabs,mel_gintx,tgt_info)
c        call set_dependency(fopt_r12_v0cabs,mel_rint,tgt_info)      
cc dbg
c        call opt_parameters(-1,parameters,ncat,nint)
c        call set_rule(fopt_r12_v0cabs,ttype_frm,OPTIMIZE,
c     &                labels,ncat+nint+1,1,
c     &                parameters,1,tgt_info)
c      endif

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
      call set_dependency(fopt_r12_xcabs,mel_ff,tgt_info)
c      call set_dependency(fopt_r12_xcabs,mel_rinba,tgt_info)
      call set_dependency(fopt_r12_xcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_xcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

c      if(r12fix)then
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = fopt_r12_x1cabs
c        labels(2) = form_r12_x1cabs
c        ncat = 1
c        nint = 0
c        call add_target(fopt_r12_x1cabs,ttype_frm,.false.,tgt_info)
c        call set_dependency(fopt_r12_x1cabs,form_r12_x1cabs,tgt_info)
c        call set_dependency(fopt_r12_x1cabs,mel_x1_def,tgt_info)
cc        call set_dependency(fopt_r12_x1cabs,mel_x_def,tgt_info)
cc dbg
c        call set_dependency(fopt_r12_x1cabs,mel_ff,tgt_info)
c        call set_dependency(fopt_r12_x1cabs,mel_rint,tgt_info)      
cc dbg
c        call opt_parameters(-1,parameters,ncat,nint)
c        call set_rule(fopt_r12_x1cabs,ttype_frm,OPTIMIZE,
c     &                labels,ncat+nint+1,1,
c     &                parameters,1,tgt_info)
c      endif

      ! set B
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_bcabs
      labels(2) = form_r12_bcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_bcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_bcabs,form_r12_bcabs,tgt_info)
      call set_dependency(fopt_r12_bcabs,mel_b_def,tgt_info)
      call set_dependency(fopt_r12_bcabs,mel_rttr,tgt_info)
      call set_dependency(fopt_r12_bcabs,mel_ham,tgt_info)
c      call set_dependency(fopt_r12_bcabs,mel_rinba,tgt_info)      
      call set_dependency(fopt_r12_bcabs,mel_rint,tgt_info)      
      if (approx(1:2).ne.'A '.and.approx(1:1).ne.'C') then
        call set_dependency(fopt_r12_bcabs,mel_x_def,tgt_info)
        if (approx(1:1).eq.'B') then
          call set_dependency(form_r12_bcabs,mel_rinttilde,tgt_info)        
          call set_dependency(form_r12_bcabs,mel_rintbar,tgt_info)        
          call set_dependency(form_r12_bcabs,mel_ffbar,tgt_info)        
          call set_dependency(form_r12_bcabs,mel_rdagbar,tgt_info)        
          call set_dependency(form_r12_bcabs,mel_rintbreve,tgt_info)        
        end if
      else if (approx(1:1).eq.'C') then
        call set_dependency(form_r12_bcabs,mel_rdagbar,tgt_info)        
        call set_dependency(form_r12_bcabs,mel_rinttilde,tgt_info)        
        call set_dependency(form_r12_bcabs,mel_ffbar,tgt_info)        
        call set_dependency(form_r12_bcabs,mel_rdagbar,tgt_info)        
        call set_dependency(form_r12_bcabs,mel_rintbreve,tgt_info)        
      end if
      if (ansatz.gt.1) then
        call set_dependency(form_r12_bcabs,mel_c_def,tgt_info)        
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_bcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

c      if(r12fix)then
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = fopt_r12_b0cabs
c        labels(2) = form_r12_b0cabs
c        ncat = 1
c        nint = 0
c        call add_target(fopt_r12_b0cabs,ttype_frm,.false.,tgt_info)
c        call set_dependency(fopt_r12_b0cabs,form_r12_b0cabs,tgt_info)
c        call set_dependency(fopt_r12_b0cabs,mel_b0_def,tgt_info)
c        call set_dependency(fopt_r12_b0cabs,mel_b_def,tgt_info)
c        call opt_parameters(-1,parameters,ncat,nint)
c        call set_rule(fopt_r12_b0cabs,ttype_frm,OPTIMIZE,
c     &                labels,ncat+nint+1,1,
c     &                parameters,1,tgt_info)
c      endif

      ! set C intermediate
      ! currently approx C only
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_ccabs
      labels(2) = form_r12_ccabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_ccabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_ccabs,form_r12_ccabs,tgt_info)
      call set_dependency(fopt_r12_ccabs,mel_c_def,tgt_info)
c      call set_dependency(fopt_r12_ccabs,mel_rintc,tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_ccabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)


c dbg
      ! Fixed amplitude integrals for testing.
      if(r12fix)then
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = fopt_f12int_fix
        labels(2) = form_f12int_fix
        ncat = 1
        nint = 0
        call add_target(fopt_f12int_fix,ttype_frm,.false.,tgt_info)
        call set_dependency(fopt_f12int_fix,form_f12int_fix,tgt_info)
        call set_dependency(fopt_f12int_fix,mel_f12fix_def,tgt_info)
        call opt_parameters(-1,parameters,ncat,nint)
        call set_rule(fopt_f12int_fix,ttype_frm,OPTIMIZE,
     &                labels,ncat+nint+1,1,
     &                parameters,1,tgt_info)

      endif
c dbg

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
c      call add_target(mel_rinba,ttype_opme,.false.,tgt_info)
c      call set_dependency(mel_rinba,op_rinba,tgt_info)
c      ! (a) define
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_rinba
c      labels(2) = op_rinba
c      call me_list_parameters(-1,parameters,
c     &     0,0,1,0,0)
c      call set_rule(mel_rinba,ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)
c      ! (b) import
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_rinba
c      call import_parameters(-1,parameters,env_type)
c      call set_rule(mel_rinba,ttype_opme,IMPORT,
c     &              labels,1,1,
c     &              parameters,1,tgt_info)

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
      call add_target(mel_ff,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ff,op_ff,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ff
      labels(2) = op_ff
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_ff,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ff
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_ff,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! {R12^2}BAR integrals
      call add_target(mel_ffbar,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ffbar,op_ffbar,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ffbar
      labels(2) = op_ffbar
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_ffbar,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ffbar
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_ffbar,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! G.R12 integrals (for f(R12))
      call add_target(mel_gr,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_gr,op_gr,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_gr
      labels(2) = op_gr
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_gr,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_gr
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_gr,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

c      ! R12C integrals (for approx C)
c      call add_target(mel_rintc,ttype_opme,.false.,tgt_info)
c      call set_dependency(mel_rintc,op_rintc,tgt_info)
c      ! (a) define
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_rintc
c      labels(2) = op_rintc
c      call me_list_parameters(-1,parameters,
c     &     0,0,1,0,0)
c      call set_rule(mel_rintc,ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)
c      ! (b) import
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_rintc
c      call import_parameters(-1,parameters,env_type)
c      call set_rule(mel_rintc,ttype_opme,IMPORT,
c     &              labels,1,1,
c     &              parameters,1,tgt_info)

      ! R12BAR integrals
      call add_target(mel_rintbar,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_rintbar,op_rintbar,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rintbar
      labels(2) = op_rintbar
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_rintbar,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rintbar
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_rintbar,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! R12tilde integrals
      call add_target(mel_rinttilde,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_rinttilde,op_rinttilde,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rinttilde
      labels(2) = op_rinttilde
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_rinttilde,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rinttilde
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_rinttilde,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! R12BAR^+ integrals
      call add_target(mel_rdagbar,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_rdagbar,op_rdagbar,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rdagbar
      labels(2) = op_rdagbar
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_rdagbar,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rdagbar
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_rdagbar,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! R12BREVE integrals
      call add_target(mel_rintbreve,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_rintbreve,op_rintbreve,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rintbreve
      labels(2) = op_rintbreve
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_rintbreve,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rintbreve
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_rintbreve,ttype_opme,IMPORT,
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
c      ! RBAR-list
c      call add_target(meldef_rintbar,ttype_opme,.false.,tgt_info)
c      call set_dependency(meldef_rintbar,op_rintbar,tgt_info)
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_rintbar
c      labels(2) = op_rintbar
c      call me_list_parameters(-1,parameters,
c     &     0,0,1,0,0)
c      call set_rule(meldef_rintbar,ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)
c
c      ! RTILDE-list
c      call add_target(meldef_rinttilde,ttype_opme,.false.,tgt_info)
c      call set_dependency(meldef_rinttilde,op_rinttilde,tgt_info)
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_rinttilde
c      labels(2) = op_rinttilde
c      call me_list_parameters(-1,parameters,
c     &     0,0,1,0,0)
c      call set_rule(meldef_rinttilde,ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)

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

c      if(r12fix)then
c        ! V0-list
c        call add_target(mel_v0_def,ttype_opme,.false.,tgt_info)
c        call set_dependency(mel_v0_def,op_v0_inter,tgt_info)
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = mel_v0_inter
c        labels(2) = op_v0_inter
c        call me_list_parameters(-1,parameters,
c     &       0,0,1,0,0)
c        call set_rule(mel_v0_def,ttype_opme,DEF_ME_LIST,
c     &                labels,2,1,
c     &                parameters,1,tgt_info)
c      endif

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

c      if(r12fix)then
c        ! X1-list
c        call add_target(mel_x1_def,ttype_opme,.false.,tgt_info)
c        call set_dependency(mel_x1_def,op_x1_inter,tgt_info)
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = mel_x1_inter
c        labels(2) = op_x1_inter
c        call me_list_parameters(-1,parameters,
c     &       0,0,1,0,0)
c        call set_rule(mel_x1_def,ttype_opme,DEF_ME_LIST,
c     &                labels,2,1,
c     &                parameters,1,tgt_info)
c      endif

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

c      if(r12fix)then
c      ! B0-list
c        call add_target(mel_b0_def,ttype_opme,.false.,tgt_info)
c        call set_dependency(mel_b0_def,op_b0_inter,tgt_info)
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = mel_b0_inter
c        labels(2) = op_b0_inter
c        call me_list_parameters(-1,parameters,
c     &       0,0,1,0,0)
c        call set_rule(mel_b0_def,ttype_opme,DEF_ME_LIST,
c     &                labels,2,1,
c     &                parameters,1,tgt_info)
c      endif

      ! C-list
      call add_target(mel_c_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_c_def,op_c_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_c_inter
      labels(2) = op_c_inter
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_c_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

c dbg
      ! Test contractions of integrals with fixed amplitudes.
      if(r12fix)then
        call add_target(mel_f12fix_def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_f12fix_def,op_f12int_fix,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_f12fix
        labels(2) = op_f12int_fix
        call me_list_parameters(-1,parameters,
     &       0,0,1,0,0)
        call set_rule(mel_f12fix_def,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)

      endif
c dbg

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

      if(.not.r12fix)then
        ! diagonal of B(ij) for testing
        call add_target(mel_b_dia,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_b_dia,op_diar12,tgt_info)
        call set_dependency(mel_b_dia,eval_r12_inter,tgt_info)
        call set_dependency(mel_b_dia,mel_ham,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_b_dia
        labels(2) = op_diar12
        call me_list_parameters(-1,parameters,
     &       0,0,1,0,0)
        call set_rule(mel_b_dia,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        labels(1) = mel_b_dia   ! output
        labels(2) = mel_ham     ! input
        labels(3) = mel_b_inter ! input
        labels(4) = mel_x_inter ! input
        call set_rule(mel_b_dia,ttype_opme,PRECONDITIONER,
     &                labels,4,1,
     &                parameters,1,tgt_info)
      
        ! X^-1 for testing
        call add_target(mel_x_inv,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_x_inv,op_diar12,tgt_info)
        call set_dependency(mel_x_inv,eval_r12_inter,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_x_inv
        labels(2) = op_x_inter
        call me_list_parameters(-1,parameters,
     &       0,0,1,0,0)
        call set_rule(mel_x_inv,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        labels(1) = mel_x_inv   ! output
        labels(2) = mel_x_inter ! input
        call set_rule(mel_x_inv,ttype_opme,INVERT,
     &                labels,2,1,
     &                parameters,1,tgt_info)
      
      endif

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
      call set_dependency(eval_r12_inter,mel_ff,tgt_info)
      call set_dependency(eval_r12_inter,mel_v_def,tgt_info)
c      call set_dependency(eval_r12_inter,mel_vbar_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_x_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_b_def,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_vcabs,tgt_info)
c      call set_dependency(eval_r12_inter,fopt_r12_vbcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_xcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_bcabs,tgt_info)
c      if(r12fix)then
c        call set_dependency(eval_r12_inter,fopt_r12_v0cabs,tgt_info)
c        call set_dependency(eval_r12_inter,fopt_r12_x1cabs,tgt_info)
c        call set_dependency(eval_r12_inter,fopt_r12_b0cabs,tgt_info)
c      endif
      if (ansatz.ne.1)
     &     call set_dependency(eval_r12_inter,fopt_r12_ccabs,tgt_info)
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
      if (ansatz.ne.1) then
        labels(1) = fopt_r12_ccabs
        call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      end if
      labels(1) = fopt_r12_bcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

c      if(r12fix)then
c        labels(1) = fopt_r12_v0cabs
c        call set_rule(eval_r12_inter,ttype_opme,EVAL,
c     &       labels,1,0,
c     &       parameters,0,tgt_info)
c        labels(1) = fopt_r12_x1cabs
c        call set_rule(eval_r12_inter,ttype_opme,EVAL,
c     &       labels,1,0,
c     &       parameters,0,tgt_info)
c        labels(1) = fopt_r12_b0cabs
c        call set_rule(eval_r12_inter,ttype_opme,EVAL,
c     &       labels,1,0,
c     &       parameters,0,tgt_info)
c      endif

c dbg
      if(r12fix)then
        call add_target('EVALINTS',ttype_gen,.false.,tgt_info)
        call set_dependency('EVALINTS',mel_rint,tgt_info)
        call set_dependency('EVALINTS',mel_c12def,tgt_info)
        call set_dependency('EVALINTS',mel_f12fix,tgt_info)
        call set_dependency('EVALINTS',fopt_f12int_fix,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = fopt_f12int_fix
        call set_rule('EVALINTS',ttype_opme,EVAL,
     &       labels,1,0,
     &       parameters,0,tgt_info)
      endif
c dbg

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
