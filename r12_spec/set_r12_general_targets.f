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
     &     occ_def(ngastp,2,20), ndef, mode
      logical ::
     &     needed, r12fix, extend, truncate
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
      call get_argument_value('method.R12','extend',ival=mode)
      call get_argument_value('method.R12','truncate',lval=truncate)

      ! actual processing moved to set_r12f_general_targets
      extend = mode.gt.0

      if(extend.and..not.r12fix)
     &     call quit(1,'set_r12_general_targets',
     &     'Extension only valid for MP2-R12 with fixed C12')

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
      if(.not.extend)then
c        min_rank = 2  ! 1 is a possibility 
        call r12gem_parameters(-1,parameters,
     &                         0,2,2,ansatz)
        call set_rule(op_r12,ttype_op,DEF_R12GEMINAL,
     &                op_r12,1,1,
     &                parameters,1,tgt_info)
      else
        occ_def = 0
        ! 1
        occ_def(IEXTR,1,1) = 2
        occ_def(IHOLE,2,1) = 2
        ! 2
        occ_def(IEXTR,1,2) = 2
        occ_def(IPART,2,2) = 1
        occ_def(IHOLE,2,2) = 1
        ndef = 2
        if(ansatz.gt.1)then
          ! 3
          occ_def(IEXTR,1,3) = 1
          occ_def(IPART,1,3) = 1
          occ_def(IHOLE,2,3) = 2
          ! 4
          occ_def(IEXTR,1,4) = 1
          occ_def(IPART,1,4) = 1
          occ_def(IPART,2,4) = 1
          occ_def(IHOLE,2,4) = 1
          ndef = 4
        endif

        call op_from_occ_parameters(-1,parameters,2,
     &                              occ_def,ndef,1,ndef)
        call set_rule(op_r12,ttype_op,DEF_OP_FROM_OCC,
     &                op_r12,1,1,
     &                parameters,2,tgt_info)
      endif

      ! Only need coefficients if optimising the R12 contribution.
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
     &     op_c12,.true.)       ! <- dagger=.true.
      call set_rule(op_cba,ttype_op,CLONE_OP,
     &              op_cba,1,1,
     &              parameters,1,tgt_info)

      if(extend)then
        ! T1' operators for extended MP2-F12.
        call add_target(op_cex,ttype_op,.false.,tgt_info)
        call xop_parameters(-1,parameters,
     &       .false.,1,1,0,2)
        call set_rule(op_cex,ttype_op,DEF_EXCITATION,
     &                op_cex,1,1,
     &                parameters,1,tgt_info)

        ! The Lagrangian multipliers.
        call add_target(op_cexbar,ttype_op,.false.,tgt_info)
        call set_dependency(op_cexbar,op_cex,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_cex,.true.) ! <- dagger=.true.
        call set_rule(op_cexbar,ttype_op,CLONE_OP,
     &                op_cexbar,1,1,
     &                parameters,1,tgt_info)
      endif


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
     &     0,min_rank,2,0,4) ! 4: two externals for 2-el ops
      call set_rule(op_rint,ttype_op,DEF_R12INT,
     &              op_rint,1,1,
     &              parameters,1,tgt_info)
      
      ! (pq)_frozen/(pq)_ae block of 2e-Hamiltonian
      call add_target(op_g_x,ttype_op,.false.,tgt_info)
      min_rank = 2 
c      call hop_parameters(-1,parameters,
c     &     min_rank,2,3,.true.)        ! 1-external at most
c      call set_rule(op_g_x,ttype_op,DEF_HAMILTONIAN,
c     &     op_g_x,1,1,
c     &     parameters,1,tgt_info)
      call r12int_parameters(-1,parameters,
     &     2,min_rank,2,0,3)
c dbg
c     &     2,min_rank,2,0,5)
c dbg
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
     &     0,min_rank,2,0,2)
      call set_rule(op_ttr,ttype_op,DEF_R12INT,
     &              op_ttr,1,1,
     &              parameters,1,tgt_info)
      
      ! A(f+k) modified integrals r12bar
      call add_target(op_rintbar,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
     &     0,min_rank,2,0,2)
      call set_rule(op_rintbar,ttype_op,DEF_R12INT,
     &              op_rintbar,1,1,
     &              parameters,1,tgt_info)
      
      ! C(f+k) modified integrals r12bar+
      call add_target(op_rdagbar,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
     &     0,min_rank,2,0,2)
      call set_rule(op_rdagbar,ttype_op,DEF_R12INT,
     &              op_rdagbar,1,1,
     &              parameters,1,tgt_info)
      
      ! C(f) modified integrals r12breve
      call add_target(op_rintbreve,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
     &     0,min_rank,2,0,2)
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
      
      ! commutator integrals <kl|r12[T1+T2,r12]|ij>
      call add_target(op_rttr,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,2,2,0,2)
      call set_rule(op_rttr,ttype_op,DEF_R12INTERM,
     &              op_rttr,1,1,
     &              parameters,1,tgt_info)

      ! (G.R)^{ij}_{pq}
      call add_target(op_gr,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,2,2) = 2
      ! 2
      occ_def(IHOLE,1,3) = 1
      occ_def(IPART,1,3) = 1
      occ_def(IHOLE,2,4) = 2
      ! 3
      occ_def(IPART,1,5) = 2
      occ_def(IHOLE,2,6) = 2
      ! 4
      occ_def(IHOLE,1,7) = 1
      occ_def(IEXTR,1,7) = 1
      occ_def(IHOLE,2,8) = 2
      ! 5
      occ_def(IPART,1,9) = 1
      occ_def(IEXTR,1,9) = 1
      occ_def(IHOLE,2,10) = 2
      ! 6
      occ_def(IHOLE,1,11) = 2
      occ_def(IHOLE,2,12) = 1
      occ_def(IPART,2,12) = 1
      ! 7
      occ_def(IHOLE,1,13) = 2
      occ_def(IPART,2,14) = 2

      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,7,2,14)
      call set_rule(op_gr,ttype_op,DEF_OP_FROM_OCC,
     &              op_gr,1,1,
     &              parameters,2,tgt_info)

c      call set_dependency(op_gr,op_v_inter,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_v_inter,.false.) ! <- dagger=.false.
c      call set_rule(op_gr,ttype_op,CLONE_OP,
c     &              op_gr,1,1,
c     &              parameters,1,tgt_info)

      ! V^{ij}_{pq}
      call add_target(op_v_inter,ttype_op,.false.,tgt_info)
c      call xop_parameters(-1,parameters,
c     &     .false.,2,2,0,2)
c      call set_rule(op_v_inter,ttype_op,DEF_R12INTERM,
c     &              op_v_inter,1,1,
c     &              parameters,1,tgt_info)
      occ_def = 0
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,2,2) = 2
      ! 2
      occ_def(IHOLE,1,3) = 1
      occ_def(IPART,1,3) = 1
      occ_def(IHOLE,2,4) = 2
      ! 3
      occ_def(IPART,1,5) = 2
      occ_def(IHOLE,2,6) = 2
      ! 4
      occ_def(IHOLE,1,7) = 1
      occ_def(IEXTR,1,7) = 1
      occ_def(IHOLE,2,8) = 2
      ! 5
      occ_def(IPART,1,9) = 1
      occ_def(IEXTR,1,9) = 1
      occ_def(IHOLE,2,10) = 2

      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,5,2,10)
      call set_rule(op_v_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_v_inter,1,1,
     &              parameters,2,tgt_info)
      
      ! V3
      call add_target(op_v3_inter,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,2,1) = 1
      occ_def(IPART,1,2) = 1
      occ_def(IHOLE,2,2) = 2
      ! 2
      occ_def(IHOLE,1,3) = 2
      occ_def(IPART,2,3) = 1
      occ_def(IPART,1,4) = 1
      occ_def(IHOLE,2,4) = 2
      ! 3
      occ_def(IHOLE,1,5) = 1
      occ_def(IPART,1,5) = 1
      occ_def(IHOLE,2,5) = 1
      occ_def(IPART,1,6) = 1
      occ_def(IHOLE,2,6) = 2
      ! 4
      occ_def(IHOLE,1,7) = 1
      occ_def(IPART,1,7) = 1
      occ_def(IPART,2,7) = 1
      occ_def(IPART,1,8) = 1
      occ_def(IHOLE,2,8) = 2
      ! 5
      occ_def(IPART,1,9) = 2
      occ_def(IHOLE,2,9) = 1
      occ_def(IPART,1,10) = 1
      occ_def(IHOLE,2,10) = 2
      ! 6
      occ_def(IPART,1,11) = 2
      occ_def(IPART,2,11) = 1
      occ_def(IPART,1,12) = 1
      occ_def(IHOLE,2,12) = 2

      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,6,2,12)
      call set_rule(op_v3_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_v3_inter,1,1,
     &              parameters,2,tgt_info)

      ! V4
      call add_target(op_v4_inter,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IPART,1,2) = 1
      occ_def(IHOLE,2,2) = 2
      occ_def(IPART,1,3) = 1
      occ_def(IHOLE,2,3) = 2
      ! 2
      occ_def(IHOLE,1,4) = 1
      occ_def(IPART,1,4) = 1
      occ_def(IPART,1,5) = 1
      occ_def(IHOLE,2,5) = 2
      occ_def(IPART,1,6) = 1
      occ_def(IHOLE,2,6) = 2
      ! 3
      occ_def(IPART,1,7) = 2
      occ_def(IPART,1,8) = 1
      occ_def(IHOLE,2,8) = 2
      occ_def(IPART,1,9) = 1
      occ_def(IHOLE,2,9) = 2

      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,3,3,9)
      call set_rule(op_v4_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_v4_inter,1,1,
     &              parameters,2,tgt_info)

      ! B intermediate
      call add_target(op_b_inter,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,2,2,0,2)
      call set_rule(op_b_inter,ttype_op,DEF_R12INTERM,
     &              op_b_inter,1,1,
     &              parameters,1,tgt_info)

      ! R12^{2} integrals
      call add_target(op_ff,ttype_op,.false.,tgt_info)
c      if (approx(1:1).eq.'A') then
      if(is_keyword_set('method.CC').gt.0.and..not.truncate)then
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
        call set_rule(op_ff,ttype_op,DEF_OP_FROM_OCC,
     &                op_ff,1,1,
     &                parameters,2,tgt_info)
      else
        call set_dependency(op_ff,op_b_inter,tgt_info)
        call cloneop_parameters(-1,parameters,
     &       op_b_inter,.false.) ! <- dagger=.false.
        call set_rule(op_ff,ttype_op,CLONE_OP,
     &                op_ff,1,1,
     &                parameters,1,tgt_info)
      end if

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

c dbg
      call add_target(op_x_test,ttype_op,.false.,tgt_info)
      call set_dependency(op_x_test,op_x_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_x_inter,.false.) ! <- dagger=.false.
      call set_rule(op_x_test,ttype_op,CLONE_OP,
     &              op_x_test,1,1,
     &              parameters,1,tgt_info)
c dbg

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
      occ_def = 0
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,2,2) = 2

      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,1,2,2)
      call set_rule(op_p_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_p_inter,1,1,
     &              parameters,2,tgt_info)

c dbg
      call add_target(op_p_test,ttype_op,.false.,tgt_info)
      call set_dependency(op_p_test,op_p_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_p_inter,.false.) ! <- dagger=.false.
      call set_rule(op_p_test,ttype_op,CLONE_OP,
     &              op_p_test,1,1,
     &              parameters,1,tgt_info)
c dbg

      ! R12^{2}*G12 integrals
      call add_target(op_ffg,ttype_op,.false.,tgt_info)
      call set_dependency(op_ffg,op_p_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     op_p_inter,.false.)  ! <- dagger=.false.
      call set_rule(op_ffg,ttype_op,CLONE_OP,
     &              op_ffg,1,1,
     &              parameters,1,tgt_info)

      ! P3 intermediate:
      ! G-part
      call add_target(op_p3g_inter,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IPART,2,1) = 1
      occ_def(IHOLE,1,2) = 1
      occ_def(IHOLE,2,3) = 2

      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,1,3,3)
      call set_rule(op_p3g_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_p3g_inter,1,1,
     &              parameters,2,tgt_info)

      ! F-part
      call add_target(op_p3f_inter,ttype_op,.false.,tgt_info)
      call set_dependency(op_p3f_inter,op_p3g_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_p3g_inter,.false.) ! <- dagger=.false.
      call set_rule(op_p3f_inter,ttype_op,CLONE_OP,
     &              op_p3f_inter,1,1,
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
c      ! 3
c      occ_def(IHOLE,1,7) = 2
c      occ_def(IPART,1,8) = 1
c      occ_def(IHOLE,2,8) = 1
c      occ_def(IHOLE,2,9) = 2
c      ! 4
c      occ_def(IHOLE,1,10) = 2
c      occ_def(IPART,1,11) = 1
c      occ_def(IPART,2,11) = 1
c      occ_def(IHOLE,2,12) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,2,3,6)
      call set_rule(op_z_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_z_inter,1,1,
     &              parameters,2,tgt_info)

c dbg
      call add_target(op_z_test,ttype_op,.false.,tgt_info)
      call set_dependency(op_z_test,op_z_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_z_inter,.false.) ! <- dagger=.false.
      call set_rule(op_z_test,ttype_op,CLONE_OP,
     &              op_z_test,1,1,
     &              parameters,1,tgt_info)
c dbg

      ! Z4 intermediate
      call add_target(op_z4_inter,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,1,2) = 2
      occ_def(IHOLE,2,3) = 2
      occ_def(IHOLE,2,4) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,1,4,4)
      call set_rule(op_z4_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_z4_inter,1,1,
     &              parameters,2,tgt_info)

      ! K4 intermediate
      call add_target(op_k4_inter,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IPART,2,1) = 1
      occ_def(IHOLE,1,2) = 2
      occ_def(IPART,2,2) = 1
      occ_def(IHOLE,2,3) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,1,3,3)
      call set_rule(op_k4_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_k4_inter,1,1,
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

      ! formal definition of V3
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_v3int
      labels(2) = op_v3_inter
      labels(3) = op_r12
      labels(4) = op_ham
      call add_target(form_r12_v3int,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_v3int,op_v3_inter,tgt_info)
      call set_dependency(form_r12_v3int,op_ham,tgt_info)
      call set_dependency(form_r12_v3int,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vint,0,'V3')
      call set_rule(form_r12_v3int,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! formal definition of V4
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_v4int
      labels(2) = op_v4_inter
      labels(3) = op_r12
      labels(4) = op_ham
      call add_target(form_r12_v4int,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_v4int,op_v4_inter,tgt_info)
      call set_dependency(form_r12_v4int,op_ham,tgt_info)
      call set_dependency(form_r12_v4int,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_v4int,0,'V4')
      call set_rule(form_r12_v4int,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

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

c dbg
      ! formal definition of X-test
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_x_test
      labels(2) = op_x_test
      labels(3) = op_r12
      labels(4) = op_r12
      call add_target(form_x_test,ttype_frm,.false.,tgt_info)
      call set_dependency(form_x_test,op_x_test,tgt_info)
      call set_dependency(form_x_test,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_x_test,0,'rxr')
      call set_rule(form_x_test,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_x_test
      labels(2) = form_x_test
      labels(3) = op_r12
      labels(4) = op_rint
      call set_dependency(form_x_test,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_x_test,1,'---')
      call set_rule(form_x_test,ttype_frm,REPLACE,
     &     labels,4,1,
     &     parameters,2,tgt_info)
c dbg

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

      ! CABS approximation to B
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_r12_bcabs
      labels(2) = op_b_inter
      labels(3) = op_rint
      labels(4) = op_ttr
      labels(5) = op_rttr
      nlab = 5
      call add_target(form_r12_bcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_bcabs,op_b_inter,tgt_info)
      call set_dependency(form_r12_bcabs,op_rttr,tgt_info)
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
      call form_parameters(-1,
     &     parameters,2,title_r12_bcabs,ansatz,'B '//approx)
      approx(12:12) = ' ' ! unset flag
      call set_rule(form_r12_bcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,nlab,1,
     &              parameters,2,tgt_info)

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
     &     parameters,2,title_r12_pint,0,'P')
      call set_rule(form_r12_pint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to P
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_pcabs
      labels(2) = op_p_inter
      labels(3) = op_rint
      labels(4) = op_g_x
      labels(5) = op_gr
      labels(6) = op_v_inter
      labels(7) = op_ffg
      call add_target(form_r12_pcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_pcabs,op_p_inter,tgt_info)
      call set_dependency(form_r12_pcabs,op_ffg,tgt_info)
      call set_dependency(form_r12_pcabs,op_gr,tgt_info)
      call set_dependency(form_r12_pcabs,op_g_x,tgt_info)
      call set_dependency(form_r12_pcabs,op_rint,tgt_info)
      call set_dependency(form_r12_pcabs,op_v_inter,tgt_info)
      approx(12:12) = 'S' ! set symmetrization flag
      call form_parameters(-1,
     &     parameters,2,title_r12_pcabs,ansatz,'P '//approx)
      approx(12:12) = ' ' ! unset flag
      call set_rule(form_r12_pcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,7,1,
     &              parameters,2,tgt_info)

c dbg
      ! formal definition of P-test
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_p_test
      labels(2) = op_p_test
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_p_test,ttype_frm,.false.,tgt_info)
      call set_dependency(form_p_test,op_p_test,tgt_info)
      call set_dependency(form_p_test,op_ham,tgt_info)
      call set_dependency(form_p_test,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_p_test,0,'P')
      call set_rule(form_p_test,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_p_test
      labels(2) = form_p_test
      labels(3) = op_r12
      labels(4) = op_rint
      call set_dependency(form_p_test,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_p_test,1,'---')
      call set_rule(form_p_test,ttype_frm,REPLACE,
     &     labels,4,1,
     &     parameters,2,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_p_test
      labels(2) = form_p_test
      labels(3) = op_ham
      labels(4) = op_g_x
      call set_dependency(form_p_test,op_g_x,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_p_test,1,'---')
      call set_rule(form_p_test,ttype_frm,REPLACE,
     &     labels,4,1,
     &     parameters,2,tgt_info)
c dbg

      ! Formal definition of P3F
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_p3fint
      labels(2) = op_p3f_inter
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_r12_p3fint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_p3fint,op_p3f_inter,tgt_info)
      call set_dependency(form_r12_p3fint,op_ham,tgt_info)
      call set_dependency(form_r12_p3fint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_p3fint,0,'P3F')
      call set_rule(form_r12_p3fint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to P3F
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_p3fcabs
      labels(2) = op_p3f_inter
      labels(3) = op_rint
      labels(4) = op_ham
      call add_target(form_r12_p3fcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_p3fcabs,op_p3f_inter,tgt_info)
      call set_dependency(form_r12_p3fcabs,op_rint,tgt_info)
      call set_dependency(form_r12_p3fcabs,op_ham,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_p3fcabs,ansatz,'PF '//approx)
      call set_rule(form_r12_p3fcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! Formal definition of P3G
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_p3gint
      labels(2) = op_p3g_inter
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_r12_p3gint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_p3gint,op_p3g_inter,tgt_info)
      call set_dependency(form_r12_p3gint,op_ham,tgt_info)
      call set_dependency(form_r12_p3gint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_p3gint,0,'P3G')
      call set_rule(form_r12_p3gint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to P3G.
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_p3gcabs
      labels(2) = op_p3g_inter
      labels(3) = op_rint
      labels(4) = op_v_inter
      call add_target(form_r12_p3gcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_p3gcabs,op_p3g_inter,tgt_info)
      call set_dependency(form_r12_p3gcabs,op_rint,tgt_info)
      call set_dependency(form_r12_p3gcabs,op_v_inter,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_p3gcabs,ansatz,'PG '//approx)
      call set_rule(form_r12_p3gcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,4,1,
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
     &     parameters,2,title_r12_zint,0,'Z')
      call set_rule(form_r12_zint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to Z.
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_zcabs
      labels(2) = op_z_inter
      labels(3) = op_rint
      labels(4) = op_g_x
      labels(5) = op_ff
      call add_target(form_r12_zcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_zcabs,op_z_inter,tgt_info)
      call set_dependency(form_r12_zcabs,op_ff,tgt_info)
      call set_dependency(form_r12_zcabs,op_rint,tgt_info)
      call set_dependency(form_r12_zcabs,op_g_x,tgt_info)
c      approx(12:12) = 'S' ! set symmetrization flag
      call form_parameters(-1,
     &     parameters,2,title_r12_zcabs,ansatz,'Z '//approx)
c      approx(12:12) = ' ' ! unset flag
      call set_rule(form_r12_zcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)

c dbg
      ! formal definition of Z-test
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_z_test
      labels(2) = op_z_test
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_z_test,ttype_frm,.false.,tgt_info)
      call set_dependency(form_z_test,op_z_test,tgt_info)
      call set_dependency(form_z_test,op_ham,tgt_info)
      call set_dependency(form_z_test,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_z_test,0,'Z')
      call set_rule(form_z_test,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_z_test
      labels(2) = form_z_test
      labels(3) = op_r12
      labels(4) = op_rint
      call set_dependency(form_z_test,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_z_test,1,'---')
      call set_rule(form_z_test,ttype_frm,REPLACE,
     &     labels,4,1,
     &     parameters,2,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_z_test
      labels(2) = form_z_test
      labels(3) = op_ham
      labels(4) = op_g_x
      call set_dependency(form_z_test,op_g_x,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_z_test,1,'---')
      call set_rule(form_z_test,ttype_frm,REPLACE,
     &     labels,4,1,
     &     parameters,2,tgt_info)
c dbg

      ! Formal definition of Z4
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_z4int
      labels(2) = op_z4_inter
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      labels(6) = op_r12
      call add_target(form_r12_z4int,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_z4int,op_z4_inter,tgt_info)
      call set_dependency(form_r12_z4int,op_r12,tgt_info)
      call set_dependency(form_r12_z4int,op_ham,tgt_info)
      call set_dependency(form_r12_z4int,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_z4int,0,'Z4')
      call set_rule(form_r12_z4int,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,6,1,
     &              parameters,2,tgt_info)

      ! Formal definition of K4
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_k4int
      labels(2) = op_k4_inter
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target(form_r12_k4int,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_k4int,op_k4_inter,tgt_info)
      call set_dependency(form_r12_k4int,op_r12,tgt_info)
      call set_dependency(form_r12_k4int,op_ham,tgt_info)
      call set_dependency(form_r12_k4int,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_k4int,0,'K4')
      call set_rule(form_r12_k4int,ttype_frm,DEF_R12INTM_FORMAL,
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
      call set_dependency(fopt_r12_vcabs,mel_gr,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_v_def,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_gintx,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_vcabs,ttype_frm,OPTIMIZE,
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
      call set_dependency(fopt_r12_xcabs,mel_ff,tgt_info)
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
      call set_dependency(fopt_r12_bcabs,mel_rttr,tgt_info)
      call set_dependency(fopt_r12_bcabs,mel_ham,tgt_info)
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

      ! set P
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_pcabs
      labels(2) = form_r12_pcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_pcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_pcabs,form_r12_pcabs,tgt_info)
      call set_dependency(fopt_r12_pcabs,mel_p_def,tgt_info)
      call set_dependency(fopt_r12_pcabs,mel_v_def,tgt_info)
      call set_dependency(fopt_r12_pcabs,mel_rint,tgt_info)      
      call set_dependency(fopt_r12_pcabs,mel_gintx,tgt_info)      
      call set_dependency(fopt_r12_pcabs,mel_ffg,tgt_info)      
      call set_dependency(fopt_r12_pcabs,mel_gr,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_pcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! set P3F
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_p3fcabs
      labels(2) = form_r12_p3fcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_p3fcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_p3fcabs,form_r12_p3fcabs,tgt_info)
      call set_dependency(fopt_r12_p3fcabs,mel_p3f_def,tgt_info)
      call set_dependency(fopt_r12_p3fcabs,mel_ham,tgt_info)
      call set_dependency(fopt_r12_p3fcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_p3fcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! set P3G
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_p3gcabs
      labels(2) = form_r12_p3gcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_p3gcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_p3gcabs,form_r12_p3gcabs,tgt_info)
      call set_dependency(fopt_r12_p3gcabs,mel_p3g_def,tgt_info)
      call set_dependency(fopt_r12_p3gcabs,mel_v_def,tgt_info)
      call set_dependency(fopt_r12_p3gcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_p3gcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! set Z
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_zcabs
      labels(2) = form_r12_zcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_zcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_zcabs,form_r12_zcabs,tgt_info)
      call set_dependency(fopt_r12_zcabs,mel_z_def,tgt_info)
      call set_dependency(fopt_r12_zcabs,mel_rint,tgt_info)      
      call set_dependency(fopt_r12_zcabs,mel_gintx,tgt_info)      
      call set_dependency(fopt_r12_zcabs,mel_ff,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_zcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

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
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_ccabs,ttype_frm,OPTIMIZE,
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

      ! R12^2*G integrals
      call add_target(mel_ffg,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ffg,op_ffg,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ffg
      labels(2) = op_ffg
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_ffg,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ffg
      call import_parameters(-1,parameters,env_type)
      call set_rule(mel_ffg,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

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

c dbg
      ! X-test
      call add_target(mel_x_test_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_x_test_def,op_x_test,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_x_test
      labels(2) = op_x_test
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_x_test_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
c dbg

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

      ! P-list
      call add_target(mel_p_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_p_def,op_p_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_p_inter
      labels(2) = op_p_inter
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_p_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

c dbg
      ! P-test
      call add_target(mel_p_test_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_p_test_def,op_p_test,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_p_test
      labels(2) = op_p_test
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_p_test_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
c dbg

      ! P3F-list
      call add_target(mel_p3f_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_p3f_def,op_p3f_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_p3f_inter
      labels(2) = op_p3f_inter
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_p3f_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! P3G-list
      call add_target(mel_p3g_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_p3g_def,op_p3g_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_p3g_inter
      labels(2) = op_p3g_inter
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_p3g_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! Z-list
      call add_target(mel_z_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_z_def,op_z_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_z_inter
      labels(2) = op_z_inter
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_z_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

c dbg
      ! Z-test
      call add_target(mel_z_test_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_z_test_def,op_z_test,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_z_test
      labels(2) = op_z_test
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_z_test_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
c dbg

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

      ! B^-1 for "diagonal"
      call add_target(mel_b_inv,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_b_inv,op_b_inv,tgt_info)
      call set_dependency(mel_b_inv,eval_r12_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_b_inv
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
c dbg
c      call add_target(eval_r12_inter,ttype_gen,.true.,tgt_info)
c dbg
      call set_dependency(eval_r12_inter,mel_ham,tgt_info)
      call set_dependency(eval_r12_inter,mel_rint,tgt_info)
      call set_dependency(eval_r12_inter,mel_gintx,tgt_info)
      call set_dependency(eval_r12_inter,mel_ttr,tgt_info)
      call set_dependency(eval_r12_inter,mel_ff,tgt_info)
      call set_dependency(eval_r12_inter,mel_gr,tgt_info)
      call set_dependency(eval_r12_inter,mel_v_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_x_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_b_def,tgt_info)
      if(is_keyword_set('method.CC').gt.0.and..not.truncate)then
        call set_dependency(eval_r12_inter,mel_ffg,tgt_info)
        call set_dependency(eval_r12_inter,mel_p_def,tgt_info)
        call set_dependency(eval_r12_inter,mel_p3g_def,tgt_info)
        call set_dependency(eval_r12_inter,fopt_r12_pcabs,tgt_info)
c        call set_dependency(eval_r12_inter,mel_z_def,tgt_info)
c        call set_dependency(eval_r12_inter,fopt_r12_zcabs,tgt_info)
      endif
      call set_dependency(eval_r12_inter,fopt_r12_vcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_xcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_bcabs,tgt_info)
      if (ansatz.ne.1)then
        call set_dependency(eval_r12_inter,fopt_r12_ccabs,tgt_info)
c        call set_dependency(eval_r12_inter,fopt_r12_p3fcabs,tgt_info)
        if(is_keyword_set('method.CC').gt.0.and..not.truncate)then
          call set_dependency(eval_r12_inter,fopt_r12_p3gcabs,tgt_info)
        endif
      endif
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_vcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      labels(1) = fopt_r12_xcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      if (ansatz.ne.1) then
        labels(1) = fopt_r12_ccabs
        call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

c        labels(1) = fopt_r12_p3fcabs
c        call set_rule(eval_r12_inter,ttype_opme,EVAL,
c     &     labels,1,0,
c     &     parameters,0,tgt_info)

        if(is_keyword_set('method.CC').gt.0.and..not.truncate)then
          labels(1) = fopt_r12_p3gcabs
          call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)
        endif
      end if

      labels(1) = fopt_r12_bcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      if(is_keyword_set('method.CC').gt.0.and..not.truncate)then
        labels(1) = fopt_r12_pcabs
        call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &       labels,1,0,
     &       parameters,0,tgt_info)

c        labels(1) = fopt_r12_zcabs
c        call set_rule(eval_r12_inter,ttype_opme,EVAL,
c     &       labels,1,0,
c     &       parameters,0,tgt_info)
      endif

c dbg
      ! set X-test opt
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_x_test
      labels(2) = form_x_test
      ncat = 1
      nint = 0
      call add_target(fopt_x_test,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_x_test,form_x_test,tgt_info)
      call set_dependency(fopt_x_test,mel_x_test_def,tgt_info)
      call set_dependency(fopt_x_test,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_x_test,ttype_frm,OPTIMIZE,
     &     labels,ncat+nint+1,1,
     &     parameters,1,tgt_info)

c      call add_target(eval_x_test,ttype_gen,.true.,tgt_info)
      call add_target(eval_x_test,ttype_gen,.false.,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_x_test
      call set_dependency(eval_x_test,mel_rint,tgt_info)
      call set_dependency(eval_x_test,mel_x_test_def,tgt_info)
      call set_dependency(eval_x_test,fopt_x_test,tgt_info)
      call set_rule(eval_x_test,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      ! set P-test opt
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_p_test
      labels(2) = form_p_test
      ncat = 1
      nint = 0
      call add_target(fopt_p_test,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_p_test,form_p_test,tgt_info)
      call set_dependency(fopt_p_test,mel_p_test_def,tgt_info)
      call set_dependency(fopt_p_test,mel_rint,tgt_info)      
      call set_dependency(fopt_p_test,mel_gintx,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_p_test,ttype_frm,OPTIMIZE,
     &     labels,ncat+nint+1,1,
     &     parameters,1,tgt_info)

c      call add_target(eval_p_test,ttype_gen,.true.,tgt_info)
      call add_target(eval_p_test,ttype_gen,.false.,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_p_test
      call set_dependency(eval_p_test,mel_rint,tgt_info)
      call set_dependency(eval_p_test,mel_gintx,tgt_info)
      call set_dependency(eval_p_test,mel_p_test_def,tgt_info)
      call set_dependency(eval_p_test,fopt_p_test,tgt_info)
      call set_rule(eval_p_test,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      ! set Z-test opt
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_z_test
      labels(2) = form_z_test
      ncat = 1
      nint = 0
      call add_target(fopt_z_test,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_z_test,form_z_test,tgt_info)
      call set_dependency(fopt_z_test,mel_z_test_def,tgt_info)
      call set_dependency(fopt_z_test,mel_rint,tgt_info)      
      call set_dependency(fopt_z_test,mel_gintx,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_z_test,ttype_frm,OPTIMIZE,
     &     labels,ncat+nint+1,1,
     &     parameters,1,tgt_info)

c      call add_target(eval_z_test,ttype_gen,.true.,tgt_info)
      call add_target(eval_z_test,ttype_gen,.false.,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_z_test
      call set_dependency(eval_z_test,mel_rint,tgt_info)
      call set_dependency(eval_z_test,mel_gintx,tgt_info)
      call set_dependency(eval_z_test,mel_z_test_def,tgt_info)
      call set_dependency(eval_z_test,fopt_z_test,tgt_info)
      call set_rule(eval_z_test,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
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
