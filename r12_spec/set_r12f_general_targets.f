*----------------------------------------------------------------------*
      subroutine set_r12f_general_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
*     set targets needed in all kinds of R12 calculations based on 
*     fixed geminals
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
     &     min_rank, max_rank, ansatz, n_pp, ndef,
     &     min_rank_tp, min_rank_tpp,
     &     isim, ncat, nint, icnt, nlab, irank, idef,
     &     isym, ms, msc, sym_arr(8), extend, r12op,
     &     occ_def(ngastp,2,20),
     &     ntp_min, ntp_max, ntpp_min, ntpp_max
      logical ::
     &     needed, r12fix, set_tp, set_tpp
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

      msc = +1  ! assuming closed shell
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
      call get_argument_value('method.R12','min_tp',ival=min_rank_tp)
      call get_argument_value('method.R12','min_tpp',ival=min_rank_tpp)
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','extend',ival=extend)
      call get_argument_value('method.R12','r12op',ival=r12op)

      n_pp = 0  ! number of particle-particle interaction in R12
      set_tp = .false.
      set_tpp = .false.
      select case(extend)
      case(1) 
        ! T1'
        set_tp = .true.
        ntp_min=1
        ntp_max=1
        n_pp=1
      case(2)
        ! T0'
        set_tp = .true.
        ntp_min=0
        ntp_max=0
        n_pp=0
      case(3)
        ! T0' + T1'
        set_tp = .true.
        ntp_min=0
        ntp_max=1
        n_pp=1
      case(4)
        ! T1' + T2' (for CC)
        set_tp = .true.
        ntp_min=1
        ntp_max=2
        n_pp=1
      case(5)
        ! T1' + T2'
        set_tp = .true.
        ntp_min=1
        ntp_max=2
        n_pp=2
      case(6)
        ! T1' 
        set_tp = .true.
        ntp_min=1
        ntp_max=1
        n_pp=2
      case default 
        set_tp = .false.
        ntp_min=0
        ntp_max=0
        n_pp=0
      end select

      if (r12op.gt.0.and.extend.gt.0)
     &     call quit(1,'set_r12f_general_targets',
     &     'use either r12op or extend')
      if (r12op.gt.0) then
        set_tp  = .false.
        set_tpp = .false.
        ntp_min=0
        ntp_max=0
        ntpp_min=0
        ntpp_max=0
      end if
      select case(r12op)
      case(1)
        ! T' operators (singly p-connected to R12)
        set_tp = .true.
        ntp_min=min_rank_tp
        ntp_max=max_rank-1
        n_pp=1
      case(2)
        ! T'' operators (doubly p-connected to R12)
        set_tpp = .true.
        ntpp_min=min_rank_tpp
        ntpp_max=max_rank
        n_pp=2
      case(3,4)
        ! T' + T'' operators
        set_tp = .true.
        ntp_min=min_rank_tp
        ntp_max=max_rank-1
        n_pp=1
        set_tpp = .true.
        ntpp_min=min_rank_tpp
        ntpp_max=max_rank
        n_pp=2
      end select

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

      select case(trim(K_appr))
      case('none')
        write(luout,*) 'no approximations wrt. Xchange made'
      case('HY1')
        write(luout,*) 'Y contribution omitted'
        approx(8:10) = 'HY1'
      case('HY2')
        write(luout,*) 'Y contribution approx with 1 CABS index'
        approx(8:10) = 'HY2'
      case default
        call quit(0,'set_r12_general_targets',
     &       'K_appr unknown: "'//trim(K_appr)//'"')
      end select

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! the formal R12 geminal: P12 r12 X12
      ! need to be extended ...
      call add_target(op_r12,ttype_op,.false.,tgt_info)
c      min_rank = 2  ! 1 is a possibility 
      call r12gem_parameters(-1,parameters,
     &                   n_pp,2,2,ansatz)
      call set_rule(op_r12,ttype_op,DEF_R12GEMINAL,
     &              op_r12,1,1,
     &              parameters,1,tgt_info)

      if (set_tp) then
        ! T1' operators for extended CC/MP2-F12.
        call add_target(op_cex,ttype_op,.false.,tgt_info)
c        occ_def = 0
c        ndef = ntp_max-ntp_min+1
c        irank = ntp_min-1
c        do idef = 1, ndef
c          irank = irank+1
c          occ_def(IPART,1,(idef-1)*2+1) = irank-1
c          occ_def(IHOLE,2,(idef-1)*2+1) = irank
c          occ_def(IPART,1,(idef-1)*2+2) = 1
c        end do
c        call op_from_occ_parameters(-1,parameters,2,
c     &       occ_def,ndef,2,ndef)
c        call set_rule(op_cex,ttype_op,DEF_OP_FROM_OCC,
c     &       op_cex,1,1,
c     &       parameters,2,tgt_info)
        call xop_parameters(-1,parameters,
     &       .false.,ntp_min,ntp_max,0,ntp_max+2)
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

        call add_target(op_diar12,ttype_op,.false.,tgt_info)
        call set_dependency(op_diar12,op_cex,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_cex,.false.) ! <- dagger=.false.
        call set_rule(op_diar12,ttype_op,CLONE_OP,
     &                op_diar12,1,1,
     &                parameters,1,tgt_info)

      endif
c dbg
      print *,'in r12f_gen: set_tp = ',set_tp
c dbg

      if (set_tpp) then
        ! T1'' operators for extended CC/MP2-F12.
        call add_target(op_cexx,ttype_op,.false.,tgt_info)
        call xop_parameters(-1,parameters,
     &       .false.,ntpp_min,ntpp_max,0,ntpp_max+2)
        call set_rule(op_cexx,ttype_op,DEF_EXCITATION,
     &                op_cexx,1,1,
     &                parameters,1,tgt_info)

        ! The Lagrangian multipliers.
        call add_target(op_cexxbar,ttype_op,.false.,tgt_info)
        call set_dependency(op_cexxbar,op_cexx,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_cexx,.true.) ! <- dagger=.true.
        call set_rule(op_cexxbar,ttype_op,CLONE_OP,
     &                op_cexxbar,1,1,
     &                parameters,1,tgt_info)

      endif

      ! Now: the operators associated with the actual R12 integrals:
      !  <pq'|r12|ij> 
      call add_target(op_rint,ttype_op,.false.,tgt_info)
      
      call r12int_parameters(-1,parameters,
c     &     .false.,min_rank,2,0,2)
     &      n_pp,min_rank,2,0,4) ! 4: two externals for 2-el ops
      call set_rule(op_rint,ttype_op,DEF_R12INT,
     &              op_rint,1,1,
     &              parameters,1,tgt_info)
      
      ! (pq)_frozen/(pq)_ae block of 2e-Hamiltonian
      call add_target(op_g_x,ttype_op,.false.,tgt_info)
c      min_rank = 2 
      call r12int_parameters(-1,parameters,
     &     2,2,2,0,2)
      call set_rule(op_g_x,ttype_op,DEF_R12INT,
     &              op_g_x,1,1,
     &              parameters,1,tgt_info)
      
      ! commutator integrals <ab|[T1+T2,r12]|ij>
      call add_target(op_ttr,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
     &     n_pp,min_rank,2,0,2)
      call set_rule(op_ttr,ttype_op,DEF_R12INT,
     &              op_ttr,1,1,
     &              parameters,1,tgt_info)
      
      ! A(f+k) modified integrals r12bar
      call add_target(op_rintbar,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
     &     n_pp,min_rank,2,0,2)
      call set_rule(op_rintbar,ttype_op,DEF_R12INT,
     &              op_rintbar,1,1,
     &              parameters,1,tgt_info)
      
      ! C(f+k) modified integrals r12bar+
      call add_target(op_rdagbar,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
     &     n_pp,min_rank,2,0,2)
      call set_rule(op_rdagbar,ttype_op,DEF_R12INT,
     &              op_rdagbar,1,1,
     &              parameters,1,tgt_info)
      
      ! C(f) modified integrals r12breve
      call add_target(op_rintbreve,ttype_op,.false.,tgt_info)
      call r12int_parameters(-1,parameters,
     &     n_pp,min_rank,2,0,2)
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
      occ_def = 0
      if (n_pp.ge.0) then
        ndef = 1
        occ_def(IHOLE,1,1) = 2
        occ_def(IHOLE,2,2) = 2
      end if
      if (n_pp.ge.1) then
        ndef = 4
        occ_def(IHOLE,1,3) = 2
        occ_def(IHOLE,2,4) = 1
        occ_def(IPART,2,4) = 1

        occ_def(IHOLE,1,5) = 1
        occ_def(IPART,1,5) = 1
        occ_def(IHOLE,2,6) = 2

        occ_def(IHOLE,1,7) = 1
        occ_def(IPART,1,7) = 1
        occ_def(IHOLE,2,8) = 1
        occ_def(IPART,2,8) = 1
      end if      
      if (n_pp.ge.2) then
        ndef = 9
        occ_def(IHOLE,1,9) = 2
        occ_def(IPART,2,10) = 2

        occ_def(IPART,1,11) = 2
        occ_def(IHOLE,2,12) = 2

        occ_def(IHOLE,1,13) = 1
        occ_def(IPART,1,13) = 1
        occ_def(IPART,2,14) = 2

        occ_def(IPART,1,15) = 2
        occ_def(IHOLE,2,16) = 1
        occ_def(IPART,2,16) = 1

        occ_def(IPART,1,17) = 2
        occ_def(IPART,2,18) = 2
      end if
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,2,ndef)
      call set_rule(op_rttr,ttype_op,DEF_OP_FROM_OCC,
     &              op_rttr,1,1,
     &              parameters,2,tgt_info)

      ! (G.R)^{ij}_{pq}
      call add_target(op_gr,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! n_pp == 0:
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,2,2) = 2

      occ_def(IHOLE,1,3) = 1
      occ_def(IPART,1,3) = 1
      occ_def(IHOLE,2,4) = 2

      occ_def(IPART,1,5) = 2
      occ_def(IHOLE,2,6) = 2
      ! n_pp == 1:
      occ_def(IHOLE,1,7) = 2
      occ_def(IHOLE,2,8) = 1
      occ_def(IPART,2,8) = 1

      occ_def(IHOLE,1,9) = 1
      occ_def(IPART,1,9) = 1
      occ_def(IHOLE,2,10) = 1
      occ_def(IPART,2,10) = 1

      occ_def(IPART,1,11) = 2
      occ_def(IHOLE,2,12) = 1
      occ_def(IPART,2,12) = 1
      ! n_pp == 2:
      occ_def(IHOLE,1,13) = 2
      occ_def(IPART,2,14) = 2

c      occ_def(IHOLE,1,15) = 1
c      occ_def(IPART,1,15) = 1
c      occ_def(IPART,2,16) = 2

      occ_def(IPART,1,17) = 2
      occ_def(IPART,2,18) = 2

      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,3*(n_pp+1),2,3*(n_pp+1))
      call set_rule(op_gr,ttype_op,DEF_OP_FROM_OCC,
     &              op_gr,1,1,
     &              parameters,2,tgt_info)

      ! R12^{2} integrals
      call add_target(op_ff,ttype_op,.false.,tgt_info)
      call set_dependency(op_ff,op_rttr,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rttr,.false.) ! <- dagger=.false.
      call set_rule(op_ff,ttype_op,CLONE_OP,
     &     op_ff,1,1,
     &     parameters,1,tgt_info)

      ! {R12^2}BAR integrals
      call add_target(op_ffbar,ttype_op,.false.,tgt_info)
      call set_dependency(op_ffbar,op_rttr,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rttr,.false.) ! <- dagger=.false.
      call set_rule(op_ffbar,ttype_op,CLONE_OP,
     &              op_ffbar,1,1,
     &              parameters,1,tgt_info)

            
      ! V^{ij}_{pq}
      call add_target(op_v_inter,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! for n_pp >= 0
      if (n_pp.ge.0) then
        ndef = 6
        ! block 1 -> scalar
        occ_def(IHOLE,1,2) = 1
        occ_def(IHOLE,2,2) = 1

        occ_def(IPART,1,3) = 1
        occ_def(IHOLE,2,3) = 1

        occ_def(IHOLE,1,4) = 2
        occ_def(IHOLE,2,4) = 2

        occ_def(IHOLE,1,5) = 1
        occ_def(IPART,1,5) = 1
        occ_def(IHOLE,2,5) = 2

        occ_def(IPART,1,6) = 2
        occ_def(IHOLE,2,6) = 2
      end if
      ! for n_pp >= 1
      if (n_pp.ge.1) then
        ndef = 11
        occ_def(IHOLE,1,7) = 1
        occ_def(IPART,2,7) = 1

        occ_def(IPART,1,8) = 1
        occ_def(IPART,2,8) = 1

        occ_def(IHOLE,1,9) = 2
        occ_def(IHOLE,2,9) = 1
        occ_def(IPART,2,9) = 1

        occ_def(IHOLE,1,10) = 1
        occ_def(IPART,1,10) = 1
        occ_def(IHOLE,2,10) = 1
        occ_def(IPART,2,10) = 1

        occ_def(IPART,1,11) = 2
        occ_def(IHOLE,2,11) = 1
        occ_def(IPART,2,11) = 1
      end if
      ! for n_pp >= 2
      if (n_pp.ge.2) then
        ndef = 14
        occ_def(IHOLE,1,12) = 2
        occ_def(IPART,2,12) = 2

        occ_def(IHOLE,1,13) = 1
        occ_def(IPART,1,13) = 1
        occ_def(IPART,2,13) = 2

        occ_def(IPART,1,14) = 2
        occ_def(IPART,2,14) = 2
      end if
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,ndef)
      call set_rule(op_v_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_v_inter,1,1,
     &              parameters,2,tgt_info)

      ! V' intermediate
      call add_target(op_vp_inter,ttype_op,.false.,tgt_info)
      call set_dependency(op_vp_inter,op_v_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_v_inter,.false.) ! <- dagger=.false.
      call set_rule(op_vp_inter,ttype_op,CLONE_OP,
     &              op_vp_inter,1,1,
     &              parameters,1,tgt_info)
            
      ! B intermediate
      occ_def = 0
      if (n_pp.ge.0) then
        ndef = 1
      end if
      if (n_pp.ge.1) then
        ndef = 7
        occ_def(IHOLE,1,2) = 1
        occ_def(IPART,2,2) = 1

        occ_def(IPART,1,3) = 1
        occ_def(IHOLE,2,3) = 1

        occ_def(IPART,1,4) = 1
        occ_def(IPART,2,4) = 1

        occ_def(IHOLE,1,5) = 1
        occ_def(IPART,1,5) = 1
        occ_def(IHOLE,2,5) = 2

        occ_def(IHOLE,1,6) = 2
        occ_def(IHOLE,2,6) = 1
        occ_def(IPART,2,6) = 1

        occ_def(IHOLE,1,7) = 1
        occ_def(IPART,1,7) = 1
        occ_def(IHOLE,2,7) = 1
        occ_def(IPART,2,7) = 1
      end if
      if (n_pp.ge.2) then
        occ_def(IHOLE,1,ndef+1) = 2
        occ_def(IPART,2,ndef+1) = 2

        occ_def(IPART,1,ndef+2) = 2
        occ_def(IHOLE,2,ndef+2) = 2

        occ_def(IHOLE,1,ndef+3) = 1
        occ_def(IPART,1,ndef+3) = 1
        occ_def(IPART,2,ndef+3) = 2

        occ_def(IPART,1,ndef+4) = 2
        occ_def(IHOLE,2,ndef+4) = 1
        occ_def(IPART,2,ndef+4) = 1

        occ_def(IPART,1,ndef+5) = 2
        occ_def(IPART,2,ndef+5) = 2
        ndef = ndef + 5
      end if
      call add_target(op_b_inter,ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,ndef)
      call set_rule(op_b_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_b_inter,1,1,
     &              parameters,2,tgt_info)

      ! Bhole intermediate
      call add_target(op_bh_inter,ttype_op,.false.,tgt_info)
      call set_dependency(op_bh_inter,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
      call set_rule(op_bh_inter,ttype_op,CLONE_OP,
     &              op_bh_inter,1,1,
     &              parameters,1,tgt_info)

      ! X intermediate
      call add_target(op_x_inter,ttype_op,.false.,tgt_info)
      call set_dependency(op_x_inter,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
      call set_rule(op_x_inter,ttype_op,CLONE_OP,
     &              op_x_inter,1,1,
     &              parameters,1,tgt_info)

      ! X' intermediate
      call add_target(op_xp_inter,ttype_op,.false.,tgt_info)
      ndef = 4
      occ_def = 0
      occ_def(IPART,1,1) = 1
      occ_def(IPART,2,1) = 1

      occ_def(IHOLE,1,2) = 1
      occ_def(IPART,1,2) = 1
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1

      occ_def(IHOLE,1,3) = 1
      occ_def(IPART,1,3) = 1
      occ_def(IPART,2,3) = 2

      occ_def(IPART,1,4) = 2
      occ_def(IHOLE,2,4) = 1
      occ_def(IPART,2,4) = 1

      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,ndef)
      call set_rule(op_xp_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_xp_inter,1,1,
     &              parameters,2,tgt_info)

      ! Xh intermediate
      call add_target(op_xh_inter,ttype_op,.false.,tgt_info)
      ndef = 4
      occ_def = 0
      occ_def(IHOLE,1,1) = 1
      occ_def(IHOLE,2,2) = 1

      occ_def(IHOLE,1,3) = 2
      occ_def(IHOLE,2,4) = 1
      occ_def(IPART,2,4) = 1

      occ_def(IHOLE,1,5) = 1
      occ_def(IPART,1,5) = 1
      occ_def(IHOLE,2,6) = 2

      occ_def(IHOLE,1,7) = 1
      occ_def(IPART,1,7) = 1
      occ_def(IHOLE,2,8) = 1
      occ_def(IPART,2,8) = 1

      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,2,ndef*2)
      call set_rule(op_xh_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_xh_inter,1,1,
     &              parameters,2,tgt_info)

      ! C intermediate
      call add_target(op_c_inter,ttype_op,.false.,tgt_info)
      occ_def = 0
      ! n_pp == 0:
      occ_def(IPART,1,1) = 2
      occ_def(IHOLE,2,1) = 2
      ! n_pp == 1:
      occ_def(IPART,1,2) = 2
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1
      ! n_pp == 2:
      occ_def(IPART,1,3) = 2
      occ_def(IPART,2,3) = 2

      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,n_pp+1,1,6)
      call set_rule(op_c_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_c_inter,1,1,
     &              parameters,2,tgt_info)

c      ! Z^{ijp}_{klm} intermediate (needed for CC)
cc      call add_target(op_z_inter,ttype_op,.false.,tgt_info)
cc dbg
c      call add_target(op_z_inter,ttype_op,.true.,tgt_info)
cc dbg
c      occ_def = 0
c      occ_def(IHOLE,1,1) = 3
c      occ_def(IHOLE,2,1) = 3
c      occ_def(IHOLE,1,2) = 3
c      occ_def(IHOLE,2,2) = 2
c      occ_def(IPART,2,2) = 1
c      ndef = 2
c      call op_from_occ_parameters(-1,parameters,2,
c     &     occ_def,ndef,1,ndef)
c      call set_rule(op_z_inter,ttype_op,DEF_OP_FROM_OCC,
c     &              op_z_inter,1,1,
c     &              parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*
      ! formal definition of V
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_vint
      labels(2) = op_v_inter
      labels(3) = op_r12
      labels(4) = op_ham
      call add_target(form_r12_vint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_vint,op_v_inter,tgt_info)
      call set_dependency(form_r12_vint,op_ham,tgt_info)
      call set_dependency(form_r12_vint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vint,0,'V')
      call set_rule(form_r12_vint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to V
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_vcabs
      labels(2) = op_v_inter
      labels(3) = op_g_x !op_ham
      labels(4) = op_rint
      labels(5) = op_gr   
      ! F12: op_gr
      call add_target(form_r12_vcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_vcabs,op_v_inter,tgt_info)
      call set_dependency(form_r12_vcabs,op_gr,tgt_info)
      call set_dependency(form_r12_vcabs,op_g_x,tgt_info)
      call set_dependency(form_r12_vcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vcabs,ansatz,'V '//approx)
      call set_rule(form_r12_vcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! formal definition of X
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xint
      labels(2) = op_x_inter
      labels(3) = op_r12
      labels(4) = op_r12
      call add_target(form_r12_xint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xint,op_x_inter,tgt_info)
      call set_dependency(form_r12_xint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xint,0,'X')
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

      ! formal definition of Xh
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xhint
      labels(2) = op_xh_inter
      labels(3) = op_r12
      labels(4) = op_r12
      call add_target(form_r12_xhint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xhint,op_xh_inter,tgt_info)
      call set_dependency(form_r12_xhint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xhint,0,'X ')
      call set_rule(form_r12_xhint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)
c dbg
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule(form_r12_xhint,ttype_frm,PRINT_FORMULA,
     &              labels,2,1,
     &              parameters,2,tgt_info)
c dbg

      ! CABS approximation to Xh
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xhcabs
      labels(2) = op_xh_inter
      labels(3) = op_rint
      labels(4) = op_rint
      labels(5) = op_ff
      call add_target(form_r12_xhcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xhcabs,op_xh_inter,tgt_info)
      call set_dependency(form_r12_xhcabs,op_ff,tgt_info)
      call set_dependency(form_r12_xhcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xcabs,ansatz,'XH'//approx)
      call set_rule(form_r12_xhcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      ! formal definition of X'
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xpint
      labels(2) = op_xp_inter
      labels(3) = op_r12
      labels(4) = op_r12
      call add_target(form_r12_xpint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xpint,op_xp_inter,tgt_info)
      call set_dependency(form_r12_xpint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xpint,0,'X''')
      call set_rule(form_r12_xpint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to X'
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_xpcabs
      labels(2) = form_r12_xpint
      labels(3) = op_r12
      labels(4) = op_rint
      call add_target(form_r12_xpcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xcabs,op_xp_inter,tgt_info)
      call set_dependency(form_r12_xcabs,op_r12,tgt_info)
      call set_dependency(form_r12_xcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xpcabs,1,'---')
      call set_rule(form_r12_xpcabs,ttype_frm,REPLACE,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! formal definition of B
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_bint
      labels(2) = op_b_inter
      labels(3) = op_r12
      labels(4) = op_ham
      call add_target(form_r12_bint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_bint,op_b_inter,tgt_info)
      call set_dependency(form_r12_bint,op_ham,tgt_info)
      call set_dependency(form_r12_bint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_bint,0,'Bp')
      call set_rule(form_r12_bint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
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

      ! formal definition of Bh
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_bhint
      labels(2) = op_bh_inter
      labels(3) = op_r12
      labels(4) = op_ham
      call add_target(form_r12_bhint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_bhint,op_bh_inter,tgt_info)
      call set_dependency(form_r12_bhint,op_r12,tgt_info)
      call set_dependency(form_r12_bhint,op_ham,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'title missing',0,'Bh')
      call set_rule(form_r12_bhint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to Bh
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_bhcabs
      labels(2) = op_bh_inter
      labels(3) = op_rint
      labels(4) = op_rint
      labels(5) = op_ff
      labels(6) = '-'
      labels(7) = op_ham
      call add_target(form_r12_bhcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_bhcabs,op_x_inter,tgt_info)
      call set_dependency(form_r12_bhcabs,op_ff,tgt_info)
      call set_dependency(form_r12_bhcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'title missing',ansatz,'BH'//approx)
      call set_rule(form_r12_bhcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,7,1,
     &              parameters,2,tgt_info)

      ! formal definition of C
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_cint
      labels(2) = op_c_inter
      labels(3) = op_r12
      labels(4) = op_ham
      call add_target(form_r12_cint,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_cint,op_c_inter,tgt_info)
      call set_dependency(form_r12_cint,op_r12,tgt_info)
      call set_dependency(form_r12_cint,op_ham,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_cint,0,'C')
      call set_rule(form_r12_cint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to C
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_ccabs
      labels(2) = op_c_inter
      if (approx(1:1).eq.'B') then
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
      call set_rule(form_r12_ccabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,nlab,1,
     &              parameters,2,tgt_info)

c      ! Formal definition of Z.
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_r12_zint
c      labels(2) = op_z_inter
c      labels(3) = op_r12
c      labels(4) = op_ham
c      labels(5) = op_r12
c      call add_target(form_r12_zint,ttype_frm,.true.,tgt_info)
c      call set_dependency(form_r12_zint,op_b_inter,tgt_info)
c      call set_dependency(form_r12_zint,op_ham,tgt_info)
c      call set_dependency(form_r12_zint,op_r12,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_r12_zint,0,'Z')
c      call set_rule(form_r12_zint,ttype_frm,DEF_R12INTM_FORMAL,
c     &              labels,5,1,
c     &              parameters,2,tgt_info)

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

      ! set Xh
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_xhcabs
      labels(2) = form_r12_xhcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_xhcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_xhcabs,form_r12_xhcabs,tgt_info)
      call set_dependency(fopt_r12_xhcabs,mel_xh_def,tgt_info)
      call set_dependency(fopt_r12_xhcabs,mel_ff,tgt_info)
      call set_dependency(fopt_r12_xhcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_xhcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! set X'
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_xpcabs
      labels(2) = form_r12_xpcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_xpcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_xpcabs,form_r12_xpcabs,tgt_info)
      call set_dependency(fopt_r12_xpcabs,mel_xp_def,tgt_info)
      call set_dependency(fopt_r12_xpcabs,mel_ff,tgt_info)
      call set_dependency(fopt_r12_xpcabs,mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_xpcabs,ttype_frm,OPTIMIZE,
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

      ! set Bhole
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_bhcabs
      labels(2) = form_r12_bhcabs
      ncat = 1
      nint = 0
      call add_target(fopt_r12_bhcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_r12_bhcabs,form_r12_bhcabs,tgt_info)
      call set_dependency(fopt_r12_bhcabs,mel_bh_def,tgt_info)
      call set_dependency(fopt_r12_bhcabs,mel_ff,tgt_info)
      call set_dependency(fopt_r12_bhcabs,mel_rint,tgt_info)      
      call set_dependency(fopt_r12_bhcabs,mel_ham,tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_bhcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
c dbg
c         labels(1) = fopt_r12_bhcabs
c         labels(2) = fopt_r12_bhcabs
c         call modify_parameters(-1,
c     &       parameters,4,(/2,2,2,1/),4)
c         call set_rule(fopt_r12_bhcabs,ttype_frm,MODIFY_FACTORIZATION,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)        
c dbg

      ! set C
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_rint,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rint
      call import_parameters(-1,parameters,'F12_INT',env_type)
      call set_rule(mel_rint,ttype_opme,IMPORT,
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_gintx,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_gintx
      call import_parameters(-1,parameters,'G_INT',env_type)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_ttr,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ttr
      call import_parameters(-1,parameters,'TTF_INT',env_type)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_rttr,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rttr
      call import_parameters(-1,parameters,'FTF_INT',env_type)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_ff,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ff
      call import_parameters(-1,parameters,'FF_INT',env_type)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_ffbar,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ffbar
      call import_parameters(-1,parameters,'FFBAR_INT',env_type)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_gr,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_gr
      call import_parameters(-1,parameters,'FG_INT',env_type)
      call set_rule(mel_gr,ttype_opme,IMPORT,
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_rintbar,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rintbar
      call import_parameters(-1,parameters,'F12BAR_INT',env_type)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_rinttilde,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rinttilde
      call import_parameters(-1,parameters,'F12TLD_INT',env_type)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_rdagbar,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rdagbar
      call import_parameters(-1,parameters,'FDGBAR_INT',env_type)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_rintbreve,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_rintbreve
      call import_parameters(-1,parameters,'F12BRV_INT',env_type)
      call set_rule(mel_rintbreve,ttype_opme,IMPORT,
     &              labels,1,1,
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
     &     msc,0,1,0,0,.false.)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_x_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! Xh-list
      call add_target(mel_xh_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_xh_def,op_xh_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_xh_inter
      labels(2) = op_xh_inter
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_xh_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! X'-list
      call add_target(mel_xp_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_xp_def,op_xp_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_xp_inter
      labels(2) = op_xp_inter
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_xp_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! B-list
      call add_target(mel_b_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_b_def,op_b_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_b_inter
      labels(2) = op_b_inter
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_b_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! Bhole-list
      call add_target(mel_bh_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_bh_def,op_bh_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_bh_inter
      labels(2) = op_bh_inter
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_bh_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! C-list
      call add_target(mel_c_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_c_def,op_c_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_c_inter
      labels(2) = op_c_inter
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_c_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      
      if (set_tp) then
        ! diagonal preconditioner
        call add_target(mel_diar12,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_diar12,op_diar12,tgt_info)
        call set_dependency(mel_diar12,eval_r12_inter,tgt_info)
        call set_dependency(mel_diar12,mel_ham,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_diar12
        labels(2) = op_diar12
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_diar12,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        labels(1) = mel_diar12   ! output
        labels(2) = mel_ham     ! input
        labels(3) = mel_b_inter ! input
        labels(4) = mel_x_inter ! input
        call set_rule(mel_diar12,ttype_opme,PRECONDITIONER,
     &                labels,4,1,
     &                parameters,1,tgt_info)
      end if
*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      ! test
      call add_target(eval_r12_inter,ttype_gen,.false.,tgt_info)
      call set_dependency(eval_r12_inter,mel_ham,tgt_info)
      call set_dependency(eval_r12_inter,mel_rint,tgt_info)
      call set_dependency(eval_r12_inter,mel_gintx,tgt_info)
      call set_dependency(eval_r12_inter,mel_ttr,tgt_info)
      call set_dependency(eval_r12_inter,mel_ff,tgt_info)
      call set_dependency(eval_r12_inter,mel_v_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_x_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_xh_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_b_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_bh_def,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_vcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_xcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_xhcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_bcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_bhcabs,tgt_info)
      if (ansatz.ne.1)
     &     call set_dependency(eval_r12_inter,fopt_r12_ccabs,tgt_info)
c      if (ansatz.ne.1)
c     &     call set_dependency(eval_r12_inter,fopt_r12_xpcabs,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_vcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = fopt_r12_xcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = fopt_r12_xhcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      if (ansatz.ne.1) then
        labels(1) = fopt_r12_ccabs
        call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
c        labels(1) = fopt_r12_xpcabs
c        call set_rule(eval_r12_inter,ttype_opme,EVAL,
c     &     labels,1,0,
c     &     parameters,0,tgt_info)
      end if
      labels(1) = fopt_r12_bcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = fopt_r12_bhcabs
      call set_rule(eval_r12_inter,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      return

      end
