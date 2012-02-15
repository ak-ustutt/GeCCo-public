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
     &     occ_def(ngastp,2,60), ndef, mode, trunc_type,
     &     ninproj, navoid, nconnect, nreplace, 
     &     connect(20), idx_sv(20), iblkmin(20),
     &     iblkmax(20)
      logical ::
     &     needed, r12fix, extend, truncate, frozen,frozen_old, pz_eval,
     &     pf12_trunc
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character(20) ::
     &     approx, F_appr, K_appr, Z_appr, shell_typ

      character(*), intent(in) ::
     &     env_type

*----------------------------------------------------------------------*
      if (iprlvl.gt.0)
     &     write(luout,*) 'setting general targets for R12 ...'

      msc = +1 ! assuming closed shell
*----------------------------------------------------------------------*
*     read input
*----------------------------------------------------------------------*
      ! set approx string
      approx(1:20) = ' '
      F_appr(1:20) = ' '
      K_appr(1:20) = ' '
      Z_appr(1:20) = ' '
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','F_appr',str=F_appr)
      call get_argument_value('method.R12','K_appr',str=K_appr)
      call get_argument_value('method.R12','Z_appr',str=Z_appr)
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','pz_eval',lval=pz_eval)
      call get_argument_value('method.R12','extend',ival=mode)
      call get_argument_value('method.R12','trunc',ival=trunc_type)
      truncate = trunc_type.ge.0
      if (is_keyword_set('method.truncate').gt.0) then
        truncate = is_keyword_set('method.truncate').gt.0
        call get_argument_value('method.truncate','trunc_type',
     &       ival=trunc_type)
      end if
      pf12_trunc = truncate.and.trunc_type.eq.0
      ! Frozen core?
      frozen_old = .false.
      frozen = .false.
      shell_typ(1:len(shell_typ)) = ' '
      if(is_keyword_set('orb_space.shell').gt.0)then
        call get_argument_value('orb_space.shell','type',
     &       str=shell_typ)
        frozen = trim(shell_typ).eq.'frozen'
      endif

      select case(trim(Z_appr))
      case('direct')
        write(luout,*) 'direct RI evaluation of Z intermediate'
        approx(14:17) = 'DRCT'
      case('none','J2K3')
        write(luout,*) 'no approximations to Z intermediate made'
        approx(14:17) = 'J2K3'
      case default
        if (Z_appr(1:1).ne.'J'.or.Z_appr(3:3).ne.'K'.or.
     &      (Z_appr(2:2).ne.'0'.and.
     &       Z_appr(2:2).ne.'1'.and.
     &       Z_appr(2:2).ne.'2').or.
     &      (Z_appr(4:4).ne.'0'.and.
     &       Z_appr(4:4).ne.'1'.and.
     &       Z_appr(4:4).ne.'2'.and.
     &       Z_appr(4:4).ne.'3')) then
          call quit(0,'set_r12_general_targets',
     &       'Z_appr unknown: "'//trim(Z_appr)//'"')
        end if
        write(luout,*) 'approximation to Z intermediate: ',trim(Z_appr)
        approx(14:17) = Z_appr(1:4)
      end select


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
     &                        occ_def,ndef,1,(/     0,     0/),ndef)
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
c      min_rank = 2 
cc      call hop_parameters(-1,parameters,
cc     &     min_rank,2,3,.true.)        ! 1-external at most
cc      call set_rule(op_g_x,ttype_op,DEF_HAMILTONIAN,
cc     &     op_g_x,1,1,
cc     &     parameters,1,tgt_info)
c      call r12int_parameters(-1,parameters,
c     &     2,min_rank,2,0,3)
cc dbg
cc     &     2,min_rank,2,0,5)
cc dbg
c      call set_rule(op_g_x,ttype_op,DEF_R12INT,
c     &              op_g_x,1,1,
c     &              parameters,1,tgt_info)
      ndef = 16
      call set_gxx(occ_def)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/      2,     0/),ndef)
      call set_rule(op_g_x,ttype_op,DEF_OP_FROM_OCC,
     &              op_g_x,1,1,
     &              parameters,2,tgt_info)



      ! ae/ae blocks of 2e-Hamilt. (formal)
      call add_target('G-XX',ttype_op,.false.,tgt_info)
      ndef = 16
      call set_gxx(occ_def)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/      2,      2/),ndef)
      call set_rule('G-XX',ttype_op,DEF_OP_FROM_OCC,
     &              'G-XX',1,1,
     &              parameters,2,tgt_info)

      ! (pq)_frozen/(iq)_ae block of 2e-Hamiltonian (actual)
      call add_target('G-Acore',ttype_op,.false.,tgt_info)
c      min_rank = 2 
      occ_def = 0
      ndef = 9
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,2,1) = 2
      ! 2
      occ_def(IHOLE,1,2) = 2
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1
      ! 3
      occ_def(IHOLE,1,3) = 2
      occ_def(IHOLE,2,3) = 1
      occ_def(IEXTR,2,3) = 1
      ! 4
      occ_def(IHOLE,1,4) = 1
      occ_def(IPART,1,4) = 1
      occ_def(IHOLE,2,4) = 2
      ! 5
      occ_def(IHOLE,1,5) = 1
      occ_def(IPART,1,5) = 1
      occ_def(IHOLE,2,5) = 1
      occ_def(IPART,2,5) = 1
      ! 6
      occ_def(IHOLE,1,6) = 1
      occ_def(IPART,1,6) = 1
      occ_def(IHOLE,2,6) = 1
      occ_def(IEXTR,2,6) = 1
      ! 7
      occ_def(IPART,1,7) = 2
      occ_def(IHOLE,2,7) = 2
      ! 8
      occ_def(IPART,1,8) = 2
      occ_def(IHOLE,2,8) = 1
      occ_def(IPART,2,8) = 1
      ! 9
      occ_def(IPART,1,9) = 2
      occ_def(IHOLE,2,9) = 1
      occ_def(IEXTR,2,9) = 1
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,      2/),ndef)
      call set_rule('G-Acore',ttype_op,DEF_OP_FROM_OCC,
     &              'G-Acore',1,1,
     &              parameters,2,tgt_info)

      ! auxiliary HX|HX list
      call add_target('H-ext',ttype_op,.false.,tgt_info)
c      min_rank = 2 
      occ_def = 0
      ndef = 1
      ! 1
      occ_def(IHOLE,1,1) = 1
      occ_def(IEXTR,1,1) = 1
      occ_def(IHOLE,2,1) = 1
      occ_def(IEXTR,2,1) = 1
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,      2/),ndef)
      call set_rule('H-ext',ttype_op,DEF_OP_FROM_OCC,
     &              'H-ext',1,1,
     &              parameters,2,tgt_info)

      ! (jq)_ae/(ir)_ae block of 2e-Hamiltonian
      call add_target('G-CAcore',ttype_op,.false.,tgt_info)
c      min_rank = 2 
      occ_def = 0
      ndef = 12
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,2,1) = 2
      ! 2
      occ_def(IHOLE,1,2) = 2
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1
      ! 3
      occ_def(IHOLE,1,3) = 2
      occ_def(IHOLE,2,3) = 1
      occ_def(IEXTR,2,3) = 1
      ! 4
      occ_def(IHOLE,1,4) = 1
      occ_def(IPART,1,4) = 1
      occ_def(IHOLE,2,4) = 2
      ! 5
      occ_def(IHOLE,1,5) = 1
      occ_def(IPART,1,5) = 1
      occ_def(IHOLE,2,5) = 1
      occ_def(IPART,2,5) = 1
      ! 6
      occ_def(IHOLE,1,6) = 1
      occ_def(IPART,1,6) = 1
      occ_def(IHOLE,2,6) = 1
      occ_def(IEXTR,2,6) = 1
      ! 7
      occ_def(IHOLE,1,7) = 1
      occ_def(IEXTR,1,7) = 1
      occ_def(IHOLE,2,7) = 2
      ! 8
      occ_def(IHOLE,1,8) = 1
      occ_def(IEXTR,1,8) = 1
      occ_def(IHOLE,2,8) = 1
      occ_def(IPART,2,8) = 1
      ! 9
      occ_def(IHOLE,1,9) = 1
      occ_def(IEXTR,1,9) = 1
      occ_def(IHOLE,2,9) = 1
      occ_def(IEXTR,2,9) = 1
      ! 10
      occ_def(IHOLE,1,10) = 2
      occ_def(IPART,2,10) = 2
      ! 11
      occ_def(IHOLE,1,11) = 1
      occ_def(IPART,1,11) = 1
      occ_def(IPART,2,11) = 2
      ! 12
      occ_def(IHOLE,1,12) = 1
      occ_def(IEXTR,1,12) = 1
      occ_def(IPART,2,12) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/      2,      2/),ndef)
      call set_rule('G-CAcore',ttype_op,DEF_OP_FROM_OCC,
     &              'G-CAcore',1,1,
     &              parameters,2,tgt_info)

      
      ! i,(a/x) block of Fock
      call add_target('F-X',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 2
      occ_def(IHOLE,1,1) = 1
      occ_def(IPART,2,1) = 1
      occ_def(IHOLE,1,2) = 1
      occ_def(IEXTR,2,2) = 1
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/      2,     0/),ndef)
      call set_rule('F-X',ttype_op,DEF_OP_FROM_OCC,
     &              'F-X',1,1,
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
      occ_def = 0
      ndef = 2
      ! 1
      occ_def(IHOLE,1,1) = 1
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,1) = 2
      ! 2
      occ_def(IHOLE,1,2) = 1
      occ_def(IEXTR,1,2) = 1
      occ_def(IHOLE,2,2) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/      2,     0/),ndef)
      call set_rule(op_rintbreve,ttype_op,DEF_OP_FROM_OCC,
     &              op_rintbreve,1,1,
     &              parameters,2,tgt_info)
      
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
     &     occ_def,7,2,(/   0,0,  0,  0/),14)
      call set_rule(op_gr,ttype_op,DEF_OP_FROM_OCC,
     &              op_gr,1,1,
     &              parameters,2,tgt_info)

      ! extended variant (formal)
      call add_target(op_gr_x,ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,7,2,(/   2, 0,   2,     0/),14)
      call set_rule(op_gr_x,ttype_op,DEF_OP_FROM_OCC,
     &              op_gr_x,1,1,
     &              parameters,2,tgt_info)

      ! core list
      ndef = 5
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
      ! 4
      occ_def(IHOLE,1,7) = 2
      occ_def(IHOLE,2,8) = 1
      occ_def(IPART,2,8) = 1
      ! 5
      occ_def(IHOLE,1,9) = 2
      occ_def(IPART,2,10) = 2
      call add_target('G.R-Ccore',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,2,(/  2, 0, 2, 0/),2*ndef)
      call set_rule('G.R-Ccore',ttype_op,DEF_OP_FROM_OCC,
     &              'G.R-Ccore',1,1,
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

      ndef = 1
      if (is_keyword_set('method.CC').gt.0.and.pf12_trunc) ndef = 3
      if (is_keyword_set('method.CC').gt.0.and..not.pf12_trunc) ndef = 5
      
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,2,(/  0,   0,  0,   0/),10)
      call set_rule(op_v_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_v_inter,1,1,
     &              parameters,2,tgt_info)
      
      ! extended variant (formal)
      call add_target(op_v_x,ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,5,2,(/  2, 0,  2,  0/),10)
      call set_rule(op_v_x,ttype_op,DEF_OP_FROM_OCC,
     &              op_v_x,1,1,
     &              parameters,2,tgt_info)

      ! core contributions containing part
      ndef = 3
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

      call add_target('V-Ccore',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,2,(/  2, 0,  2,  0/),2*ndef)
      call set_rule('V-Ccore',ttype_op,DEF_OP_FROM_OCC,
     &              'V-Ccore',1,1,
     &              parameters,2,tgt_info)
      
c dbg
c      call add_target(op_v_test,ttype_op,.false.,tgt_info)
cc      call xop_parameters(-1,parameters,
cc     &     .false.,2,2,0,2)
cc      call set_rule(op_v_inter,ttype_op,DEF_R12INTERM,
cc     &              op_v_inter,1,1,
cc     &              parameters,1,tgt_info)
c      occ_def = 0
c      ! 1
c      occ_def(IHOLE,1,1) = 2
c      occ_def(IHOLE,2,2) = 2
c      ! 2
c      occ_def(IHOLE,1,3) = 1
c      occ_def(IPART,1,3) = 1
c      occ_def(IHOLE,2,4) = 2
c      ! 3
c      occ_def(IPART,1,5) = 2
c      occ_def(IHOLE,2,6) = 2
c      ! 4
c      occ_def(IHOLE,1,7) = 1
c      occ_def(IEXTR,1,7) = 1
c      occ_def(IHOLE,2,8) = 2
c      ! 5
c      occ_def(IPART,1,9) = 1
c      occ_def(IEXTR,1,9) = 1
c      occ_def(IHOLE,2,10) = 2
c
c      call op_from_occ_parameters(-1,parameters,2,
c     &     occ_def,5,2,(/     0,     0/),10)
c      call set_rule(op_v_test,ttype_op,DEF_OP_FROM_OCC,
c     &              op_v_test,1,1,
c     &              parameters,2,tgt_info)
c dbg
      
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
      call set_dependency(op_ff,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     op_b_inter,.false.)  ! <- dagger=.false.
      call set_rule(op_ff,ttype_op,CLONE_OP,
     &                op_ff,1,1,
     &                parameters,1,tgt_info)

c      if(is_keyword_set('method.CC').gt.0.and.(.not.truncate
c     &     .or.(truncate.and.trunc_type.gt.0)))then
        ! we also need this operator
        call add_target('R.R-X',ttype_op,.false.,tgt_info)
        occ_def = 0
        ! 1
        occ_def(IHOLE,1,1)  = 2
        occ_def(IHOLE,2,2)  = 2
        ! 2
        occ_def(IHOLE,1,3)  = 1
        occ_def(IPART,1,3)  = 1
        occ_def(IHOLE,2,4)  = 2
        ! 3
        occ_def(IHOLE,1,5)  = 1
        occ_def(IEXTR,1,5)  = 1
        occ_def(IHOLE,2,6)  = 2
        ! 4
        occ_def(IHOLE,1,7)  = 2
        occ_def(IHOLE,2,8)  = 1
        occ_def(IPART,2,8)  = 1
c        ! 5
        occ_def(IHOLE,1,9)  = 2
        occ_def(IHOLE,2,10) = 1
        occ_def(IEXTR,2,10) = 1
        call op_from_occ_parameters(-1,parameters,2,
c     &       occ_def,5,2,(/     0,      2/),10)
     &       occ_def,5,2,(/      2,     0/),10)
        call set_rule('R.R-X',ttype_op,DEF_OP_FROM_OCC,
     &                'R.R-X',1,1,
     &                parameters,2,tgt_info)
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

c dbg
c      call add_target(op_x_test,ttype_op,.false.,tgt_info)
c      call set_dependency(op_x_test,op_x_inter,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_x_inter,.false.) ! <- dagger=.false.
c      call set_rule(op_x_test,ttype_op,CLONE_OP,
c     &              op_x_test,1,1,
c     &              parameters,1,tgt_info)
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
     &     occ_def,1,1,(/     0,     0/),6)
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
     &     occ_def,1,2,(/     0,     0/),2)
      call set_rule(op_p_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_p_inter,1,1,
     &              parameters,2,tgt_info)

c dbg
c      call add_target(op_p_test,ttype_op,.false.,tgt_info)
c      call set_dependency(op_p_test,op_p_inter,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_p_inter,.false.) ! <- dagger=.false.
c      call set_rule(op_p_test,ttype_op,CLONE_OP,
c     &              op_p_test,1,1,
c     &              parameters,1,tgt_info)
c dbg

      ! R12^{2}*G12 integrals
      call add_target(op_ffg,ttype_op,.false.,tgt_info)
      call set_dependency(op_ffg,op_p_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     op_p_inter,.false.)  ! <- dagger=.false.
      call set_rule(op_ffg,ttype_op,CLONE_OP,
     &              op_ffg,1,1,
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
     &     occ_def,2,3,(/ 0,0,0,0,0,0/),6)
      call set_rule(op_z_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_z_inter,1,1,
     &              parameters,2,tgt_info)


      if(.not.frozen_old)then
        ! Non-anti-symmetrised Hamiltonian integrals.
        call add_target(op_g_z,ttype_op,.false.,tgt_info)
        occ_def = 0
        ! 1
        occ_def(IHOLE,1,1) = 1
        occ_def(IHOLE,2,1) = 1
        occ_def(IHOLE,1,2) = 1
        occ_def(IHOLE,2,2) = 1
        ! 2
        occ_def(IHOLE,1,3) = 1
        occ_def(IPART,2,3) = 1
        occ_def(IHOLE,1,4) = 1
        occ_def(IHOLE,2,4) = 1
        ! 3
        occ_def(IHOLE,1,5) = 1
        occ_def(IEXTR,2,5) = 1
        occ_def(IHOLE,1,6) = 1
        occ_def(IHOLE,2,6) = 1
        ! 4
        occ_def(IHOLE,1,7) = 1
        occ_def(IHOLE,2,7) = 1
        occ_def(IHOLE,1,8) = 1
        occ_def(IPART,2,8) = 1
        ! 5
        occ_def(IHOLE,1,9) = 1
        occ_def(IPART,2,9) = 1
        occ_def(IHOLE,1,10) = 1
        occ_def(IPART,2,10) = 1
        ! 6
        occ_def(IHOLE,1,11) = 1
        occ_def(IEXTR,2,11) = 1
        occ_def(IHOLE,1,12) = 1
        occ_def(IPART,2,12) = 1
        ! 7
        occ_def(IHOLE,1,13) = 1
        occ_def(IHOLE,2,13) = 1
        occ_def(IPART,1,14) = 1
        occ_def(IHOLE,2,14) = 1
        ! 8
        occ_def(IHOLE,1,15) = 1
        occ_def(IPART,2,15) = 1
        occ_def(IPART,1,16) = 1
        occ_def(IHOLE,2,16) = 1
        ! 9
        occ_def(IHOLE,1,17) = 1
        occ_def(IEXTR,2,17) = 1
        occ_def(IPART,1,18) = 1
        occ_def(IHOLE,2,18) = 1
        ! 10
        occ_def(IHOLE,1,19) = 1
        occ_def(IHOLE,2,19) = 1
        occ_def(IPART,1,20) = 1
        occ_def(IPART,2,20) = 1
        ! 11
        occ_def(IHOLE,1,21) = 1
        occ_def(IPART,2,21) = 1
        occ_def(IPART,1,22) = 1
        occ_def(IPART,2,22) = 1
        ! 12
        occ_def(IHOLE,1,23) = 1
        occ_def(IEXTR,2,23) = 1
        occ_def(IPART,1,24) = 1
        occ_def(IPART,2,24) = 1
        ! 13
        occ_def(IHOLE,1,25) = 1
        occ_def(IHOLE,2,25) = 1
        occ_def(IHOLE,1,26) = 1
        occ_def(IEXTR,2,26) = 1
        ! 14
        occ_def(IHOLE,1,27) = 1
        occ_def(IPART,2,27) = 1
        occ_def(IHOLE,1,28) = 1
        occ_def(IEXTR,2,28) = 1
        ! 15
        occ_def(IHOLE,1,29) = 1
        occ_def(IEXTR,2,29) = 1
        occ_def(IHOLE,1,30) = 1
        occ_def(IEXTR,2,30) = 1
        ! 16
        occ_def(IHOLE,1,31) = 1
        occ_def(IHOLE,2,31) = 1
        occ_def(IPART,1,32) = 1
        occ_def(IEXTR,2,32) = 1
        ! 17
        occ_def(IHOLE,1,33) = 1
        occ_def(IPART,2,33) = 1
        occ_def(IPART,1,34) = 1
        occ_def(IEXTR,2,34) = 1
        ! 18
        occ_def(IHOLE,1,35) = 1
        occ_def(IEXTR,2,35) = 1
        occ_def(IPART,1,36) = 1
        occ_def(IEXTR,2,36) = 1
        ! 19
        occ_def(IHOLE,1,37) = 1
        occ_def(IHOLE,2,37) = 1
        occ_def(IEXTR,1,38) = 1
        occ_def(IHOLE,2,38) = 1
        ! 20
        occ_def(IHOLE,1,39) = 1
        occ_def(IPART,2,39) = 1
        occ_def(IEXTR,1,40) = 1
        occ_def(IHOLE,2,40) = 1
        ! 21
        occ_def(IHOLE,1,41) = 1
        occ_def(IEXTR,2,41) = 1
        occ_def(IEXTR,1,42) = 1
        occ_def(IHOLE,2,42) = 1
        ! 22
        occ_def(IHOLE,1,43) = 1
        occ_def(IHOLE,2,43) = 1
        occ_def(IEXTR,1,44) = 1
        occ_def(IPART,2,44) = 1
        ! 23
        occ_def(IHOLE,1,45) = 1
        occ_def(IPART,2,45) = 1
        occ_def(IEXTR,1,46) = 1
        occ_def(IPART,2,46) = 1
        ! 24
        occ_def(IHOLE,1,47) = 1
        occ_def(IEXTR,2,47) = 1
        occ_def(IEXTR,1,48) = 1
        occ_def(IPART,2,48) = 1
        ! 25
        occ_def(IHOLE,1,49) = 1
        occ_def(IHOLE,2,49) = 1
        occ_def(IEXTR,1,50) = 1
        occ_def(IEXTR,2,50) = 1
        ! 26
        occ_def(IHOLE,1,51) = 1
        occ_def(IPART,2,51) = 1
        occ_def(IEXTR,1,52) = 1
        occ_def(IEXTR,2,52) = 1
        call op_from_occ_parameters(-1,parameters,2,
     &       occ_def,26,2,(/2,2,2,2/),52)
        call set_rule(op_g_z,ttype_op,DEF_OP_FROM_OCC,
     &                op_g_z,1,1,
     &                parameters,2,tgt_info)
      endif

c dbg - test
      call add_target('Z0',ttype_op,.false.,tgt_info)
      ndef = 2
      occ_def = 0
      occ_def(IHOLE,1,2) = 1
      occ_def(IPART,2,2) = 1
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,     0/),6)
      call set_rule('Z0',ttype_op,DEF_OP_FROM_OCC,
     &              'Z0',1,1,
     &              parameters,2,tgt_info)      
c dbg

c dbg
      call add_target(op_z_test,ttype_op,.false.,tgt_info)
      call set_dependency(op_z_test,op_z_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_z_inter,.false.) ! <- dagger=.false.
      call set_rule(op_z_test,ttype_op,CLONE_OP,
     &              op_z_test,1,1,
     &              parameters,1,tgt_info)
c dbg

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
c      ! R-bar intermediate
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_r12_r12bar
c      labels(2) = op_rintbar
c      labels(3) = op_rint
c      labels(4) = op_rintx
c      labels(5) = op_hartree
c      ! for GBC: pass op_exchange instead
c      call add_target(form_r12_r12bar,ttype_frm,.false.,tgt_info)
c      call set_dependency(form_r12_r12bar,op_rintbar,tgt_info)
c      call set_dependency(form_r12_r12bar,op_rint,tgt_info)
c      call set_dependency(form_r12_r12bar,op_rintx,tgt_info)
c      call set_dependency(form_r12_r12bar,op_hartree,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_r12_rbar,ansatz,'RB'//approx)
c      call set_rule(form_r12_r12bar,ttype_frm,DEF_R12INTM_CABS,
c     &              labels,5,1,
c     &              parameters,2,tgt_info)
      
c      ! R-tilde intermediate
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_r12_r12tilde
c      labels(2) = op_rinttilde
c      labels(3) = op_rint
c      labels(4) = op_rintx
c      labels(5) = op_exchange
c      call add_target(form_r12_r12tilde,ttype_frm,.false.,tgt_info)
c      call set_dependency(form_r12_r12tilde,op_rinttilde,tgt_info)
c      call set_dependency(form_r12_r12tilde,op_rint,tgt_info)
c      call set_dependency(form_r12_r12tilde,op_rintx,tgt_info)
c      call set_dependency(form_r12_r12tilde,op_exchange,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_r12_rtilde,ansatz,'RT'//approx)
c      call set_rule(form_r12_r12tilde,ttype_frm,DEF_R12INTM_CABS,
c     &              labels,5,1,
c     &              parameters,2,tgt_info)

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
      ! replace formal G-X operator by H and 
      ! additional G-Acore list (containing core contributions)
      call set_dependency(form_r12_vcabs,op_ham,tgt_info)
      nint = 1
      labels(1) = form_r12_vcabs
      labels(2) = form_r12_vcabs
      labels(3) = op_g_x//'^+'
      labels(4) = op_ham
      if (frozen) then
        call set_dependency(form_r12_vcabs,'G-Acore',tgt_info)
        nint = 2
        labels(5) = op_g_x//'^+'
        labels(6) = 'G-Acore'
c      else if (t1ext.eq.1) then
c        call set_dependency(form_r12_vcabs,'H-ext',tgt_info)
c        nint = 2
c        labels(3) = op_g_x//'^+'
c        labels(4) = 'H-ext'        
c        labels(5) = op_g_x//'^+'
c        labels(6) = op_ham
      end if
      call form_parameters(-1,
     &     parameters,2,title_r12_vcabs,nint,'---')
      call set_rule(form_r12_vcabs,ttype_frm,REPLACE,
     &              labels,2+nint*2,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to V-Ccore
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'V-Ccore-CABS'
      labels(2) = 'V-Ccore'
      labels(3) = 'G-XX'
      labels(4) = op_rint
      labels(5) = op_gr_x   
      call add_target('V-Ccore-CABS',ttype_frm,.false.,tgt_info)
      call set_dependency('V-Ccore-CABS',op_v_inter,tgt_info)
      call set_dependency('V-Ccore-CABS',op_gr,tgt_info)
      call set_dependency('V-Ccore-CABS','G-XX',tgt_info)
      call set_dependency('V-Ccore-CABS',op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'V-Ccore-CABS',ansatz,'V '//approx)
      call set_rule('V-Ccore-CABS',ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      ! replace formal G-X operator by H and 
      ! additional G-Acore list (containing core contributions)
      call set_dependency('V-Ccore-CABS',op_gr_x,tgt_info)
      call set_dependency('V-Ccore-CABS',op_ham,tgt_info)      
      call set_dependency('V-Ccore-CABS',op_ham,tgt_info)
      call set_dependency('V-Ccore-CABS','G-Acore',tgt_info)
      call set_dependency('V-Ccore-CABS','G-CAcore',tgt_info)
      nint = 5
      labels(1) = 'V-Ccore-CABS'
      labels(2) = 'V-Ccore-CABS'
      labels(3) = op_gr_x
      labels(4) = op_gr
      labels(5) = op_gr_x
      labels(6) = 'G.R-Ccore'
      labels(7)  = 'G-XX^+'
      labels(8)  = op_ham
      labels(9)  = 'G-XX^+'
      labels(10) = 'G-Acore'
      labels(11) = 'G-XX^+'
      labels(12) = 'G-CAcore'
      call form_parameters(-1,
     &     parameters,2,title_r12_vcabs,nint,'---')
      call set_rule('V-Ccore-CABS',ttype_frm,REPLACE,
     &              labels,2+nint*2,1,
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
c      ! formal definition of X-test
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_x_test
c      labels(2) = op_x_test
c      labels(3) = op_r12
c      labels(4) = op_r12
c      call add_target(form_x_test,ttype_frm,.false.,tgt_info)
c      call set_dependency(form_x_test,op_x_test,tgt_info)
c      call set_dependency(form_x_test,op_r12,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_x_test,0,'rxr')
c      call set_rule(form_x_test,ttype_frm,DEF_R12INTM_FORMAL,
c     &              labels,4,1,
c     &              parameters,2,tgt_info)
c
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_x_test
c      labels(2) = form_x_test
c      labels(3) = op_r12
c      labels(4) = op_rint
c      call set_dependency(form_x_test,op_rint,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_x_test,1,'---')
c      call set_rule(form_x_test,ttype_frm,REPLACE,
c     &     labels,4,1,
c     &     parameters,2,tgt_info)
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

      if(.not.frozen_old)then
        ! CABS approximation to P
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_r12_pcabs
        labels(2) = op_p_inter
        labels(3) = op_rint
        labels(4) = op_rint  ! dummy, unused
c        labels(5) = op_gr
c        labels(6) = op_v_inter
        labels(5) = op_gr_x
        labels(6) = op_v_x
        labels(7) = op_ffg
        call add_target(form_r12_pcabs,ttype_frm,.false.,tgt_info)
        call set_dependency(form_r12_pcabs,op_p_inter,tgt_info)
        call set_dependency(form_r12_pcabs,op_ffg,tgt_info)
c        call set_dependency(form_r12_pcabs,op_gr,tgt_info)
        call set_dependency(form_r12_pcabs,op_gr_x,tgt_info)
c        call set_dependency(form_r12_pcabs,op_g_x,tgt_info)
        call set_dependency(form_r12_pcabs,op_rint,tgt_info)
c        call set_dependency(form_r12_pcabs,op_v_inter,tgt_info)
        call set_dependency(form_r12_pcabs,op_v_x,tgt_info)
        approx(12:12) = 'S'     ! set symmetrization flag
        call form_parameters(-1,
     &       parameters,2,title_r12_pcabs,ansatz,'P '//approx)
        approx(12:12) = ' '     ! unset flag
        call set_rule(form_r12_pcabs,ttype_frm,DEF_R12INTM_CABS,
     &                labels,7,1,
     &                parameters,2,tgt_info)
        ! replace formal operators by G.R-X and V-X
        call set_dependency(form_r12_pcabs,op_gr,tgt_info)
        call set_dependency(form_r12_pcabs,op_v_inter,tgt_info)
        nint = 2
        labels(1) = form_r12_pcabs
        labels(2) = form_r12_pcabs
        labels(3) = op_gr_x
        labels(4) = op_gr
c     labels(3) = op_gr_x//'^+'
c     labels(4) = op_gr//'^+'
        labels(5) = op_v_x//'^+'
        labels(6) = op_v_inter//'^+'
        if (frozen) then
          nint = 4
          call set_dependency(form_r12_pcabs,'G.R-Ccore',tgt_info)
          call set_dependency(form_r12_pcabs,'V-Ccore',tgt_info)
c     labels(7) = op_gr_x//'^+'
c     labels(8) = 'G.R-Ccore^+'
          labels(7) = op_gr_x
          labels(8) = 'G.R-Ccore'
          labels(9) = op_v_x//'^+'
          labels(10) = 'V-Ccore^+'        
        end if
        call form_parameters(-1,
     &       parameters,2,title_r12_pcabs,nint,'---')
        call set_rule(form_r12_pcabs,ttype_frm,REPLACE,
     &       labels,2+nint*2,1,
     &       parameters,2,tgt_info)

      endif

c dbg
c      ! formal definition of P-test
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_p_test
c      labels(2) = op_p_test
c      labels(3) = op_r12
c      labels(4) = op_ham
c      labels(5) = op_r12
c      call add_target(form_p_test,ttype_frm,.false.,tgt_info)
c      call set_dependency(form_p_test,op_p_test,tgt_info)
c      call set_dependency(form_p_test,op_ham,tgt_info)
c      call set_dependency(form_p_test,op_r12,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_p_test,0,'P')
c      call set_rule(form_p_test,ttype_frm,DEF_R12INTM_FORMAL,
c     &              labels,5,1,
c     &              parameters,2,tgt_info)
c
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_p_test
c      labels(2) = form_p_test
c      labels(3) = op_r12
c      labels(4) = op_rint
c      call set_dependency(form_p_test,op_rint,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_p_test,1,'---')
c      call set_rule(form_p_test,ttype_frm,REPLACE,
c     &     labels,4,1,
c     &     parameters,2,tgt_info)
c
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_p_test
c      labels(2) = form_p_test
c      labels(3) = op_ham
c      labels(4) = op_g_x
c      call set_dependency(form_p_test,op_g_x,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_p_test,1,'---')
c      call set_rule(form_p_test,ttype_frm,REPLACE,
c     &     labels,4,1,
c     &     parameters,2,tgt_info)
c dbg

      ! formal definition of Z
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_zint
      labels(2) = op_z_inter
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      if(truncate.and.(trunc_type.eq.2.or.trunc_type.eq.3))then
        call add_target(form_r12_zint,ttype_frm,.true.,tgt_info)
      else
        call add_target(form_r12_zint,ttype_frm,.false.,tgt_info)
      endif
      call set_dependency(form_r12_zint,op_z_inter,tgt_info)
      call set_dependency(form_r12_zint,op_r12,tgt_info)
      call set_dependency(form_r12_zint,op_ham,tgt_info)
      call set_dependency(form_r12_zint,op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_zint,0,'Z')
      call set_rule(form_r12_zint,ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)

      if(.not.frozen_old)then
        ! CABS approximation to Z.
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_r12_zcabs
        labels(2) = op_z_inter
        labels(3) = op_rint
        labels(4) = op_g_z
        labels(5) = 'R.R-X'
        call add_target(form_r12_zcabs,ttype_frm,.false.,tgt_info)
        call set_dependency(form_r12_zcabs,op_z_inter,tgt_info)
        call set_dependency(form_r12_zcabs,'R.R-X',tgt_info)
        call set_dependency(form_r12_zcabs,op_rint,tgt_info)
        call set_dependency(form_r12_zcabs,op_g_z,tgt_info)
c        approx(12:12) = 'S' ! set symmetrization flag
        call form_parameters(-1,
     &       parameters,2,title_r12_zcabs,ansatz,'Z '//approx)
c       approx(12:12) = ' ' ! unset flag
        call set_rule(form_r12_zcabs,ttype_frm,DEF_R12INTM_CABS,
     &                labels,5,1,
     &                parameters,2,tgt_info)
      endif

c dbg - test
      call add_target('Z0TEST-FRM',ttype_frm,.false.,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z0TEST-FRM'
      labels(2) = 'Z0'
      labels(3) = 'Z0'
      labels(4) = op_z_inter
      labels(5) = op_z_inter
      labels(6) = op_z_inter
      labels(7) = 'Z0'
      idx_sv(1:5) = (/1,2,2,2,1/)
      iblkmin(1:5) = -1
      iblkmax(1:5) = -1
      nconnect = 0
      connect  = 0
      navoid = 0
      ninproj = 0
      call set_dependency('Z0TEST-FRM','Z0',tgt_info)
      call set_dependency('Z0TEST-FRM',op_z_inter,tgt_info)
      call expand_parameters(-1,
     &     parameters,3,
     &     'XXX',5,idx_sv,iblkmin,iblkmax,
     &     connect,nconnect,
     &     0,navoid,
     &     0,ninproj)
      call set_rule('Z0TEST-FRM',ttype_frm,EXPAND_OP_PRODUCT,
     &     labels,7,1,
     &     parameters,3,tgt_info)
c dbg

c dbg
c      ! formal definition of Z-test
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_z_test
c      labels(2) = op_z_test
c      labels(3) = op_r12
cc      labels(4) = op_ham ! Full
c      labels(4) = op_g_z ! Coulomb or exchange
c      labels(5) = op_r12
c      call add_target(form_z_test,ttype_frm,.false.,tgt_info)
c      call set_dependency(form_z_test,op_z_test,tgt_info)
cc      call set_dependency(form_z_test,op_ham,tgt_info) ! Full
c      call set_dependency(form_z_test,op_g_z,tgt_info) ! Coulomb or exchnage
c      call set_dependency(form_z_test,op_r12,tgt_info)
c      call form_parameters(-1,
cc     &     parameters,2,title_z_test,0,'Z') ! Full
c     &     parameters,2,title_z_test,0,'ZT') ! Coulomb or exchange.
c      call set_rule(form_z_test,ttype_frm,DEF_R12INTM_FORMAL,
c     &              labels,5,1,
c     &              parameters,2,tgt_info)
c
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_z_test
c      labels(2) = form_z_test
c      labels(3) = op_r12
c      labels(4) = op_rint
c      call set_dependency(form_z_test,op_rint,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_z_test,1,'---')
c      call set_rule(form_z_test,ttype_frm,REPLACE,
c     &     labels,4,1,
c     &     parameters,2,tgt_info)

      ! Only if full.
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_z_test
c      labels(2) = form_z_test
c      labels(3) = op_ham
c      labels(4) = op_g_x
c      call set_dependency(form_z_test,op_g_x,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_z_test,1,'---')
c      call set_rule(form_z_test,ttype_frm,REPLACE,
c     &     labels,4,1,
c     &     parameters,2,tgt_info)
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
c      call set_dependency(fopt_r12_vcabs,mel_gintx,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_ham,tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_rint,tgt_info)      
      if (frozen)
     &     call set_dependency(fopt_r12_vcabs,'G-Ac-INT',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_r12_vcabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! set V-Ccore
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'V-Ccore-OPT'
      labels(2) = 'V-Ccore-CABS'
      ncat = 1
      nint = 0
      call add_target('V-Ccore-OPT',ttype_frm,.false.,tgt_info)
      call set_dependency('V-Ccore-OPT','V-Ccore-CABS',tgt_info)
      call set_dependency('V-Ccore-OPT','G-Ac-INT',tgt_info)
      call set_dependency('V-Ccore-OPT','G-CAc-INT',tgt_info)
      call set_dependency('V-Ccore-OPT','G.R-Cc-INT',tgt_info)
      call set_dependency('V-Ccore-OPT','DEF-V-Cc-INTER',tgt_info)
      call set_dependency('V-Ccore-OPT',mel_rint,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('V-Ccore-OPT',ttype_frm,OPTIMIZE,
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

      if(.not.frozen_old)then
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
c        call set_dependency(fopt_r12_pcabs,mel_gintx,tgt_info)      
        call set_dependency(fopt_r12_pcabs,mel_ffg,tgt_info)      
        call set_dependency(fopt_r12_pcabs,mel_gr,tgt_info)      
        if (frozen) then
          call set_dependency(fopt_r12_pcabs,'G.R-Cc-INT',tgt_info)      
          call set_dependency(fopt_r12_pcabs,'DEF-V-Cc-INTER',tgt_info)      
        end if
        call opt_parameters(-1,parameters,ncat,nint)
        call set_rule(fopt_r12_pcabs,ttype_frm,OPTIMIZE,
     &                labels,ncat+nint+1,1,
     &                parameters,1,tgt_info)
      endif


      if(.not.frozen_old)then
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
        call set_dependency(fopt_r12_zcabs,mel_gintz,tgt_info)      
        call set_dependency(fopt_r12_zcabs,'R.R-INTX',tgt_info)      
        call opt_parameters(-1,parameters,ncat,nint)
        call set_rule(fopt_r12_zcabs,ttype_frm,OPTIMIZE,
     &                labels,ncat+nint+1,1,
     &                parameters,1,tgt_info)

      endif        

      ! set Z (direct eval)
      call add_target('ZINT_R12_DIR',ttype_frm,.false.,tgt_info)
      call set_dependency('ZINT_R12_DIR',form_r12_zint,tgt_info)
      call set_dependency('ZINT_R12_DIR',mel_z_def,tgt_info)
      call set_dependency('ZINT_R12_DIR',mel_rint,tgt_info)
      call set_dependency('ZINT_R12_DIR',mel_ham,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ZINT_R12_DIR'
      labels(2) = form_r12_zint
      labels(3) = op_r12
      labels(4) = op_rint
      labels(5) = op_r12//'^+'
      labels(6) = op_rint//'^+'
      nint = 2
      call form_parameters(-1,
     &     parameters,2,'Z direct',nint,'---')
      call set_rule('ZINT_R12_DIR',ttype_frm,REPLACE,
     &     labels,2+nint*2,1,
     &     parameters,2,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ZINT_R12_DIR'
      labels(2) = 'ZINT_R12_DIR'
      ncat = 1
      nint = 0
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('ZINT_R12_DIR',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)


c dbg - test
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z0TEST-OPT'
      labels(2) = 'Z0TEST-FRM'
      ncat = 1
      nint = 0
      call add_target('Z0TEST-OPT',ttype_frm,.false.,tgt_info)
      call set_dependency('Z0TEST-OPT','Z0TEST-FRM',tgt_info)
      call set_dependency('Z0TEST-OPT',mel_z_def,tgt_info)
      call set_dependency('Z0TEST-OPT','DEF-Z0TEST',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('Z0TEST-OPT',ttype_frm,OPTIMIZE,
     &                labels,ncat+nint+1,1,
     &                parameters,1,tgt_info)
c dbg 

      ! set C intermediate
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

c      ! R12integrals (extended list)
c      call add_target(mel_rintx,ttype_opme,.false.,tgt_info)
c      call set_dependency(mel_rintx,op_rintx,tgt_info)
c      ! (a) define
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_rintx
c      labels(2) = op_rintx
c      call me_list_parameters(-1,parameters,
c     &     msc,0,1,0,0,.false.)
c      call set_rule(mel_rintx,ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)
c      ! (b) import
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_rintx
c      call import_parameters(-1,parameters,'F12_INT',env_type)
c      call set_rule(mel_rintx,ttype_opme,IMPORT,
c     &              labels,1,1,
c     &              parameters,1,tgt_info)

      ! special one-electron integral list
      call add_target('F-X-INT',ttype_opme,.false.,tgt_info)
      call set_dependency('F-X-INT','F-X',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'F-X-INT'
      labels(2) = 'F-X'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('F-X-INT',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'F-X-INT'
      call import_parameters(-1,parameters,'F_INT',env_type)
      call set_rule('F-X-INT',ttype_opme,IMPORT,
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

      ! new special two-electron integral list
      call add_target('G-Ac-INT',ttype_opme,.false.,tgt_info)
      call set_dependency('G-Ac-INT','G-Acore',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'G-Ac-INT'
      labels(2) = 'G-Acore'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('G-Ac-INT',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'G-Ac-INT'
      call import_parameters(-1,parameters,'G_INT',env_type)
      call set_rule('G-Ac-INT',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! new special two-electron integral list (extended list with HX|HX block)
      call add_target('H-ext-INT',ttype_opme,.false.,tgt_info)
      call set_dependency('H-ext-INT','H-ext',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'H-ext-INT'
      labels(2) = 'H-ext'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('H-ext-INT',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'H-ext-INT'
      call import_parameters(-1,parameters,'G_INT',env_type)
      call set_rule('H-ext-INT',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! new additional two-electron integral list
      call add_target('G-CAc-INT',ttype_opme,.false.,tgt_info)
      call set_dependency('G-CAc-INT','G-CAcore',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'G-CAc-INT'
      labels(2) = 'G-CAcore'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('G-CAc-INT',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'G-CAc-INT'
      call import_parameters(-1,parameters,'G_INT',env_type)
      call set_rule('G-CAc-INT',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      if(.not.frozen_old)then
        ! special two-electron integral list 2
        call add_target(mel_gintz,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_gintz,op_g_z,tgt_info)
        ! (a) define
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_gintz
        labels(2) = op_g_z
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.true.)
c     &       0,0,1,0,0,.false.)
        call set_rule(mel_gintz,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        ! (b) import
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_gintz
        call import_parameters(-1,parameters,'G_INT',env_type)
        call set_rule(mel_gintz,ttype_opme,IMPORT,
     &                labels,1,1,
     &                parameters,1,tgt_info)
      endif

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

      ! R12^2 integrals II
      call add_target('R.R-INTX',ttype_opme,.false.,tgt_info)
      call set_dependency('R.R-INTX','R.R-X',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'R.R-INTX'
      labels(2) = 'R.R-X'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('R.R-INTX',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'R.R-INTX'
      call import_parameters(-1,parameters,'FF_INT',env_type)
      call set_rule('R.R-INTX',ttype_opme,IMPORT,
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

      call add_target('G.R-Cc-INT',ttype_opme,.false.,tgt_info)
      call set_dependency('G.R-Cc-INT','G.R-Ccore',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'G.R-Cc-INT'
      labels(2) = 'G.R-Ccore'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('G.R-Cc-INT',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'G.R-Cc-INT'
      call import_parameters(-1,parameters,'FG_INT',env_type)
      call set_rule('G.R-Cc-INT',ttype_opme,IMPORT,
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_ffg,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ffg
      call import_parameters(-1,parameters,'FFG_INT',env_type)
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
c      ! OLD:
c      ! (b) import
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_rintbreve
c      call import_parameters(-1,parameters,'F12BRV_INT',env_type)
c      call set_rule(mel_rintbreve,ttype_opme,IMPORT,
c     &              labels,1,1,
c     &              parameters,1,tgt_info)
      ! (b) define formula
      call set_dependency(mel_rintbreve,'F-X',tgt_info)
      call set_dependency(mel_rintbreve,'F-X-INT',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'BBRV_FRM'
      labels(2) = op_rintbreve
      labels(3) = op_rint
      labels(4) = '-'
      labels(5) = 'F-X'
      call form_parameters(-1,
     &     parameters,2,'R12BREVE',3,'RV')
      call set_rule(mel_rintbreve,ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule(mel_rintbreve,ttype_frm,PRINT_FORMULA,
     &              labels,1,0,
     &              parameters,2,tgt_info)
      ! (c) optimize formula
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'BBRV_OPT'
      labels(2) = 'BBRV_FRM'
      call opt_parameters(-1,parameters,1,0)
      call set_rule(mel_rintbreve,ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (d) evaluate
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'BBRV_OPT'
      call set_rule(mel_rintbreve,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      ! Exchange integrals, K.
      call add_target(mel_exchange,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_exchange,op_exchange,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_exchange
      labels(2) = op_exchange
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_exchange,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_exchange
      call import_parameters(-1,parameters,'K_INT',env_type)
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
     &     msc,0,1,0,0,.false.)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_v_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! V-Ccore-list
      call add_target('DEF-V-Cc-INTER',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF-V-Cc-INTER','V-Ccore',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'V-Cc-INTER'
      labels(2) = 'V-Ccore'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF-V-Cc-INTER',ttype_opme,DEF_ME_LIST,
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

c dbg
c      ! X-test
c      call add_target(mel_x_test_def,ttype_opme,.false.,tgt_info)
c      call set_dependency(mel_x_test_def,op_x_test,tgt_info)
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_x_test
c      labels(2) = op_x_test
c      call me_list_parameters(-1,parameters,
c     &     msc,0,1,0,0,.false.)
c      call set_rule(mel_x_test_def,ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)
c dbg

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

      ! P-list
      if(frozen_old)then
        call add_target(mel_p_def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_p_def,op_p_inter,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_p_def
        labels(2) = op_p_inter
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_p_def,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
        ! import
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_p_def
        call import_parameters(-1,parameters,'P_LIST',env_type)
        call set_rule(mel_p_def,ttype_opme,IMPORT,
     &                labels,1,1,
     &                parameters,1,tgt_info)
      else
        call add_target(mel_p_def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_p_def,op_p_inter,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_p_inter
        labels(2) = op_p_inter
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_p_def,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
      endif

c dbg
c      ! P-test
c      call add_target(mel_p_test_def,ttype_opme,.false.,tgt_info)
c      call set_dependency(mel_p_test_def,op_p_test,tgt_info)
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_p_test
c      labels(2) = op_p_test
c      call me_list_parameters(-1,parameters,
c     &     msc,0,1,0,0,.false.)
c      call set_rule(mel_p_test_def,ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)
c dbg

      ! Z-list
      if(frozen_old)then
        call add_target(mel_z_def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_z_def,op_z_inter,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_z_def
        labels(2) = op_z_inter
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_z_def,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
        ! import
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_z_def
        call import_parameters(-1,parameters,'Z_LIST',env_type)
        call set_rule(mel_z_def,ttype_opme,IMPORT,
     &                labels,1,1,
     &                parameters,1,tgt_info)
      else
        call add_target(mel_z_def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_z_def,op_z_inter,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_z_inter
        labels(2) = op_z_inter
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_z_def,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
      endif

c dbg - test
      call add_target('DEF-Z0TEST',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF-Z0TEST','Z0',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z0TEST'
      labels(2) = 'Z0'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF-Z0TEST',ttype_opme,DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)
c dbg

c dbg
      ! Z-test
      call add_target(mel_z_test_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_z_test_def,op_z_test,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_z_test
      labels(2) = op_z_test
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
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
     &     msc,0,1,0,0,.false.)
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
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_b_inv,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1) = mel_b_inv   ! output
      labels(2) = mel_b_inter ! input
      call form_parameters(-1,parameters,2,
     &     '---',0,'inverse')
      call set_rule(mel_b_inv,ttype_opme,INVERT,
     &              labels,2,1,
     &              parameters,2,tgt_info)

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
     &       msc,0,1,0,0,.false.)
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
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_x_inv,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        labels(1) = mel_x_inv   ! output
        labels(2) = mel_x_inter ! input
      call form_parameters(-1,parameters,2,
     &     '---',0,'inverse')
        call set_rule(mel_x_inv,ttype_opme,INVERT,
     &                labels,2,1,
     &                parameters,2,tgt_info)
      
      endif

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      ! test
c      if(.not.pz_eval)then
        call add_target(eval_r12_inter,ttype_gen,.false.,tgt_info)
c      else
c        call add_target(eval_r12_inter,ttype_gen,.true.,tgt_info)
c      endif
      call set_dependency(eval_r12_inter,mel_ham,tgt_info)
      call set_dependency(eval_r12_inter,mel_rint,tgt_info)
      call set_dependency(eval_r12_inter,mel_gintx,tgt_info)
      call set_dependency(eval_r12_inter,mel_ttr,tgt_info)
      call set_dependency(eval_r12_inter,mel_ff,tgt_info)
      call set_dependency(eval_r12_inter,mel_gr,tgt_info)
c      call set_dependency(eval_r12_inter,mel_v_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_x_def,tgt_info)
      call set_dependency(eval_r12_inter,mel_b_def,tgt_info)
      call set_dependency(eval_r12_inter,'EVAL_V',tgt_info)
      if(is_keyword_set('method.CC').gt.0.and.(.not.truncate
     &     .or.(truncate.and.trunc_type.gt.0)))then
        call set_dependency(eval_r12_inter,'EVAL_P',tgt_info)
        call set_dependency(eval_r12_inter,'EVAL_Z',tgt_info)
c        call set_dependency(eval_r12_inter,mel_ffg,tgt_info)
c        if(.not.frozen_old)then
c          call set_dependency(eval_r12_inter,mel_p_def,tgt_info)
c          call set_dependency(eval_r12_inter,mel_z_def,tgt_info)
c          call set_dependency(eval_r12_inter,fopt_r12_pcabs,tgt_info)
c          call set_dependency(eval_r12_inter,fopt_r12_zcabs,tgt_info)
c          call set_dependency(eval_r12_inter,mel_gintz,tgt_info)
c        endif
      endif
c      call set_dependency(eval_r12_inter,fopt_r12_vcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_xcabs,tgt_info)
      call set_dependency(eval_r12_inter,fopt_r12_bcabs,tgt_info)
      if (ansatz.ne.1)then
        call set_dependency(eval_r12_inter,fopt_r12_ccabs,tgt_info)
        call set_dependency(eval_r12_inter,mel_c_def,tgt_info)
      endif
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = fopt_r12_vcabs
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

c      if(is_keyword_set('method.CC').gt.0.and.(.not.truncate
c     &     .or.(truncate.and.trunc_type.gt.0)))then
c        if(.not.frozen_old)then
c          labels(1) = fopt_r12_pcabs
c          call set_rule(eval_r12_inter,ttype_opme,EVAL,
c     &         labels,1,0,
c     &         parameters,0,tgt_info)
c          
c          call set_dependency(eval_r12_inter,'R.R-INTX',tgt_info)
c          labels(1) = fopt_r12_zcabs
c          call set_rule(eval_r12_inter,ttype_opme,EVAL,
c     &         labels,1,0,
c     &         parameters,0,tgt_info)
c
c        endif
c      endif

      call add_target('EVAL_V',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_V',fopt_r12_vcabs,tgt_info)
c      call set_dependency('EVAL_V',mel_ham,tgt_info)
c      call set_dependency('EVAL_V',mel_v_def,tgt_info)
c      call set_dependency('EVAL_V',mel_rint,tgt_info)
c      call set_dependency('EVAL_V',mel_gr,tgt_info)
      labels(1) = fopt_r12_vcabs
      call set_rule('EVAL_V',ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)

c      call add_target('EVAL_P',ttype_gen,pz_eval,tgt_info)
      call add_target('EVAL_P',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_P',fopt_r12_pcabs,tgt_info)
      call set_dependency('EVAL_P',mel_p_def,tgt_info)
      call set_dependency('EVAL_P',mel_rint,tgt_info)
      call set_dependency('EVAL_P',mel_ffg,tgt_info)
      call set_dependency('EVAL_P',mel_gr,tgt_info)
c      call set_dependency('EVAL_P',eval_r12_inter,tgt_info)
      if (frozen) then
        call set_dependency('EVAL_P','V-Ccore-OPT',tgt_info)

        labels(1) = 'V-Ccore-OPT'
        call set_rule('EVAL_P',ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)
      end if
      labels(1) = fopt_r12_pcabs
      call set_rule('EVAL_P',ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)

c      call add_target('EVAL_Z',ttype_gen,pz_eval,tgt_info)
      call add_target('EVAL_Z',ttype_gen,pz_eval,tgt_info)
      call set_dependency('EVAL_Z',mel_z_def,tgt_info)
      call set_dependency('EVAL_Z',mel_rint,tgt_info)
      call set_dependency('EVAL_Z','R.R-INTX',tgt_info)
      if (approx(14:17).eq.'DRCT') then
        call set_dependency('EVAL_Z',mel_ham,tgt_info)
        call set_dependency('EVAL_Z','ZINT_R12_DIR',tgt_info)
      else
        call set_dependency('EVAL_Z',mel_gintz,tgt_info)
        call set_dependency('EVAL_Z',fopt_r12_zcabs,tgt_info)
      end if
      if (approx(14:17).eq.'DRCT') then
        labels(1) = 'ZINT_R12_DIR'
      else
        labels(1) = fopt_r12_zcabs
      end if
      call set_rule('EVAL_Z',ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)
c dbg -- test
c      call set_dependency('EVAL_Z','Z0TEST-OPT',tgt_info)
c      labels(1) = 'Z0TEST-OPT'
c      call set_rule('EVAL_Z',ttype_opme,EVAL,
c     &         labels,1,0,
c     &         parameters,0,tgt_info)
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

      subroutine set_gxx(occ_def)

      implicit none

      integer, intent(out) ::
     &     occ_def(ngastp,2,*)

      integer ::
     &     idxc, idxa, occ_c(ngastp), occ_a(ngastp)

      occ_def(1:ngastp,1:2,1:4*4) = 0
      ndef = 0
      do idxc = 1, 4
        occ_c = 0
        if (idxc.eq.1) occ_c(IHOLE) = 2
        if (idxc.eq.2) occ_c(IHOLE) = 1
        if (idxc.eq.2) occ_c(IPART) = 1
        if (idxc.eq.3) occ_c(IHOLE) = 1
        if (idxc.eq.3) occ_c(IEXTR) = 1
        if (idxc.eq.4) occ_c(IPART) = 2
        do idxa = 1, 4
          occ_a = 0
          if (idxa.eq.1) occ_a(IHOLE) = 2
          if (idxa.eq.2) occ_a(IHOLE) = 1
          if (idxa.eq.2) occ_a(IPART) = 1
          if (idxa.eq.3) occ_a(IHOLE) = 1
          if (idxa.eq.3) occ_a(IEXTR) = 1
          if (idxa.eq.4) occ_a(IPART) = 2

          ndef = ndef+1
          occ_def(1:ngastp,1,ndef) = occ_c(1:ngastp)
          occ_def(1:ngastp,2,ndef) = occ_a(1:ngastp)
          
        end do
      end do

      end subroutine set_gxx

      end
