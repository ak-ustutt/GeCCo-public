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
      include 'ifc_targets.h'

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
     &     min_rank_tp, min_rank_tpp, nblk, nj,
     &     isim, ncat, nint, icnt, nlab, irank, idef,
     &     isym, ms, msc, sym_arr(8), extend, r12op,
     &     occ_def(ngastp,2,60), vring_mode,
     &     ntp_min, ntp_max, ntpp_min, ntpp_max, t1ext, trunc_type
      logical ::
     &     needed, r12fix, set_tp, set_tpp, truncate, set_RT2T2, CC,
     &     pf12_trunc, frozen, pz_eval, use_CS, xsp_opt1, active_orbs,
     &     semi_r12
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character(20) ::
     &     approx, approx2, F_appr, K_appr, Z_appr, Z2_appr, shell_typ
      character(len=256) ::
     &     descr

      character(*), intent(in) ::
     &     env_type

*----------------------------------------------------------------------*
      if (iprlvl.gt.0)
     &     write(luout,*) 'setting general targets for F12 (SP) ...'

      msc = +1  ! assuming closed shell
*----------------------------------------------------------------------*
*     read input
*----------------------------------------------------------------------*
      ! set approx string
      approx(1:20) = ' '
      approx2(1:20) = ' '
      F_appr(1:20) = ' '
      K_appr(1:20) = ' '
      Z_appr(1:20) = ' '
      Z2_appr(1:20) = ' '
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','F_appr',str=F_appr)
      call get_argument_value('method.R12','K_appr',str=K_appr)
      call get_argument_value('method.R12','Z_appr',str=Z_appr)
      call get_argument_value('method.R12','Z2_appr',str=Z2_appr)
      call get_argument_value('method.R12','min_tp',ival=min_rank_tp)
      call get_argument_value('method.R12','min_tpp',ival=min_rank_tpp)
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','pz_eval',lval=pz_eval)
      call get_argument_value('method.R12','extend',ival=extend)
      call get_argument_value('method.R12','r12op',ival=r12op)
      call get_argument_value('method.R12','T1ext',ival=t1ext)
      call get_argument_value('method.R12','trunc',ival=trunc_type)
      call get_argument_value('method.R12','vring',ival=vring_mode)
      call get_argument_value('method.R12','use_CS',lval=use_CS)
      call get_argument_value('method.R12','xsp1',lval=xsp_opt1)
      call get_argument_value('method.R12','semi_r12',lval=semi_r12)
      truncate = trunc_type.ge.0
      if (is_keyword_set('method.truncate').gt.0) then
        truncate = is_keyword_set('method.truncate').gt.0
        call get_argument_value('method.truncate','trunc_type',
     &       ival=trunc_type)
      end if
      CC = is_keyword_set('method.CC').gt.0.or.
     &     is_keyword_set('method.MRCC').gt.0
      ! new defaults for MRCC:
      if (is_keyword_set('method.MRCC').gt.0) then
        if (iprlvl.gt.0) 
     &      write(luout,*) 'MRCC-F12: switching on (F12*) approx.'
        vring_mode=3
        use_CS=.true.
      end if

      pf12_trunc = truncate.and.trunc_type.eq.0
      ! Frozen core?
      frozen = .false.
      shell_typ(1:) = ' '
      if(is_keyword_set('orb_space.shell').gt.0)then
        call get_argument_value('orb_space.shell','type',
     &       str=shell_typ)
        frozen = trim(shell_typ).eq.'frozen'
      endif

      call get_argument_value('method.CCPT','RT2T2',
     &     lval=set_RT2T2)

      ! active orbitals present?
      active_orbs = orb_info%nactt_hpv(IVALE).gt.0

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
        if (set_RT2T2) n_pp=2
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

      approx2 = approx
      select case(trim(Z2_appr))
      case('as-Z')
        ! do nothing, same as for Z
      case('direct')
        write(luout,*) 'direct RI evaluation of Z2 intermediate'
        approx2(14:17) = 'DRCT'
      case('none','J2K3')
        write(luout,*) 'no approximations to Z2 intermediate made'
        approx2(14:17) = 'J2K3'
      case default
        if (Z2_appr(1:1).ne.'J'.or.Z2_appr(3:3).ne.'K'.or.
     &      (Z2_appr(2:2).ne.'0'.and.
     &       Z2_appr(2:2).ne.'1'.and.
     &       Z2_appr(2:2).ne.'2').or. 
     &      (Z2_appr(4:4).ne.'0'.and.
     &       Z2_appr(4:4).ne.'1'.and.
     &       Z2_appr(4:4).ne.'2'.and.
     &       Z2_appr(4:4).ne.'3')) then
          call quit(0,'set_r12_general_targets',
     &       'Z2_appr unknown: "'//trim(Z2_appr)//'"')
        end if
        write(luout,*) 'approximation to Z2 intermediate: ',
     &      trim(Z2_appr)
        approx2(14:17) = Z2_appr(1:4)
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
c     &       occ_def,ndef,2,(/  0, 0,  0, 0/),ndef)
c        call set_rule(op_cex,ttype_op,DEF_OP_FROM_OCC,
c     &       op_cex,1,1,
c     &       parameters,2,tgt_info)
        if (.not.xsp_opt1) then
          ! the usual route
          call xop_parameters(-1,parameters,
     &       .false.,ntp_min,ntp_max,0,ntp_max+2)
          call set_rule(op_cex,ttype_op,DEF_EXCITATION,
     &                op_cex,1,1,
     &                parameters,1,tgt_info)
        else
          ! define formal op_cex for Lagrange-build
          call xop_parameters(-1,parameters,
     &       .false.,ntp_min,ntp_max,0,ntp_max+2)
          call set_rule(op_cex,ttype_op,DEF_EXCITATION,
     &                'T12FML',1,1,
     &                parameters,1,tgt_info)
          ! and the actual one with only ntp_max = 1
          call xop_parameters(-1,parameters,
     &       .false.,ntp_min,1,0,ntp_max+2)
          call set_rule(op_cex,ttype_op,DEF_EXCITATION,
     &                op_cex,1,1,
     &                parameters,1,tgt_info)
        end if

        ! The Lagrangian multipliers.
        call add_target(op_cexbar,ttype_op,.false.,tgt_info)
        call set_dependency(op_cexbar,op_cex,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_cex,.true.) ! <- dagger=.true.
        call set_rule(op_cexbar,ttype_op,CLONE_OP,
     &                op_cexbar,1,1,
     &                parameters,1,tgt_info)

        if (xsp_opt1) then
          call cloneop_parameters(-1,parameters,
     &                          'T12FML',.true.) ! <- dagger=.true.
          call set_rule(op_cexbar,ttype_op,CLONE_OP,
     &                'T12FBAR',1,1,
     &                parameters,1,tgt_info)
        end if

        call add_target(op_diar12,ttype_op,.false.,tgt_info)
        call set_dependency(op_diar12,op_cex,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_cex,.false.) ! <- dagger=.false.
        call set_rule(op_diar12,ttype_op,CLONE_OP,
     &                op_diar12,1,1,
     &                parameters,1,tgt_info)

      endif

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
      
      ! (pq)_frozen/(pq)_ae block of 2e-Hamiltonian (formal)
      call add_target(op_g_x,ttype_op,.false.,tgt_info)
cc      min_rank = 2 
c      call r12int_parameters(-1,parameters,
c     &     2,2,2,0,2)
c      call set_rule(op_g_x,ttype_op,DEF_R12INT,
c     &              op_g_x,1,1,
c     &              parameters,1,tgt_info)
      ndef = 16
      if (active_orbs) ndef = 64
      call set_gxx(occ_def)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/      2,     0/),ndef)
      call set_rule(op_g_x,ttype_op,DEF_OP_FROM_OCC,
     &              op_g_x,1,1,
     &              parameters,2,tgt_info)



      ! ae/ae blocks of 2e-Hamilt. (formal)
      call add_target('G-XX',ttype_op,.false.,tgt_info)
      ndef = 16
      if (active_orbs) ndef = 64
      call set_gxx(occ_def)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/      2,      2/),ndef)
      call set_rule('G-XX',ttype_op,DEF_OP_FROM_OCC,
     &              'G-XX',1,1,
     &              parameters,2,tgt_info)

      ! (pq)_frozen/(iq)_ae block of 2e-Hamiltonian (actual)
c      call add_target('G-Acore',ttype_op,.false.,tgt_info)
      call add_target2('G-Acore',.false.,tgt_info)
c      min_rank = 2 
      descr = '[HP][HP],H[HPX]'
      if (active_orbs) descr = '[HPV][HPV],H[HPVX]'

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
c      call op_from_occ_parameters(-1,parameters,2,
c     &     occ_def,ndef,1,(/     0,      2/),ndef)
c      call set_rule('G-Acore',ttype_op,DEF_OP_FROM_OCC,
c     &              'G-Acore',1,1,
c     &              parameters,2,tgt_info)
      call set_rule2('G-Acore',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('G-Acore',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'G-Acore'/))
      call set_arg('G-Acore',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('G-Acore',DEF_OP_FROM_OCC,'CORE',2,tgt_info,
     &             val_int=(/0,2/))
      call set_arg('G-Acore',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)

      ! auxiliary HX|HX list
      !call add_target('H-ext',ttype_op,.false.,tgt_info)
      call add_target2('H-ext',.false.,tgt_info)
      descr = 'HX,HX'
      if (active_orbs) descr = '[HV]X,[HV]X'
c      min_rank = 2 
c      occ_def = 0
c      ndef = 1
c      ! 1
c      occ_def(IHOLE,1,1) = 1
c      occ_def(IEXTR,1,1) = 1
c      occ_def(IHOLE,2,1) = 1
c      occ_def(IEXTR,2,1) = 1
c      call op_from_occ_parameters(-1,parameters,2,
c     &     occ_def,ndef,1,(/     0,      2/),ndef)
c      call set_rule('H-ext',ttype_op,DEF_OP_FROM_OCC,
c     &              'H-ext',1,1,
c     &              parameters,2,tgt_info)
      call set_rule2('H-ext',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('H-ext',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'H-ext'/))
      call set_arg('H-ext',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('H-ext',DEF_OP_FROM_OCC,'CORE',2,tgt_info,
     &             val_int=(/0,2/))
      call set_arg('H-ext',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)


      ! (jq)_ae/(ir)_ae block of 2e-Hamiltonian
      !call add_target('G-CAcore',ttype_op,.false.,tgt_info)
      call add_target2('G-CAcore',.false.,tgt_info)
      descr = 'H[HPX],H[HPX]|H[HPX],PP'

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
c      call op_from_occ_parameters(-1,parameters,2,
c     &     occ_def,ndef,1,(/      2,      2/),ndef)
c      call set_rule('G-CAcore',ttype_op,DEF_OP_FROM_OCC,
c     &              'G-CAcore',1,1,
c     &              parameters,2,tgt_info)
      call set_rule2('G-CAcore',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('G-CAcore',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'G-CAcore'/))
      call set_arg('G-CAcore',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('G-CAcore',DEF_OP_FROM_OCC,'CORE',2,tgt_info,
     &             val_int=(/2,2/))
      call set_arg('G-CAcore',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)

      ! i,(a/x) block of Fock
      !call add_target('F-X',ttype_op,.false.,tgt_info)
      call add_target2('F-X',.false.,tgt_info)
      descr = 'H,[PX]'
      if (active_orbs) descr = '[HV],[PX]'
      occ_def = 0
      ndef = 2
      occ_def(IHOLE,1,1) = 1
      occ_def(IPART,2,1) = 1
      occ_def(IHOLE,1,2) = 1
      occ_def(IEXTR,2,2) = 1
c      call op_from_occ_parameters(-1,parameters,2,
c     &     occ_def,ndef,1,(/      2,     0/),ndef)
c      call set_rule('F-X',ttype_op,DEF_OP_FROM_OCC,
c     &              'F-X',1,1,
c     &              parameters,2,tgt_info)
      call set_rule2('F-X',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('F-X',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'F-X'/))
      call set_arg('F-X',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('F-X',DEF_OP_FROM_OCC,'CORE',2,tgt_info,
     &             val_int=(/2,0/))
      call set_arg('F-X',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)

      ! for MR: averaged Fock:
      call add_target2('Favg',.false.,tgt_info)
      call set_rule2('Favg',DEF_HAMILTONIAN,tgt_info)
      call set_arg('Favg',DEF_HAMILTONIAN,'LABEL',1,tgt_info,
     &     val_label=(/'Favg'/))
      call set_arg('Favg',DEF_HAMILTONIAN,'MIN_RANK',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('Favg',DEF_HAMILTONIAN,'MAX_RANK',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('Favg',DEF_HAMILTONIAN,'SET_X',1,tgt_info,
     &     val_log=(/.true./))


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
      !call add_target(op_rintbreve,ttype_op,.false.,tgt_info)
      call add_target2(op_rintbreve,.false.,tgt_info)
c      call r12int_parameters(-1,parameters,
c     &     n_pp,min_rank,2,0,2)
c      call set_rule(op_rintbreve,ttype_op,DEF_R12INT,
c     &              op_rintbreve,1,1,
c     &              parameters,1,tgt_info)
      occ_def = 0
      ndef = 2
      descr = 'H[PX],HH'
      if (active_orbs) descr = '[HV][PX],[HV][HV]'
      ! 1
      occ_def(IHOLE,1,1) = 1
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,1) = 2
      ! 2
      occ_def(IHOLE,1,2) = 1
      occ_def(IEXTR,1,2) = 1
      occ_def(IHOLE,2,2) = 2
      if (n_pp.ge.1) then
        descr = 'H[PX],H[HP]'
        if (active_orbs) descr = '[HV][PX],[HV][HVP]'
        ndef = 4
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
      end if
      if (n_pp.ge.2) then
        call quit(0,'set_r12f_targets',
     &              'check this definition (for n_pp==2)')
        descr = 'H[PX],[HP][HP]'
        ndef = 4 ! should be 6 ??? never ever tested, it seems ...
        ! 5
        occ_def(IHOLE,1,5) = 1
        occ_def(IPART,1,5) = 1
        occ_def(IPART,2,5) = 2
        ! 6
        occ_def(IHOLE,1,6) = 1
        occ_def(IEXTR,1,6) = 1
        occ_def(IPART,2,6) = 2
      end if
c      call op_from_occ_parameters(-1,parameters,2,
c     &     occ_def,ndef,1,(/      2,     0/),ndef)
c      call set_rule(op_rintbreve,ttype_op,DEF_OP_FROM_OCC,
c     &              op_rintbreve,1,1,
c     &              parameters,2,tgt_info)
      call set_rule2(op_rintbreve,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_rintbreve,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/op_rintbreve/))
      call set_arg(op_rintbreve,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/1/))
      call set_arg(op_rintbreve,DEF_OP_FROM_OCC,'CORE',2,tgt_info,
     &             val_int=(/2,0/))
      call set_arg(op_rintbreve,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)
      
      ! C k modified integrals r12tilde (cloning rint -> 2ext usually)
      call add_target(op_rinttilde,ttype_op,.false.,tgt_info)
      call set_dependency(op_rinttilde,op_rint,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rint,.false.) 
      call set_rule(op_rinttilde,ttype_op,CLONE_OP,
     &              op_rinttilde,1,1,
     &              parameters,1,tgt_info)
      
      ! commutator integrals <kl|r12[T1+T2,r12]|ij>
      !call add_target(op_rttr,ttype_op,.false.,tgt_info)
      call add_target2(op_rttr,.false.,tgt_info)
      occ_def = 0
      if (n_pp.ge.0) then
        descr = 'HH,,,HH'
        if (active_orbs) descr = '[HV][HV],,,[HV][HV]'
        ndef = 1
        occ_def(IHOLE,1,1) = 2
        occ_def(IHOLE,2,2) = 2
      end if
      if (n_pp.ge.1) then
        descr = 'H[HP],,,H[HP]'
        if (active_orbs) descr = '[HV][HVP],,,[HV][HVP]'
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
        descr = '[HP][HP],,,[HP][HP]'
        if (active_orbs) descr = '[HPV][HPV],,,[HPV][HPV]'
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
C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,2,(/ 0, 0, 0, 0/),ndef)
C      call set_rule(op_rttr,ttype_op,DEF_OP_FROM_OCC,
C     &              op_rttr,1,1,
C     &              parameters,2,tgt_info)
      call set_rule2(op_rttr,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_rttr,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/op_rttr/))
      call set_arg(op_rttr,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/2/))
      call set_arg(op_rttr,DEF_OP_FROM_OCC,'CORE',4,tgt_info,
     &             val_int=(/0,0,0,0/))
      call set_arg(op_rttr,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)

      ! (G.R)^{ij}_{pq}
      !call add_target(op_gr,ttype_op,.false.,tgt_info)
      call add_target2(op_gr,.false.,tgt_info)
      occ_def = 0
      descr = '[HP][HP],,,HH'
      if (active_orbs) descr = '[HPV][HPV],,,[HV][HV]'
      ! n_pp == 0:
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,2,2) = 2

      occ_def(IHOLE,1,3) = 1
      occ_def(IPART,1,3) = 1
      occ_def(IHOLE,2,4) = 2

      occ_def(IPART,1,5) = 2
      occ_def(IHOLE,2,6) = 2

      ! n_pp == 1:
      if (n_pp.ge.1) then
        descr = '[HP][HP],,,H[HP]'
        if (active_orbs) descr = '[HPV][HPV],,,[HV][HPV]'
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
      end if
      ! n_pp == 2:
      if (n_pp.ge.2) then
        descr = '[HP][HP],,,[HP][HP]'
        occ_def(IHOLE,1,13) = 2
        occ_def(IPART,2,14) = 2
        
        call quit(1,'set_r12f_targets','check this')
        ! ????
c      occ_def(IHOLE,1,15) = 1
c      occ_def(IPART,1,15) = 1
c      occ_def(IPART,2,16) = 2

        occ_def(IPART,1,17) = 2
        occ_def(IPART,2,18) = 2
      end if

      ndef = 3*(n_pp+1)

      if (t1ext.gt.0) then        
        descr = trim(descr)//'|PX,,,HH'
        occ_def(IPART,1,2*ndef+1) = 1
        occ_def(IEXTR,1,2*ndef+1) = 1
        occ_def(IHOLE,2,2*ndef+2) = 2
        ndef = ndef+1
      end if
      if (t1ext.ge.4) then        
        descr = trim(descr)//'|XX,,,HH'
        occ_def(IEXTR,1,2*ndef+1) = 2
        occ_def(IHOLE,2,2*ndef+2) = 2
        ndef = ndef+1
      end if
      if (.not.pf12_trunc) then        
        descr = trim(descr)//'|HX,,,HH'
        occ_def(IHOLE,1,2*ndef+1) = 1
        occ_def(IEXTR,1,2*ndef+1) = 1
        occ_def(IHOLE,2,2*ndef+2) = 2
        ndef = ndef+1
        if (n_pp.ge.1) then
          descr = trim(descr)//'|HX,,,HP'
          occ_def(IHOLE,1,2*ndef+1) = 1
          occ_def(IEXTR,1,2*ndef+1) = 1
          occ_def(IPART,2,2*ndef+2) = 1
          occ_def(IHOLE,2,2*ndef+2) = 1
          ndef = ndef+1
        end if
        if (n_pp.ge.2) then
          descr = trim(descr)//'|HX,,,PP'
          occ_def(IHOLE,1,2*ndef+1) = 1
          occ_def(IEXTR,1,2*ndef+1) = 1
          occ_def(IPART,2,2*ndef+2) = 2
          ndef = ndef+1
        end if
      end if

C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,2,(/ 0, 0, 0, 0/),3*(n_pp+1))
C      call set_rule(op_gr,ttype_op,DEF_OP_FROM_OCC,
C     &              op_gr,1,1,
C     &              parameters,2,tgt_info)
      call set_rule2(op_gr,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_gr,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/op_gr/))
      call set_arg(op_gr,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/2/))
      call set_arg(op_gr,DEF_OP_FROM_OCC,'CORE',4,tgt_info,
     &             val_int=(/0,0,0,0/))
      call set_arg(op_gr,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)

      ! extended variant
! NEXT
      !call add_target(op_gr_x,ttype_op,.false.,tgt_info)
      call add_target2(op_gr_x,.false.,tgt_info)
C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,2,(/ 2, 0, 2, 0/),3*(n_pp+1))
C      call set_rule(op_gr_x,ttype_op,DEF_OP_FROM_OCC,
C     &              op_gr_x,1,1,
C     &              parameters,2,tgt_info)
      call set_rule2(op_gr_x,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_gr_x,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/op_gr_x/))
      call set_arg(op_gr_x,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/2/))
      call set_arg(op_gr_x,DEF_OP_FROM_OCC,'CORE',4,tgt_info,
     &             val_int=(/2,0,2,0/))
      ! use descr as defined above for op_gr
      call set_arg(op_gr_x,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)

      ! core contributions (C part)
      occ_def = 0
      ndef = 2
      ! n_pp == 0:
      descr = 'H[HP],,,HH'
      if (active_orbs) descr = 'H[HPV],,,[HV][HV]'
      occ_def(IHOLE,1,1) = 2
      occ_def(IHOLE,2,2) = 2

      occ_def(IHOLE,1,3) = 1
      occ_def(IPART,1,3) = 1
      occ_def(IHOLE,2,4) = 2

      ! n_pp == 1:
      if (n_pp.ge.1) then
        ndef = 4
        descr = 'H[HP],,,H[HP]'
        if (active_orbs) descr = 'H[HPV],,,[HV][HPV]'
        occ_def(IHOLE,1,5) = 2
        occ_def(IHOLE,2,6) = 1
        occ_def(IPART,2,6) = 1

        occ_def(IHOLE,1,7) = 1
        occ_def(IPART,1,7) = 1
        occ_def(IHOLE,2,8) = 1
        occ_def(IPART,2,8) = 1

      end if
      ! n_pp == 2:
      if (n_pp.ge.2) then
        ndef = 5
        descr = 'H[HP],,,H[HP]|HH,,,PP'
        if (active_orbs) descr = 'H[HPV],,,[HV][HPV]|H[HV],,,PP'
        occ_def(IHOLE,1,9) = 2
        occ_def(IPART,2,10) = 2
        
      end if
      if (.not.pf12_trunc) then        
        if (.not.active_orbs) descr = trim(descr)//'|HX,,,HH'
        if (active_orbs) descr = trim(descr)//'|HX,,,[HV][HV]'
        occ_def(IHOLE,1,2*ndef+1) = 1
        occ_def(IEXTR,1,2*ndef+1) = 1
        occ_def(IHOLE,2,2*ndef+2) = 2
        ndef = ndef+1
        if (n_pp.ge.1) then
          if (.not.active_orbs) descr = trim(descr)//'|HX,,,HP'
          if (active_orbs) descr = trim(descr)//'|HX,,,[HV]P'
          occ_def(IHOLE,1,2*ndef+1) = 1
          occ_def(IEXTR,1,2*ndef+1) = 1
          occ_def(IPART,2,2*ndef+2) = 1
          occ_def(IHOLE,2,2*ndef+2) = 1
          ndef = ndef+1
        end if
        if (n_pp.ge.2) then
          descr = trim(descr)//'|HX,,,PP'
          occ_def(IHOLE,1,2*ndef+1) = 1
          occ_def(IEXTR,1,2*ndef+1) = 1
          occ_def(IPART,2,2*ndef+2) = 2
          ndef = ndef+1
        end if
      end if      
      !call add_target('G.R-Ccore',ttype_op,.false.,tgt_info)
      call add_target2('G.R-Ccore',.false.,tgt_info)
C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,2,(/ 2, 0, 2, 0/),2*ndef)
C      call set_rule('G.R-Ccore',ttype_op,DEF_OP_FROM_OCC,
C     &              'G.R-Ccore',1,1,
C     &              parameters,2,tgt_info)
      call set_rule2('G.R-Ccore',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('G.R-Ccore',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'G.R-Ccore'/))
      call set_arg('G.R-Ccore',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/2/))
      call set_arg('G.R-Ccore',DEF_OP_FROM_OCC,'CORE',4,tgt_info,
     &             val_int=(/2,0,2,0/))
      call set_arg('G.R-Ccore',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)


      ! R12^{2} integrals
      call add_target(op_ff,ttype_op,.false.,tgt_info)
      call set_dependency(op_ff,op_rttr,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rttr,.false.) ! <- dagger=.false.
      call set_rule(op_ff,ttype_op,CLONE_OP,
     &     op_ff,1,1,
     &     parameters,1,tgt_info)

c      if (is_keyword_set('method.CC').gt.0.and.(.not.truncate
c     &     .or.(truncate.and.trunc_type.gt.0)).or.
c     &     is_keyword_set('method.CCPT')) then
      !call add_target('R.R-X',ttype_op,.false.,tgt_info)
      call add_target2('R.R-X',.false.,tgt_info)

      descr = 'H[HPX],,,HH|HH,,,H[PX]'
      if (active_orbs) descr = 'H[HPVX],,,H[HV]|H[HV],,,H[PX]' ! not yet checked
      ndef = 5
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
        ! 5
      occ_def(IHOLE,1,9)  = 2
      occ_def(IHOLE,2,10) = 1
      occ_def(IEXTR,2,10) = 1
      if (n_pp.ge.1) then
        if (.not.active_orbs)
     &     descr = trim(descr)//'|H[PX],,,HP|HP,,,HX|PP,,,HH|PX,,,HH'
        if (active_orbs) 
     &     descr = trim(descr)//
     &                    '|H[PX],,,HP|HP,,,HX|PP,,,H[HV]|PX,,,H[HV]'
          
        ndef = 10
        ! 6
        occ_def(IHOLE,1,11)  = 1
        occ_def(IPART,1,11)  = 1
        occ_def(IHOLE,2,12)  = 1
        occ_def(IPART,2,12)  = 1
        ! 7
        occ_def(IHOLE,1,13)  = 1
        occ_def(IEXTR,1,13)  = 1
        occ_def(IHOLE,2,14)  = 1
        occ_def(IPART,2,14)  = 1
        ! 8
        occ_def(IHOLE,1,15)  = 1
        occ_def(IPART,1,15)  = 1
        occ_def(IHOLE,2,16) = 1
        occ_def(IEXTR,2,16) = 1
        ! 9
        occ_def(IPART,1,17)  = 2
        occ_def(IHOLE,2,18)  = 2
        ! 10
        occ_def(IPART,1,19)  = 1
        occ_def(IEXTR,1,19)  = 1
        occ_def(IHOLE,2,20)  = 2
      end if
C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,2,(/1,0,1,0/),ndef)
C      call set_rule('R.R-X',ttype_op,DEF_OP_FROM_OCC,
C     &              'R.R-X',1,1,
C     &              parameters,2,tgt_info)
      call set_rule2('R.R-X',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('R.R-X',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'R.R-X'/))
      call set_arg('R.R-X',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/2/))
      call set_arg('R.R-X',DEF_OP_FROM_OCC,'CORE',4,tgt_info,
     &             val_int=(/1,0,1,0/))
      call set_arg('R.R-X',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)

      ! {R12^2}BAR integrals
      call add_target(op_ffbar,ttype_op,.false.,tgt_info)
      call set_dependency(op_ffbar,op_rttr,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rttr,.false.) ! <- dagger=.false.
      call set_rule(op_ffbar,ttype_op,CLONE_OP,
     &              op_ffbar,1,1,
     &              parameters,1,tgt_info)

      ! special R12^{2} integrals
      !call add_target('FF-X',ttype_op,.false.,tgt_info)
      call add_target2('FF-X',.false.,tgt_info)
      occ_def = 0
      if (n_pp.ge.0) then
        ndef = 3
        descr = 'H[HPX],,,HH'
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
      end if
      if (n_pp.ge.1) then
        ndef = 10
        descr = '[HP][HPX],,,H[HP]'
        ! 4
        occ_def(IPART,1,7) = 2
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
        occ_def(IHOLE,1,13) = 1
        occ_def(IPART,1,13) = 1
        occ_def(IHOLE,2,14) = 1
        occ_def(IPART,2,14) = 1
        ! 8
        occ_def(IHOLE,1,15) = 1
        occ_def(IEXTR,1,15) = 1
        occ_def(IHOLE,2,16) = 1
        occ_def(IPART,2,16) = 1
        ! 9
        occ_def(IPART,1,17) = 2
        occ_def(IHOLE,2,18) = 1
        occ_def(IPART,2,18) = 1
        ! 10
        occ_def(IPART,1,19) = 1
        occ_def(IEXTR,1,19) = 1
        occ_def(IHOLE,2,20) = 1
        occ_def(IPART,2,20) = 1
      end if
C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,2,(/2,0,2,0/),ndef*2)
C      call set_rule('FF-X',ttype_op,DEF_OP_FROM_OCC,
C     &              'FF-X',1,1,
C     &              parameters,2,tgt_info)
      call set_rule2('FF-X',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('FF-X',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'FF-X'/))
      call set_arg('FF-X',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/2/))
      call set_arg('FF-X',DEF_OP_FROM_OCC,'CORE',4,tgt_info,
     &             val_int=(/2,0,2,0/))
      call set_arg('FF-X',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)
      
      ! R12^{2}*G12 integrals
      call add_target(op_ffg,ttype_op,.false.,tgt_info)
      call set_dependency(op_ffg,op_rttr,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     op_rttr,.false.)  
      call set_rule(op_ffg,ttype_op,CLONE_OP,
     &              op_ffg,1,1,
     &              parameters,1,tgt_info)

      ! V^{ij}_{pq}
      !call add_target(op_v_inter,ttype_op,.false.,tgt_info)
      call add_target2(op_v_inter,.false.,tgt_info)
      occ_def = 0
      ! for n_pp >= 0
      if (n_pp.ge.0) then
        ndef = 6
        descr = ',|[HP],H|[HP][HP],HH'
        if (.not.CC) descr = ',,'
        if (active_orbs) descr = ',|[HPV],[HV]|[HPV][HPV],[HV][HV]'
        if (active_orbs.and..not.CC) descr = ',|V,V|VV,VV'
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
        descr = ',|[HP],[HP]|[HP][HP],H[HP]'
c        descr = ',|[HP],H|[HP][HP],HH|[HP],P|[HP][HP],HP'
        if (active_orbs) descr = ',|[HPV],[HPV]|[HPV][HPV],[HV][HPV]'
c        if (active_orbs) descr = ',|[HPV],[HV]|[HPV][HPV],[HV][HV]
c     & |[HPV],P|[HPV][HPV],[HV]P'
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
        descr = ',|[HP],[HP]|[HP][HP],[HP][HP]'
        ndef = 14
        occ_def(IHOLE,1,12) = 2
        occ_def(IPART,2,12) = 2

        occ_def(IHOLE,1,13) = 1
        occ_def(IPART,1,13) = 1
        occ_def(IPART,2,13) = 2

        occ_def(IPART,1,14) = 2
        occ_def(IPART,2,14) = 2
      end if
      if (t1ext.gt.0) then
        descr = trim(descr)//'|X,H|PX,HH'
        occ_def(IEXTR,1,ndef+1) = 1
        occ_def(IHOLE,2,ndef+1) = 1
        ndef = ndef+1
        occ_def(IPART,1,ndef+1) = 1
        occ_def(IEXTR,1,ndef+1) = 1
        occ_def(IHOLE,2,ndef+1) = 2
        ndef = ndef+1
      end if
      if (t1ext.ge.4) then
        descr = trim(descr)//'|XX,HH'
        occ_def(IEXTR,1,ndef+1) = 2
        occ_def(IHOLE,2,ndef+1) = 2
        ndef = ndef+1
      end if
      if (.not.pf12_trunc.and..not.frozen) then
        descr = trim(descr)//'|HX,HH'
        occ_def(IHOLE,1,ndef+1) = 1
        occ_def(IEXTR,1,ndef+1) = 1
        occ_def(IHOLE,2,ndef+1) = 2
        ndef = ndef+1
        if (n_pp.ge.1) then
          descr = trim(descr)//'|HX,HP'
          occ_def(IHOLE,1,ndef+1) = 1
          occ_def(IEXTR,1,ndef+1) = 1
          occ_def(IPART,2,ndef+1) = 1
          occ_def(IHOLE,2,ndef+1) = 1
          ndef = ndef+1
        end if
        if (n_pp.ge.2) then
          descr = trim(descr)//'|HX,PP'
          occ_def(IHOLE,1,ndef+1) = 1
          occ_def(IEXTR,1,ndef+1) = 1
          occ_def(IPART,2,ndef+1) = 2
          ndef = ndef+1
        end if
      end if
C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,1,(/0,0/),ndef)
C      call set_rule(op_v_inter,ttype_op,DEF_OP_FROM_OCC,
C     &              op_v_inter,1,1,
C     &              parameters,2,tgt_info)
      call set_rule2(op_v_inter,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_v_inter,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/op_v_inter/))
      call set_arg(op_v_inter,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/1/))
      call set_arg(op_v_inter,DEF_OP_FROM_OCC,'CORE',2,tgt_info,
     &             val_int=(/0,0/))
      call set_arg(op_v_inter,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)

      ! extended variant (for formal expansion)
      !call add_target(op_v_x,ttype_op,.false.,tgt_info)
      call add_target2(op_v_x,.false.,tgt_info)
      occ_def = 0
      ! for n_pp >= 0
      if (n_pp.ge.0) then
        ndef = 3
        descr = '[HP][HP],HH'
        occ_def(IHOLE,1,1) = 2
        occ_def(IHOLE,2,1) = 2

        occ_def(IHOLE,1,2) = 1
        occ_def(IPART,1,2) = 1
        occ_def(IHOLE,2,2) = 2

        occ_def(IPART,1,3) = 2
        occ_def(IHOLE,2,3) = 2

      end if
      ! for n_pp >= 1
      if (n_pp.ge.1) then
        ndef = 6
        descr = '[HP][HP],H[HP]'
        occ_def(IHOLE,1,4) = 2
        occ_def(IHOLE,2,4) = 1
        occ_def(IPART,2,4) = 1

        occ_def(IHOLE,1,5) = 1
        occ_def(IPART,1,5) = 1
        occ_def(IHOLE,2,5) = 1
        occ_def(IPART,2,5) = 1

        occ_def(IPART,1,6) = 2
        occ_def(IHOLE,2,6) = 1
        occ_def(IPART,2,6) = 1
      end if
      ! for n_pp >= 2
      if (n_pp.ge.2) then
        ndef = 9
        descr = '[HP][HP],[HP][HP]'
        occ_def(IHOLE,1,7) = 2
        occ_def(IPART,2,7) = 2

        occ_def(IHOLE,1,8) = 1
        occ_def(IPART,1,8) = 1
        occ_def(IPART,2,8) = 2

        occ_def(IPART,1,9) = 2
        occ_def(IPART,2,9) = 2
      end if
      if (.not.pf12_trunc) then
        descr = trim(descr)//'|HX,HH'
        occ_def(IHOLE,1,ndef+1) = 1
        occ_def(IEXTR,1,ndef+1) = 1
        occ_def(IHOLE,2,ndef+1) = 2
        ndef = ndef+1
        if (n_pp.ge.1) then
          descr = trim(descr)//'|HX,HP'
          occ_def(IHOLE,1,ndef+1) = 1
          occ_def(IEXTR,1,ndef+1) = 1
          occ_def(IPART,2,ndef+1) = 1
          occ_def(IHOLE,2,ndef+1) = 1
          ndef = ndef+1
        end if
        if (n_pp.ge.2) then
          descr = trim(descr)//'|HX,PP'
          occ_def(IHOLE,1,ndef+1) = 1
          occ_def(IEXTR,1,ndef+1) = 1
          occ_def(IPART,2,ndef+1) = 2
          ndef = ndef+1
        end if
      end if
C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,1,(/2,0/),ndef)
C      call set_rule(op_v_x,ttype_op,DEF_OP_FROM_OCC,
C     &              op_v_x,1,1,
C     &              parameters,2,tgt_info)
      call set_rule2(op_v_x,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_v_x,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/op_v_x/))
      call set_arg(op_v_x,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/1/))
      call set_arg(op_v_x,DEF_OP_FROM_OCC,'CORE',2,tgt_info,
     &             val_int=(/2,0/))
      call set_arg(op_v_x,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)


      call add_target('V-Ccore',ttype_op,.false.,tgt_info)
      occ_def = 0
      ! for n_pp >= 0
      if (n_pp.ge.0) then
        ndef = 2
        descr = 'H[HP],HH'
        occ_def(IHOLE,1,1) = 2
        occ_def(IHOLE,2,1) = 2

        occ_def(IHOLE,1,2) = 1
        occ_def(IPART,1,2) = 1
        occ_def(IHOLE,2,2) = 2
      end if
      ! for n_pp >= 1
      if (n_pp.ge.1) then
        ndef = 4
        descr = 'H[HP],H[HP]'
        occ_def(IHOLE,1,3) = 2
        occ_def(IHOLE,2,3) = 1
        occ_def(IPART,2,3) = 1

        occ_def(IHOLE,1,4) = 1
        occ_def(IPART,1,4) = 1
        occ_def(IHOLE,2,4) = 1
        occ_def(IPART,2,4) = 1
      end if
      ! for n_pp >= 2
      if (n_pp.ge.2) then
        ndef = 6
        descr = 'H[HP],[HP][HP]'
        occ_def(IHOLE,1,5) = 2
        occ_def(IPART,2,5) = 2

        occ_def(IHOLE,1,6) = 1
        occ_def(IPART,1,6) = 1
        occ_def(IPART,2,6) = 2
      end if
      if (.not.pf12_trunc) then
        descr = trim(descr)//'|HX,HH'
        occ_def(IHOLE,1,ndef+1) = 1
        occ_def(IEXTR,1,ndef+1) = 1
        occ_def(IHOLE,2,ndef+1) = 2
        ndef = ndef+1
        if (n_pp.ge.1) then
          descr = trim(descr)//'|HX,HP'
          occ_def(IHOLE,1,ndef+1) = 1
          occ_def(IEXTR,1,ndef+1) = 1
          occ_def(IPART,2,ndef+1) = 1
          occ_def(IHOLE,2,ndef+1) = 1
          ndef = ndef+1
        end if
        if (n_pp.ge.2) then
          descr = trim(descr)//'|HX,PP'
          occ_def(IHOLE,1,ndef+1) = 1
          occ_def(IEXTR,1,ndef+1) = 1
          occ_def(IPART,2,ndef+1) = 2
          ndef = ndef+1
        end if
      end if
C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,1,(/2,0/),ndef)
C      call set_rule('V-Ccore',ttype_op,DEF_OP_FROM_OCC,
C     &              'V-Ccore',1,1,
C     &              parameters,2,tgt_info)
      call set_rule2('V-Ccore',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('V-Ccore',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'V-Ccore'/))
      call set_arg('V-Ccore',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('V-Ccore',DEF_OP_FROM_OCC,'CORE',2,tgt_info,
     &             val_int=(/2,0/))
      call set_arg('V-Ccore',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &             val_str=descr)

      ! V' intermediate
      if (vring_mode.eq.2) then
        nblk = 6
        nj = 1
        occ_def = 0
        descr = 'P,[HP]|[HP]P,H[HP]'
        !
        occ_def(IPART,1,1) = 1
        occ_def(IHOLE,2,1) = 1
        !
        occ_def(IPART,1,2) = 1
        occ_def(IPART,2,2) = 1
        !
        occ_def(IPART,1,3) = 2
        occ_def(IHOLE,2,3) = 2
        !
        occ_def(IPART,1,4) = 1
        occ_def(IHOLE,2,4) = 1
        occ_def(IHOLE,1,4) = 1
        occ_def(IPART,2,4) = 1
        !
        occ_def(IHOLE,1,5) = 1
        occ_def(IPART,1,5) = 1
        occ_def(IHOLE,2,5) = 2
        !
        occ_def(IPART,1,6) = 2
        occ_def(IHOLE,2,6) = 1
        occ_def(IPART,2,6) = 1
      else if (vring_mode.eq.1) then
        descr = 'PP,HH'
        nblk = 1
        nj = 1
        occ_def = 0
        !
        occ_def(IPART,1,1) = 2
        occ_def(IHOLE,2,1) = 2
      else
        descr = 'P,H|PP,HH'
        if (active_orbs) descr = 'P,[HV]|PP,[HV][HV]|VP,HH'
        if (n_pp.eq.1) then
          descr = 'P,H|PP,HH'
          if (active_orbs) descr = 'P,[HV]|PP,[HV][HPV]|VP,[PH]H'
        end if
        nblk = 2
        nj = 1
        occ_def = 0
        occ_def(IPART,1,1) = 1
        occ_def(IHOLE,2,1) = 1
        occ_def(IPART,1,2) = 2
        occ_def(IHOLE,2,2) = 2
      end if
      call add_target2(op_vp_inter,.false.,tgt_info)
      call set_dependency(op_vp_inter,op_v_inter,tgt_info)
      call set_rule2(op_vp_inter,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_vp_inter,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/op_vp_inter/))
C      call set_arg(op_vp_inter,DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
C     &     val_int=(/nblk/))
C      call set_arg(op_vp_inter,DEF_OP_FROM_OCC,'OCC',nblk*nj,tgt_info,
C     &     val_occ=occ_def(1:ngastp,1:2,1:nblk*nj))
      call set_arg(op_vp_inter,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)
      
      ! V' intermediate
c      if () then
c        nblk = 1
c        occ_def = 0
c        occ_def(IPART,1,1) = 2
c        occ_def(IHOLE,2,1) = 2
        nblk = 1
        nj   = 2
        occ_def = 0
        occ_def(IPART,1,1) = 1
        occ_def(IHOLE,2,1) = 1
        occ_def(IPART,1,2) = 1
        occ_def(IHOLE,2,2) = 1
      descr = 'P,H,P,H'
c      else
c      end if
      call add_target2('VR2',.false.,tgt_info)
      call set_dependency('VR2',op_v_inter,tgt_info)
      call set_rule2('VR2',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('VR2',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'VR2'/))
C      call set_arg('VR2',DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
C     &     val_int=(/nblk/))
      call set_arg('VR2',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/2/))
C      call set_arg('VR2',DEF_OP_FROM_OCC,'OCC',nblk*nj,tgt_info,
C     &     val_occ=occ_def(1:ngastp,1:2,1:nblk*nj))
      call set_arg('VR2',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)
      
            
      ! B intermediate
      occ_def = 0
      if (n_pp.ge.0) then
        ndef = 1
        descr = ',,'
        if (active_orbs) descr = ',|[HV],[HV]|[HV]V,[HV]V'
c        if (active_orbs) descr = ',|V,V|VV,VV'
      end if
      if (n_pp.ge.1) then
        descr = ',|H,P|P,H|P,P|HP,HH|H[HP],HP'
c        descr = ',|P,P|HP,HP'  !for triples only
c dbg
c        if (active_orbs) descr = '?'
        if (active_orbs) descr = ',|V,V|VV,VV|P,P|[HV]P,[HV]P'
c dbgend
c        ndef = 7
c        occ_def(IHOLE,1,2) = 1
c        occ_def(IPART,2,2) = 1
c
c        occ_def(IPART,1,3) = 1
c        occ_def(IHOLE,2,3) = 1
c
c        occ_def(IPART,1,4) = 1
c        occ_def(IPART,2,4) = 1
c
c        occ_def(IHOLE,1,5) = 1
c        occ_def(IPART,1,5) = 1
c        occ_def(IHOLE,2,5) = 2
c
c        occ_def(IHOLE,1,6) = 2
c        occ_def(IHOLE,2,6) = 1
c        occ_def(IPART,2,6) = 1
c
c        occ_def(IHOLE,1,7) = 1
c        occ_def(IPART,1,7) = 1
c        occ_def(IHOLE,2,7) = 1
c        occ_def(IPART,2,7) = 1
      end if
      if (n_pp.ge.2) then
        descr = '[HP]P,HH|[HP][HP],[HP]P'
        if (active_orbs) descr = '?'
c        occ_def(IHOLE,1,ndef+1) = 2
c        occ_def(IPART,2,ndef+1) = 2
c
c        occ_def(IPART,1,ndef+2) = 2
c        occ_def(IHOLE,2,ndef+2) = 2
c
c        occ_def(IHOLE,1,ndef+3) = 1
c        occ_def(IPART,1,ndef+3) = 1
c        occ_def(IPART,2,ndef+3) = 2
c
c        occ_def(IPART,1,ndef+4) = 2
c        occ_def(IHOLE,2,ndef+4) = 1
c        occ_def(IPART,2,ndef+4) = 1
c
c        occ_def(IPART,1,ndef+5) = 2
c        occ_def(IPART,2,ndef+5) = 2
c        ndef = ndef + 5
      end if
      !call add_target(op_b_inter,ttype_op,.false.,tgt_info)
      call add_target2(op_b_inter,.false.,tgt_info)
C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,1,(/0,0/),ndef)
C      call set_rule(op_b_inter,ttype_op,DEF_OP_FROM_OCC,
C     &              op_b_inter,1,1,
C     &              parameters,2,tgt_info)
      call set_rule2(op_b_inter,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_b_inter,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/op_b_inter/))
      call set_arg(op_b_inter,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/1/))
      call set_arg(op_b_inter,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)

      ! Bhole intermediate
      call add_target(op_bh_inter,ttype_op,.false.,tgt_info)
c      if (n_pp.eq.0) then
c        descr=','
c        if (active_orbs) descr=',|V,V'
      call set_dependency(op_bh_inter,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
      call set_rule(op_bh_inter,ttype_op,CLONE_OP,
     &              op_bh_inter,1,1,
     &              parameters,1,tgt_info)
c      call set_rule2(op_bh_inter,DEF_OP_FROM_OCC,tgt_info)
c      call set_arg(op_bh_inter,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
c     &     val_label=(/op_bh_inter/))
c      call set_arg(op_bh_inter,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg(op_bh_inter,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
c     &     val_str=descr)
c      else if (n_pp.eq.1) then
c      descr=',|P,P'
c      if (active_orbs) descr=',|V,V|[HVP],P|P,[HV]'
c      call set_rule2(op_bh_inter,DEF_OP_FROM_OCC,tgt_info)
c      call set_arg(op_bh_inter,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
c     &     val_label=(/op_bh_inter/))
c      call set_arg(op_bh_inter,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg(op_bh_inter,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
c     &     val_str=descr)
c      end if

      ! X intermediate
      call add_target(op_x_inter,ttype_op,.false.,tgt_info)
      if (n_pp.eq.0) then
      call set_dependency(op_x_inter,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
      call set_rule(op_x_inter,ttype_op,CLONE_OP,
     &              op_x_inter,1,1,
     &              parameters,1,tgt_info)
c      descr=','
c      if (active_orbs) descr='V,[HV]|H,V|V[HV],VV|VV,HV'
c      call set_rule2(op_x_inter,DEF_OP_FROM_OCC,tgt_info)
c      call set_arg(op_x_inter,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
c     &     val_label=(/op_x_inter/))
c      call set_arg(op_x_inter,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg(op_x_inter,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
c     &     val_str=descr)
      else if (n_pp.eq.1) then
      descr=',|H,P|P,H|P,P|HP,HH|H[HP],HP'
c      descr='H,P|P,[HP]|HH,HP|HP,[HP]H'
      if (active_orbs) descr=',|[HVP],[HVP]|[HV][HVP],[HV][HVP]'
      call set_rule2(op_x_inter,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_x_inter,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/op_x_inter/))
      call set_arg(op_x_inter,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/1/))
      call set_arg(op_x_inter,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)
      end if

      ! X' intermediate
      !call add_target(op_xp_inter,ttype_op,.false.,tgt_info)
      call add_target2(op_xp_inter,.false.,tgt_info)
      ndef = 4
      descr = 'P,P|HP,[HP]P|PP,HP'
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

C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,1,(/0,0/),ndef)
C      call set_rule(op_xp_inter,ttype_op,DEF_OP_FROM_OCC,
C     &              op_xp_inter,1,1,
C     &              parameters,2,tgt_info)
      call set_rule2(op_xp_inter,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_xp_inter,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/op_xp_inter/))
      call set_arg(op_xp_inter,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/1/))
      call set_arg(op_xp_inter,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)

      ! Xh intermediate
      !call add_target(op_xh_inter,ttype_op,.false.,tgt_info)
      call add_target2(op_xh_inter,.false.,tgt_info)
      ndef = 4
      descr = 'H,,,H|HH,,,HP|HP,,,HH|HP,,,HP'
      if (.not.CC) descr = 'H,,,H'
      if (active_orbs.and..not.CC) 
     &   descr = 'H,,,H|H,,,V|V,,,H|HV,,,HV|HV,,,VV|VV,,,HV'
      if (active_orbs.and.n_pp.eq.1) 
     &   descr = 'H,,,[HV]|V,,,H|HV,,,[HV]V|VV,,,HV|[HV],,,P
     &      |P,,,[HV]|P[HV],,,[HV][HV]|[HV][HV],,,P[HV]'
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

      if (.not.pf12_trunc) then
        descr = trim(descr)//'|HH,,,HH'
        occ_def(IHOLE,1,9) = 2
        occ_def(IHOLE,2,10) = 2
        ndef = 5
      end if

C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,2,(/0,0,0,0/),ndef*2)
C      call set_rule(op_xh_inter,ttype_op,DEF_OP_FROM_OCC,
C     &              op_xh_inter,1,1,
C     &              parameters,2,tgt_info)
      call set_rule2(op_xh_inter,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_xh_inter,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/op_xh_inter/))
      call set_arg(op_xh_inter,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/2/))
      call set_arg(op_xh_inter,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)

      ! C intermediate
      !call add_target(op_c_inter,ttype_op,.false.,tgt_info)
      call add_target2(op_c_inter,.false.,tgt_info)
      occ_def = 0
      if (n_pp.eq.0) descr = 'PP,HH'
      if (n_pp.eq.1) descr = 'PP,H[HP]'
      if (n_pp.eq.2) descr = 'PP,[HP][HP]'
      if (active_orbs) then
        if (n_pp.eq.0) descr = 'P[PV],[HV][HV]'
        if (n_pp.eq.1) descr = 'P[PV],[HV][HV]|P[PV],[HV]P'
        if (n_pp.eq.2) descr = 'P[PV],[HPV][HPV]'
      end if
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

C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,n_pp+1,1,(/0,0/),6)
C      call set_rule(op_c_inter,ttype_op,DEF_OP_FROM_OCC,
C     &              op_c_inter,1,1,
C     &              parameters,2,tgt_info)
      call set_rule2(op_c_inter,DEF_OP_FROM_OCC,tgt_info)
      call set_arg(op_c_inter,DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/op_c_inter/))
      call set_arg(op_c_inter,DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/1/))
      call set_arg(op_c_inter,DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)

      ! C1 intermediate f^k_x r^{ax}_{ik}
      nblk = 1
      nj   = 1
      occ_def = 0
      descr = 'P,H'
      if (active_orbs) descr = 'P,[HV]'
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,1) = 1
      call add_target2('C1',.false.,tgt_info)
      call set_rule2('C1',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('C1',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'C1'/))
C      call set_arg('C1',DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
C     &     val_int=(/nblk/))
      call set_arg('C1',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/1/))
C      call set_arg('C1',DEF_OP_FROM_OCC,'OCC',nblk*nj,tgt_info,
C     &     val_occ=occ_def(1:ngastp,1:2,1:nblk*nj))
      call set_arg('C1',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)
 

      ! P intermediate
      call add_target(op_p_inter,ttype_op,.false.,tgt_info)
      call set_dependency(op_p_inter,op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_b_inter,.false.) ! <- dagger=.false.
      call set_rule(op_p_inter,ttype_op,CLONE_OP,
     &              op_p_inter,1,1,
     &              parameters,1,tgt_info)

      ! Vpx operator
      call add_target('Vpx',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 3
      occ_def(IEXTR,1,1) = 1
      occ_def(IPART,2,1) = 1
      occ_def(IHOLE,1,2) = 1
      occ_def(IEXTR,1,2) = 1
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1
      occ_def(IHOLE,1,3) = 1
      occ_def(IEXTR,1,3) = 1
      occ_def(IHOLE,2,3) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('Vpx',ttype_op,DEF_OP_FROM_OCC,
     &              'Vpx',1,1,
     &              parameters,2,tgt_info)
      
      !call add_target('G.R-X',ttype_op,.false.,tgt_info)
      call add_target2('G.R-X',.false.,tgt_info)

      occ_def = 0
      descr = 'HX,,,H[HP]'
      if (active_orbs) descr = '[HV]X,,,[HV][HVP]'
      ndef = 2
      occ_def(IHOLE,1,1) = 1
      occ_def(IEXTR,1,1) = 1
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1
      occ_def(IHOLE,1,3) = 1
      occ_def(IEXTR,1,3) = 1
      occ_def(IHOLE,2,4) = 2
C      call op_from_occ_parameters(-1,parameters,2,
C     &     occ_def,ndef,2,(/0,0,0,0/),ndef)
C      call set_rule('G.R-X',ttype_op,DEF_OP_FROM_OCC,
C     &              'G.R-X',1,1,
C     &              parameters,2,tgt_info)
      call set_rule2('G.R-X',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('G.R-X',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'G.R-X'/))
      call set_arg('G.R-X',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/2/))
      call set_arg('G.R-X',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)

      ! Z intermediate
      call add_target(op_z_inter,ttype_op,.false.,tgt_info)
c      ndef = 1
c      occ_def = 0
c      occ_def(IHOLE,1,1) = 1
c      occ_def(IPART,2,1) = 1
      ndef = 2
      occ_def = 0
      occ_def(IHOLE,1,2) = 1
      occ_def(IPART,2,2) = 1
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/0,0/),6)
      call set_rule(op_z_inter,ttype_op,DEF_OP_FROM_OCC,
     &              op_z_inter,1,1,
     &              parameters,2,tgt_info)


      ! Z2 intermediate (for R^+ R couplings)
      call add_target2('Z2-INT',.false.,tgt_info)
      occ_def = 0
      ! 1
      occ_def(IHOLE,1,1) = 2
      occ_def(IPART,2,1) = 2
      ndef = 1
      descr='HH,PP'
      if (active_orbs) descr='[HV][HV],[PV]P|[HV][HV]V,[PV]PV'
      if (max_rank.gt.3.or.(min_rank_tp.eq.1.and.max_rank.gt.2)) then
        occ_def(IHOLE,1,2) = 1
        occ_def(IPART,1,2) = 1
        occ_def(IPART,2,2) = 2
        ! 3
        occ_def(IHOLE,1,3) = 2
        occ_def(IPART,1,3) = 1
        occ_def(IHOLE,2,3) = 1
        occ_def(IPART,2,3) = 2
        ndef = 3   
        descr='HH,PP|HP,PP|HP,HP'
      end if
c      call op_from_occ_parameters(-1,parameters,2,
c     &     occ_def,ndef,1,(/0,0/),6)
c      call set_rule('Z2-INT',ttype_op,DEF_OP_FROM_OCC,
c     &              'Z2-INT',1,1,
c     &              parameters,2,tgt_info)
      call set_rule2('Z2-INT',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('Z2-INT',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'Z2-INT'/))
      call set_arg('Z2-INT',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('Z2-INT',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)

        ! Non-anti-symmetrised Hamiltonian integrals.
        call add_target(op_g_z,ttype_op,.false.,tgt_info)
c        call set_g_z_old(ndef,occ_def)
        call set_g_z(ndef,occ_def)
        call op_from_occ_parameters(-1,parameters,2,
     &       occ_def,ndef,2,(/1,1,0,0/),2*ndef)
        call set_rule(op_g_z,ttype_op,DEF_OP_FROM_OCC,
     &                op_g_z,1,1,
     &                parameters,2,tgt_info)

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
c dbg
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule(form_r12_vint,ttype_frm,PRINT_FORMULA,
c     &              labels,2,1,
c     &              parameters,2,tgt_info)
c dbg

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
      else if (t1ext.eq.1) then
        call set_dependency(form_r12_vcabs,'H-ext',tgt_info)
        nint = 2
        labels(3) = op_g_x//'^+'
        labels(4) = 'H-ext'        
        labels(5) = op_g_x//'^+'
        labels(6) = op_ham
      end if
      call form_parameters(-1,
     &     parameters,2,title_r12_vcabs,nint,'---')
      call set_rule(form_r12_vcabs,ttype_frm,REPLACE,
     &              labels,2+nint*2,1,
     &              parameters,2,tgt_info)
c dbg
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule(form_r12_vcabs,ttype_frm,PRINT_FORMULA,
c     &              labels,2,1,
c     &              parameters,2,tgt_info)
c dbg

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

      ! formal definition of VR
      call add_target2('Vring_formal',.false.,tgt_info)
      if (.not.active_orbs) then
      call set_dependency('Vring_formal',op_vp_inter,tgt_info)
      call set_dependency('Vring_formal','VR2',tgt_info)
      call set_dependency('Vring_formal',op_ham,tgt_info)
      call set_dependency('Vring_formal',op_r12,tgt_info)
      call set_rule2('Vring_formal',DEF_R12INTM_FORMAL,tgt_info)
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,'LABEL',1,tgt_info,
     &             val_label=(/'Vring_formal'/))
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,
     &                                              'INTERM',1,tgt_info,
     &             val_label=(/op_vp_inter/))
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,
     &                                           'OPERATORS',2,tgt_info,
     &             val_label=(/op_r12,op_ham/))
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,
     &                                           'ANSATZ',1,tgt_info,
     &             val_int=(/ansatz/))
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,
     &                                           'MODE',1,tgt_info,
     &             val_str='VR')
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,
     &                                           'TITLE',1,tgt_info,
     &            val_str='V(ring) intermediate, formal definition')
      ! do the same for non-symmetrized ia,jb
      call set_rule2('Vring_formal',DEF_R12INTM_FORMAL,tgt_info)
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,'LABEL',1,tgt_info,
     &             val_label=(/'Vring2_formal'/))
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,
     &                                              'INTERM',1,tgt_info,
     &             val_label=(/'VR2'/))
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,
     &                                           'OPERATORS',2,tgt_info,
     &             val_label=(/op_r12,op_ham/))
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,
     &                                           'ANSATZ',1,tgt_info,
     &             val_int=(/ansatz/))
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,
     &                                           'MODE',1,tgt_info,
     &             val_str='VR')
      call set_arg('Vring_formal',DEF_R12INTM_FORMAL,
     &                                           'TITLE',1,tgt_info,
     &            val_str='V(ring) intermediate, formal definition')
      else
      call set_dependency('Vring_formal','H',tgt_info)
      call set_dependency('Vring_formal','R12',tgt_info)
      call set_dependency('Vring_formal','U',tgt_info)
      call set_rule2('Vring_formal',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('Vring_formal',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'Vring_formal'/))
      call set_arg('Vring_formal',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'U'/))
      call set_arg('Vring_formal',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'U  ','H  ','R12','U  '/))
      call set_arg('Vring_formal',EXPAND_OP_PRODUCT,'N_DESCR',1,
     &     tgt_info,val_int=(/2/))
      call set_arg('Vring_formal',EXPAND_OP_PRODUCT,'DESCR',2,tgt_info,
     &     val_label=(/'2,3,H,X','2,,[HP][HVP],[HVX][HV]'/))
      call set_arg('Vring_formal',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,3,1/))
c dbg
c      call set_rule2('Vring_formal',PRINT_FORMULA,tgt_info)
c      call set_arg('Vring_formal',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'Vring_formal'/))
c dbgend
      end if

      call add_target2('Vring_CABS',.false.,tgt_info)
      call set_dependency('Vring_CABS',op_ham,tgt_info)
      call set_dependency('Vring_CABS','Vring_formal',tgt_info)
      call set_dependency('Vring_CABS',op_rint,tgt_info)
      call set_rule2('Vring_CABS',REPLACE,tgt_info)
      call set_arg('Vring_CABS',REPLACE,'LABEL_RES',1,tgt_info,
     &            val_label=(/'Vring_CABS'/))
      call set_arg('Vring_CABS',REPLACE,'LABEL_IN',1,tgt_info,
     &            val_label=(/'Vring_formal'/))
      call set_arg('Vring_CABS',REPLACE,'OP_LIST',2,tgt_info,
     &            val_label=(/op_r12,op_rint/))
      call set_arg('Vring_CABS',REPLACE,'TITLE',1,tgt_info,
     &            val_str='V(ring) intermediate, for evaluation')
      if(.not.active_orbs) then
        call set_rule2('Vring_CABS',REPLACE,tgt_info)
        call set_arg('Vring_CABS',REPLACE,'LABEL_RES',1,tgt_info,
     &              val_label=(/'Vring2_CABS'/))
        call set_arg('Vring_CABS',REPLACE,'LABEL_IN',1,tgt_info,
     &              val_label=(/'Vring2_formal'/))
        call set_arg('Vring_CABS',REPLACE,'OP_LIST',2,tgt_info,
     &              val_label=(/op_r12,op_rint/))
        call set_arg('Vring_CABS',REPLACE,'TITLE',1,tgt_info,
     &              val_str='V(ring) intermediate, for evaluation')
        ! for some exceptional cases (no CABS) make sure
        ! that R12 is removed 
        call set_rule2('Vring_CABS',INVARIANT,tgt_info)
        call set_arg('Vring_CABS',INVARIANT,'LABEL_RES',1,tgt_info,
     &              val_label=(/'Vring_CABS'/))
        call set_arg('Vring_CABS',INVARIANT,'LABEL_IN',1,tgt_info,
     &              val_label=(/'Vring_CABS'/))
        call set_arg('Vring_CABS',INVARIANT,'OP_RES',1,tgt_info,
     &              val_label=(/op_vp_inter/))
        call set_arg('Vring_CABS',INVARIANT,'OPERATORS',1,tgt_info,
     &              val_label=(/op_r12/))
        call set_arg('Vring_CABS',INVARIANT,'TITLE',1,tgt_info,
     &              val_str='V(ring) intermediate, for evaluation')

        call set_rule2('Vring_CABS',INVARIANT,tgt_info)
        call set_arg('Vring_CABS',INVARIANT,'LABEL_RES',1,tgt_info,
     &              val_label=(/'Vring2_CABS'/))
        call set_arg('Vring_CABS',INVARIANT,'LABEL_IN',1,tgt_info,
     &              val_label=(/'Vring2_CABS'/))
        call set_arg('Vring_CABS',INVARIANT,'OP_RES',1,tgt_info,
     &              val_label=(/'VR2'/))
        call set_arg('Vring_CABS',INVARIANT,'OPERATORS',1,tgt_info,
     &              val_label=(/op_r12/))
        call set_arg('Vring_CABS',INVARIANT,'TITLE',1,tgt_info,
     &              val_str='V(ring) intermediate, for evaluation')
      else
      ! adding terms to U
        call set_rule2('Vring_CABS',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'Vring_CABS'/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'U'/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'OPERATORS',6,
     &       tgt_info,
     &       val_label=(/'U   ','C0^+','H   ',
     &                   'R12 ','C0  ','U   '/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'N_DESCR',1,
     &       tgt_info,val_int=(/4/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'DESCR',4,tgt_info,
     &       val_label=(/'3,4,,X','2,3,,V','4,5,,V',
     &                   '3,,[HVP]V,[HVX][HV]'/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'IDX_SV',6,tgt_info,
     &       val_int=(/1,2,3,4,5,1/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'N_AVOID',1,
     &       tgt_info,val_int=(/2/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'AVOID',4,
     &       tgt_info,val_int=(/2,6,1,5/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'N_CONNECT',1,
     &       tgt_info,val_int=(/1/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'CONNECT',2,
     &       tgt_info,val_int=(/1,3/))
      ! adding terms to U again
        call set_rule2('Vring_CABS',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'Vring_CABS'/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'U'/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'OPERATORS',6,
     &       tgt_info,
     &       val_label=(/'U   ','C0^+','H   ',
     &                   'R12 ','C0  ','U   '/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'N_DESCR',1,
     &       tgt_info,val_int=(/4/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'DESCR',4,tgt_info,
     &       val_label=(/'3,4,,X','2,3,,V','4,5,,V',
     &                   '3,,[HVP]V,[HVX][HV]'/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'IDX_SV',6,tgt_info,
     &       val_int=(/1,2,3,4,5,1/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'N_AVOID',1,
     &       tgt_info,val_int=(/3/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'AVOID',6,
     &       tgt_info,val_int=(/2,6,1,5,1,3/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'N_CONNECT',1,
     &       tgt_info,val_int=(/1/))
        call set_arg('Vring_CABS',EXPAND_OP_PRODUCT,'CONNECT',2,
     &       tgt_info,val_int=(/3,6/))
        call set_rule2('Vring_CABS',REPLACE,tgt_info)
        call set_dependency('Vring_CABS','R12-INT',tgt_info)
        call set_arg('Vring_CABS',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'Vring_CABS'/))
        call set_arg('Vring_CABS',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'Vring_CABS'/))
        call set_arg('Vring_CABS',REPLACE,'OP_LIST',4,tgt_info,
     &       val_label=(/'R12      ','R12-INT  ',
     &                   'R12^+    ','R12-INT^+'/))
c dbg
c      call set_rule2('Vring_CABS',PRINT_FORMULA,tgt_info)
c      call set_arg('Vring_CABS',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'Vring_CABS'/))
c dbgend
      end if
      ! formal definition of Vpx
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Vpx_formal'
      labels(2) = 'Vpx'
      labels(3) = op_r12
      labels(4) = op_ham
      call add_target('Vpx_formal',ttype_frm,.false.,tgt_info)
      call set_dependency('Vpx_formal','Vpx',tgt_info)
      call set_dependency('Vpx_formal',op_ham,tgt_info)
      call set_dependency('Vpx_formal',op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vint,0,'V')
      call set_rule('Vpx_formal',ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)

      ! CABS approximation to V
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Vpx_CABS'
      labels(2) = 'Vpx'
      labels(3) = op_g_x !op_ham
      labels(4) = op_rint
      labels(5) = 'G.R-X'   
      ! F12: op_gr
      call add_target('Vpx_CABS',ttype_frm,.false.,tgt_info)
      call set_dependency('Vpx_CABS',op_v_inter,tgt_info)
      call set_dependency('Vpx_CABS','G.R-X',tgt_info)
      call set_dependency('Vpx_CABS',op_g_x,tgt_info)
      call set_dependency('Vpx_CABS',op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_vcabs,ansatz,'V '//approx)
      call set_rule('Vpx_CABS',ttype_frm,DEF_R12INTM_CABS,
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
c dbg
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule(form_r12_xint,ttype_frm,PRINT_FORMULA,
c     &              labels,2,1,
c     &              parameters,2,tgt_info)
c dbg

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
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule(form_r12_xcabs,ttype_frm,PRINT_FORMULA,
c     &              labels,2,1,
c     &              parameters,2,tgt_info)
c dbg

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
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule(form_r12_xhint,ttype_frm,PRINT_FORMULA,
c     &              labels,2,1,
c     &              parameters,2,tgt_info)
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
c dbg
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule(form_r12_xhcabs,ttype_frm,PRINT_FORMULA,
c     &              labels,2,1,
c     &              parameters,2,tgt_info)
c dbg

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
      labels(5) = op_r12//'^+'
      labels(6) = op_rint//'^+'
      call add_target(form_r12_xpcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_xcabs,op_xp_inter,tgt_info)
      call set_dependency(form_r12_xcabs,op_r12,tgt_info)
      call set_dependency(form_r12_xcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_xpcabs,2,'---')
      call set_rule(form_r12_xpcabs,ttype_frm,REPLACE,
     &              labels,6,1,
     &              parameters,2,tgt_info)

      ! formal definition of B
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_bint
      labels(2) = op_b_inter
      labels(3) = op_r12
      labels(4) = op_ham
      call add_target(form_r12_bint,ttype_frm,.false.,tgt_info)
      if (active_orbs)
     &    call set_dependency(form_r12_bint,'Favg',tgt_info)
      if (active_orbs) labels(4) = 'Favg'
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
        if (active_orbs) 
     &      call set_dependency(form_r12_bcabs,'Favg',tgt_info)
        if (active_orbs) labels(7) = 'Favg'
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
c dbg
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule(form_r12_bcabs,ttype_frm,PRINT_FORMULA,
c     &              labels,2,1,
c     &              parameters,2,tgt_info)
c dbg

      ! formal definition of Bh
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r12_bhint
      labels(2) = op_bh_inter
      labels(3) = op_r12
      labels(4) = op_ham
      call add_target(form_r12_bhint,ttype_frm,.false.,tgt_info)
      if (active_orbs)
     &    call set_dependency(form_r12_bhint,'Favg',tgt_info)
      if (active_orbs) labels(4) = 'Favg'
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
      if (active_orbs) labels(7) = 'Favg'
      call add_target(form_r12_bhcabs,ttype_frm,.false.,tgt_info)
      if (active_orbs) 
     &      call set_dependency(form_r12_bhcabs,'Favg',tgt_info)
      call set_dependency(form_r12_bhcabs,op_x_inter,tgt_info)
      call set_dependency(form_r12_bhcabs,op_ff,tgt_info)
      call set_dependency(form_r12_bhcabs,op_rint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'title missing',ansatz,'BH'//approx)
      call set_rule(form_r12_bhcabs,ttype_frm,DEF_R12INTM_CABS,
     &              labels,7,1,
     &              parameters,2,tgt_info)
c dbg
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule(form_r12_bhcabs,ttype_frm,PRINT_FORMULA,
c     &              labels,2,1,
c     &              parameters,2,tgt_info)
c dbg

c     for response: may define formal B_V in this way
c      ! formal definition of B_V
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_r12_bint
c      labels(2) = op_bv_inter
c      labels(3) = op_r12
c      labels(4) = op_v
c      call add_target(form_r12_bint,ttype_frm,.false.,tgt_info)
c      call set_dependency(form_r12_bint,op_bv_inter,tgt_info)
c      call set_dependency(form_r12_bint,op_v,tgt_info)
c      call set_dependency(form_r12_bint,op_r12,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,title_r12_bint,0,'Bp')
c      call set_rule(form_r12_bint,ttype_frm,DEF_R12INTM_FORMAL,
c     &              labels,4,1,
c     &              parameters,2,tgt_info)

      ! formal definition of C
      if (active_orbs) then
        call add_target2('CINT_R12',.false.,tgt_info)
        call set_dependency(form_r12_cint,op_c_inter,tgt_info)
        call set_dependency(form_r12_cint,op_r12,tgt_info)
        call set_dependency(form_r12_cint,op_ham,tgt_info)
        call set_rule2('CINT_R12',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('CINT_R12',EXPAND_OP_PRODUCT,'LABEL',1,
     &       tgt_info,val_label=(/'CINT_R12'/))
        call set_arg('CINT_R12',EXPAND_OP_PRODUCT,'OP_RES',1,
     &       tgt_info,val_label=(/'C-INT'/))
        call set_arg('CINT_R12',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &       tgt_info,
     &       val_label=(/'C-INT','H    ',
     &                   'R12  ','C-INT'/))
        call set_arg('CINT_R12',EXPAND_OP_PRODUCT,'N_DESCR',1,
     &       tgt_info,val_int=(/2/))
        call set_arg('CINT_R12',EXPAND_OP_PRODUCT,'DESCR',2,
     &       tgt_info,
     &       val_label=(/'2,3,,X','2,,[VP],X'/))
        call set_arg('CINT_R12',EXPAND_OP_PRODUCT,'IDX_SV',4,
     &       tgt_info,val_int=(/1,2,3,1/))
      else
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
      end if
c        call set_rule2('CINT_R12',REPLACE,tgt_info)
c        call set_dependency('CINT_R12','R12-INT',tgt_info)
c        call set_arg('CINT_R12',REPLACE,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'CINT_R12'/))
c        call set_arg('CINT_R12',REPLACE,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'CINT_R12'/))
c        call set_arg('CINT_R12',REPLACE,'OP_LIST',4,tgt_info,
c     &       val_label=(/'R12      ','R12-INT  ',
c     &                   'R12^+    ','R12-INT^+'/))
c
c dbg
c      call set_rule2('CINT_R12',PRINT_FORMULA,tgt_info)
c      call set_arg('CINT_R12',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'CINT_R12'/))
c dbgend

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
      ! adding terms to C_CABS
      if (active_orbs) labels(4) = 'Favg'
      if (.not.active_orbs) then
       call add_target(form_r12_ccabs,ttype_frm,.false.,tgt_info)
      else
       call add_target2(form_r12_ccabs,.false.,tgt_info)
       call set_dependency(form_r12_ccabs,'Favg',tgt_info)
      end if
      call set_dependency(form_r12_ccabs,op_c_inter,tgt_info)
      call set_dependency(form_r12_ccabs,op_rint,tgt_info)
      call set_dependency(form_r12_ccabs,op_ham,tgt_info)
      if (active_orbs) then
c        call set_rule2('CINT_R12_CABS',EXPAND_OP_PRODUCT,tgt_info)
c        call set_arg('CINT_R12_CABS',EXPAND_OP_PRODUCT,'LABEL',1,
c     &       tgt_info,val_label=(/'CINT_R12_CABS'/))
c        call set_arg('CINT_R12_CABS',EXPAND_OP_PRODUCT,'OP_RES',1,
c     &       tgt_info,val_label=(/'C-INT'/))
c        call set_arg('CINT_R12_CABS',EXPAND_OP_PRODUCT,'OPERATORS',6,
c     &       tgt_info,
c     &       val_label=(/'C-INT','C0^+ ','H    ',
c     &                   'R12  ','C0   ','C-INT'/))
c        call set_arg('CINT_R12_CABS',EXPAND_OP_PRODUCT,'N_DESCR',1,
c     &       tgt_info,val_int=(/4/))
c        call set_arg('CINT_R12_CABS',EXPAND_OP_PRODUCT,'DESCR',4,
c     &       tgt_info,
c     &       val_label=(/'3,4,,X','3,,[VP]V,XV',
c     &                   '2,3,,V','3,5,,V'/))
c        call set_arg('CINT_R12_CABS',EXPAND_OP_PRODUCT,'IDX_SV',6,
c     &       tgt_info,val_int=(/1,2,3,4,5,1/))
c        call set_arg('CINT_R12_CABS',EXPAND_OP_PRODUCT,'N_AVOID',1,
c     &       tgt_info,val_int=(/2/))
c        call set_arg('CINT_R12_CABS',EXPAND_OP_PRODUCT,'AVOID',4,
c     &       tgt_info,val_int=(/2,6,1,5/))
        call set_rule2('CINT_R12_CABS',REPLACE,tgt_info)
        call set_dependency('CINT_R12_CABS','R12-INT',tgt_info)
        call set_dependency('CINT_R12_CABS','CINT_R12',tgt_info)
        call set_arg('CINT_R12_CABS',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'CINT_R12_CABS'/))
        call set_arg('CINT_R12_CABS',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'CINT_R12'/))
        call set_arg('CINT_R12_CABS',REPLACE,'OP_LIST',6,tgt_info,
     &       val_label=(/'R12      ','R12-INT  ',
     &                   'R12^+    ','R12-INT^+',
     &                   'H        ','Favg     '/))
      else
        call form_parameters(-1,
     &       parameters,2,title_r12_ccabs,ansatz,'C '//approx)
        call set_rule(form_r12_ccabs,ttype_frm,DEF_R12INTM_CABS,
     &                labels,nlab,1,
     &                parameters,2,tgt_info)
      end if
c dbg
c      call set_rule2('CINT_R12_CABS',PRINT_FORMULA,tgt_info)
c      call set_arg('CINT_R12_CABS',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'CINT_R12_CABS'/))
c dbgend

      ! formal definition of C1
      call add_target2('C1_formal',.false.,tgt_info)
      call set_dependency('C1_formal','C1',tgt_info)
      call set_dependency('C1_formal',op_ham,tgt_info)
      call set_dependency('C1_formal',op_r12,tgt_info)
      call set_rule2('C1_formal',DEF_R12INTM_FORMAL,tgt_info)
      call set_arg('C1_formal',DEF_R12INTM_FORMAL,'LABEL',1,tgt_info,
     &             val_label=(/'C1_formal'/))
      call set_arg('C1_formal',DEF_R12INTM_FORMAL,
     &                                              'INTERM',1,tgt_info,
     &             val_label=(/'C1'/))
      call set_arg('C1_formal',DEF_R12INTM_FORMAL,
     &                                           'OPERATORS',2,tgt_info,
     &             val_label=(/op_r12,op_ham/))
      call set_arg('C1_formal',DEF_R12INTM_FORMAL,
     &                                           'ANSATZ',1,tgt_info,
     &             val_int=(/ansatz/))
      call set_arg('C1_formal',DEF_R12INTM_FORMAL,
     &                                           'MODE',1,tgt_info,
     &             val_str='C1')
      call set_arg('C1_formal',DEF_R12INTM_FORMAL,
     &                                           'TITLE',1,tgt_info,
     &            val_str='C1 intermediate, formal definition')

      call add_target2('C1_CABS',.false.,tgt_info)
      call set_dependency('C1_CABS',op_ham,tgt_info)
      call set_dependency('C1_CABS','C1_formal',tgt_info)
      call set_dependency('C1_CABS',op_rint,tgt_info)
      if (active_orbs)
     &    call set_dependency('C1_CABS','Favg',tgt_info)
      call set_rule2('C1_CABS',REPLACE,tgt_info)
      call set_arg('C1_CABS',REPLACE,'LABEL_RES',1,tgt_info,
     &            val_label=(/'C1_CABS'/))
      call set_arg('C1_CABS',REPLACE,'LABEL_IN',1,tgt_info,
     &            val_label=(/'C1_formal'/))
      if (.not.active_orbs)
     &   call set_arg('C1_CABS',REPLACE,'OP_LIST',2,tgt_info,
     &            val_label=(/op_r12,op_rint/))
      if (active_orbs)
     &   call set_arg('C1_CABS',REPLACE,'OP_LIST',4,tgt_info,
     &            val_label=(/op_r12,op_rint,op_ham,'Favg'/))
      call set_arg('C1_CABS',REPLACE,'TITLE',1,tgt_info,
     &            val_str='C1 intermediate, for evaluation')
      ! for some exceptional cases (no CABS) make sure
      ! that R12 is removed 
      call set_rule2('C1_CABS',INVARIANT,tgt_info)
      call set_arg('C1_CABS',INVARIANT,'LABEL_RES',1,tgt_info,
     &            val_label=(/'C1_CABS'/))
      call set_arg('C1_CABS',INVARIANT,'LABEL_IN',1,tgt_info,
     &            val_label=(/'C1_CABS'/))
      call set_arg('C1_CABS',INVARIANT,'OP_RES',1,tgt_info,
     &            val_label=(/'C1'/))
      call set_arg('C1_CABS',INVARIANT,'OPERATORS',1,tgt_info,
     &            val_label=(/op_r12/))
      call set_arg('C1_CABS',INVARIANT,'TITLE',1,tgt_info,
     &            val_str='C1 intermediate, for evaluation')

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
      labels(4) = op_rint ! dummy, unused
      labels(5) = op_gr_x
      labels(6) = op_v_x
      labels(7) = op_ffg
      call add_target(form_r12_pcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_pcabs,op_p_inter,tgt_info)
      call set_dependency(form_r12_pcabs,op_ffg,tgt_info)
      call set_dependency(form_r12_pcabs,op_gr_x,tgt_info)
c      call set_dependency(form_r12_pcabs,op_g_x,tgt_info)
      call set_dependency(form_r12_pcabs,op_rint,tgt_info)
      call set_dependency(form_r12_pcabs,op_v_x,tgt_info)
      approx(12:12) = 'S'       ! set symmetrization flag
      call form_parameters(-1,
     &       parameters,2,title_r12_pcabs,ansatz,'P '//approx)
      approx(12:12) = ' '       ! unset flag
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
c      labels(3) = op_gr_x//'^+'
c      labels(4) = op_gr//'^+'
      labels(5) = op_v_x//'^+'
      labels(6) = op_v_inter//'^+'
      if (frozen) then
        nint = 4
        call set_dependency(form_r12_pcabs,'G.R-Ccore',tgt_info)
        call set_dependency(form_r12_pcabs,'V-Ccore',tgt_info)
c        labels(7) = op_gr_x//'^+'
c        labels(8) = 'G.R-Ccore^+'
        labels(7) = op_gr_x
        labels(8) = 'G.R-Ccore'
        labels(9) = op_v_x//'^+'
        labels(10) = 'V-Ccore^+'        
      end if
      call form_parameters(-1,
     &     parameters,2,title_r12_pcabs,nint,'---')
      call set_rule(form_r12_pcabs,ttype_frm,REPLACE,
     &              labels,2+nint*2,1,
     &              parameters,2,tgt_info)
c dbg
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule(form_r12_pcabs,ttype_frm,PRINT_FORMULA,
c     &              labels,2,1,
c     &              parameters,2,tgt_info)
c dbg

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
      labels(4) = op_g_z
      labels(5) = 'R.R-X'
c dbg
c      call add_target(form_r12_zcabs,ttype_frm,.true.,tgt_info)
      call add_target(form_r12_zcabs,ttype_frm,.false.,tgt_info)
      call set_dependency(form_r12_zcabs,op_z_inter,tgt_info)
      call set_dependency(form_r12_zcabs,'R.R-X',tgt_info)
      call set_dependency(form_r12_zcabs,op_rint,tgt_info)
      call set_dependency(form_r12_zcabs,op_g_z,tgt_info)
      call form_parameters(-1,
     &       parameters,2,title_r12_zcabs,ansatz,'Z '//approx)
      call set_rule(form_r12_zcabs,ttype_frm,DEF_R12INTM_CABS,
     &                labels,5,1,
     &                parameters,2,tgt_info)
!      call set_rule(form_r12_zcabs,ttype_frm,TEX_FORMULA,
!     &              labels,5,1,
!     &              'ZINT.tex',1,tgt_info)

      ! Formal definition of Z2
      if (.not.active_orbs) then
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z2INT_R12'
      labels(2) = 'Z2-INT'
      labels(3) = op_r12
      labels(4) = op_ham
      labels(5) = op_r12
      call add_target('Z2INT_R12',ttype_frm,.false.,tgt_info)
      call set_dependency('Z2INT_R12','Z2-INT',tgt_info)
      call set_dependency('Z2INT_R12',op_ham,tgt_info)
      call set_dependency('Z2INT_R12',op_r12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_r12_zint,0,'Z')
      call set_rule('Z2INT_R12',ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      else
        call set_z2int_formal(active_orbs,tgt_info)
      end if
c dbg
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule('Z2INT_R12',ttype_frm,PRINT_FORMULA,
c     &              labels,1,0,
c     &              parameters,2,tgt_info)
c dbg

      ! CABS approximation to Z2
c dbg
c      call add_target('Z2-INT-CABS',ttype_frm,.true.,tgt_info)
c dbg
      call add_target('Z2-INT-CABS',ttype_frm,.false.,tgt_info)
      call set_dependency('Z2-INT-CABS','Z2-INT',tgt_info)
      call set_dependency('Z2-INT-CABS','R.R-X',tgt_info)
      call set_dependency('Z2-INT-CABS',op_g_z,tgt_info)
      call set_dependency('Z2-INT-CABS',op_rint,tgt_info)
      ! define dummy operators for 'G-Z' and 'R.R-X'
      call set_rule2('Z2-INT-CABS',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('Z2-INT-CABS',DEF_OP_FROM_OCC,'LABEL',
     &              1,tgt_info,
     &              val_label=(/'G-Zdum'/) )
      call set_arg('Z2-INT-CABS',DEF_OP_FROM_OCC,'JOIN',
     &              1,tgt_info,
     &              val_int=(/2/) )
      call set_arg('Z2-INT-CABS',DEF_OP_FROM_OCC,'CORE',
     &              4,tgt_info,
     &              val_int=(/1,1,0,0/) )
      call set_arg('Z2-INT-CABS',DEF_OP_FROM_OCC,'DESCR',
     &              1,tgt_info,
     &              val_str='[HPX],[HPX],[HP],[HP]' )
      call set_rule2('Z2-INT-CABS',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('Z2-INT-CABS',DEF_OP_FROM_OCC,'LABEL',
     &              1,tgt_info,
     &              val_label=(/'R.R-Xdum'/) )
      call set_arg('Z2-INT-CABS',DEF_OP_FROM_OCC,'JOIN',
     &              1,tgt_info,
     &              val_int=(/2/) )
      call set_arg('Z2-INT-CABS',DEF_OP_FROM_OCC,'CORE',
     &              4,tgt_info,
     &              val_int=(/1,0,1,0/) )
      call set_arg('Z2-INT-CABS',DEF_OP_FROM_OCC,'DESCR',
     &              1,tgt_info,
     &              val_str='[HPX][HP],,,[HPX][HP]' )
      
      ! generate formula for dummy operators
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z2-INT-CABS'
      labels(2) = 'Z2-INT'
      labels(3) = op_rint
      labels(4) = op_g_z
      labels(5) = 'R.R-X'
      call form_parameters(-1,
     &     parameters,2,title_r12_xcabs,ansatz,'Z '//approx2)
      call set_rule2('Z2-INT-CABS',DEF_R12INTM_CABS,tgt_info)
      call set_arg('Z2-INT-CABS',DEF_R12INTM_CABS,'LABEL',
     &              1,tgt_info,
     &              val_label=(/'Z2-INT-CABS'/) )
      call set_arg('Z2-INT-CABS',DEF_R12INTM_CABS,'INTERM',
     &              1,tgt_info,
     &              val_label=(/'Z2-INT'/) )
      call set_arg('Z2-INT-CABS',DEF_R12INTM_CABS,'OPERATORS',
     &              3,tgt_info,
     &              val_label=(/op_rint,op_g_z,'R.R-X'/) )
C     &              val_label=(/op_rint,'G-Zdum','R.R-Xdum'/) )
      call set_arg('Z2-INT-CABS',DEF_R12INTM_CABS,'ANSATZ',
     &              1,tgt_info,
     &              val_int=(/ansatz/) )
      call set_arg('Z2-INT-CABS',DEF_R12INTM_CABS,'APPROX',
     &              1,tgt_info,
     &              val_str='Z '//approx2)
      call set_arg('Z2-INT-CABS',DEF_R12INTM_CABS,'TITLE',
     &              1,tgt_info,
     &              val_str='Z2 CABS ')
      
      ! replace dummy operators by actual operators and their h.c.'s

      ! remove dummy operators from list
c dbg
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule('Z2-INT-CABS',ttype_frm,PRINT_FORMULA,
c     &              labels,1,0,
c     &              parameters,2,tgt_info)
c dbg


      ! set semi-internal R12 geminals
      if (semi_r12) then
      descr = 'VX,[VH]V|VX,HH'
      call add_target2('R12si',.false.,tgt_info)
      call set_rule2('R12si',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('R12si',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'R12si'/))
      call set_arg('R12si',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('R12si',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)
      descr = 'X,V|X,H'
      call add_target2('sR12',.false.,tgt_info)
      call set_rule2('sR12',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('sR12',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'sR12'/))
      call set_arg('sR12',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('sR12',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)

      descr = ',V,V,'
      call add_target2('DENSinv',.false.,tgt_info)
      call set_rule2('DENSinv',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('DENSinv',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'DENSinv'/))
      call set_arg('DENSinv',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/2/))
      call set_arg('DENSinv',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)
      call add_target2('DEF_ME_DENSinv',.false.,tgt_info)
      call set_dependency('DEF_ME_DENSinv','DENSinv',tgt_info)
      call set_rule2('DEF_ME_DENSinv',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_DENSinv',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_DENSinv'/))
      call set_arg('DEF_ME_DENSinv',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'DENSinv   '/))
      call set_arg('DEF_ME_DENSinv',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_DENSinv',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_DENSinv',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      call set_dependency('DEF_ME_DENSinv','EVAL_D',tgt_info)
      call set_rule2('DEF_ME_DENSinv',SCALE_COPY,tgt_info)
      call set_arg('DEF_ME_DENSinv',SCALE_COPY,'LIST_RES',1,tgt_info,
     &             val_label=(/'ME_DENSinv'/))
      call set_arg('DEF_ME_DENSinv',SCALE_COPY,'LIST_INP',1,tgt_info,
     &             val_label=(/'ME_DENS   '/))
      call set_arg('DEF_ME_DENSinv',SCALE_COPY,'FAC',1,tgt_info,
     &             val_rl8=(/1d0/))
      call set_rule2('DEF_ME_DENSinv',INVERT,tgt_info)
      call set_arg('DEF_ME_DENSinv',INVERT,'LIST_INV',1,tgt_info,
     &             val_label=(/'ME_DENSinv'/))
      call set_arg('DEF_ME_DENSinv',INVERT,'LIST',1,tgt_info,
     &             val_label=(/'ME_DENSinv'/))
      call set_arg('DEF_ME_DENSinv',INVERT,'MODE',1,tgt_info,
     &             val_str='pseudoinv')

      call add_target2('sR12-INT',.false.,tgt_info)
      call set_dependency('sR12-INT','sR12',tgt_info)
      call set_dependency('sR12-INT','R12si',tgt_info)
      call set_dependency('sR12-INT','F_DENS0',tgt_info)
      call set_dependency('sR12-INT','DEF_ME_DENSinv',tgt_info)
      call set_rule2('sR12-INT',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'sR12-INT'/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.true./))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'sR12'/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &     tgt_info,
     &     val_label=(/'sR12','DENS','R12si','DENS','sR12'/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &     val_int=(/1,2,3,2,1/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'DESCR',1,tgt_info,
     &     val_label=(/'3,,VX,VH'/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'CONNECT',4,tgt_info,
     &             val_int=(/2,3,3,4/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &             val_int=(/2,4/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/-1d0/))
      call set_rule2('sR12-INT',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'sR12-INT'/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'sR12'/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'OPERATORS',7,
     &     tgt_info,
     &     val_label=(/'sR12','DENSinv','DENS','R12si',
     &                 'DENS','DENSinv','sR12'/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'IDX_SV',7,tgt_info,
     &     val_int=(/1,2,3,4,3,2,1/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'DESCR',3,tgt_info,
     &     val_label=(/'3,4,,V','4,,VX,VV','4,5,,VV'/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'AVOID',10,tgt_info,
     &             val_int=(/2,3,5,6,2,6,3,5,2,5/))
      call set_arg('sR12-INT',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/-1d0/))

c dbg
c      call set_rule2('sR12-INT',PRINT_FORMULA,tgt_info)
c      call set_arg('sR12-INT',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'sR12-INT'/))
c dbgend

      call add_target2('Rsi-INT',.false.,tgt_info)
      call set_dependency('Rsi-INT','R12si',tgt_info)
      call set_rule2('Rsi-INT',CLONE_OP,tgt_info)
      call set_arg('Rsi-INT',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Rsi-INT'/))
      call set_arg('Rsi-INT',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'R12si'/))
      call set_arg('Rsi-INT',CLONE_OP,'ADJOINT',1,tgt_info,
     &     val_log=(/.false./))
      call add_target2('ME_Rsi',.false.,tgt_info)
      call set_dependency('ME_Rsi','Rsi-INT',tgt_info)
      call set_rule2('ME_Rsi',DEF_ME_LIST,tgt_info)
      call set_arg('ME_Rsi',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Rsi '/))
      call set_arg('ME_Rsi',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'Rsi-INT'/))
      call set_arg('ME_Rsi',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('ME_Rsi',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('ME_Rsi',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      call set_rule2('ME_Rsi',IMPORT,tgt_info)
      call set_arg('ME_Rsi',IMPORT,'LIST',1,tgt_info,
     &             val_label=(/'ME_Rsi '/))
      call set_arg('ME_Rsi',IMPORT,'TYPE',1,tgt_info,
     &             val_str='F12_INT')
      call set_arg('ME_Rsi',IMPORT,'ENV',1,tgt_info,
     &             val_str=env_type)
      end if
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
      if (frozen)
     &     call set_dependency(fopt_r12_vcabs,'G-Ac-INT',tgt_info)
      if (.not.frozen.and.t1ext.eq.1)
     &     call set_dependency(fopt_r12_vcabs,'H-ext-INT',tgt_info)
      call set_dependency(fopt_r12_vcabs,mel_rint,tgt_info)      
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

      call add_target2('Vring-EVAL',.false.,tgt_info)
      call set_dependency('Vring-EVAL','Vring_CABS',tgt_info)
      call set_dependency('Vring-EVAL',mel_vp_def,tgt_info)
      call set_dependency('Vring-EVAL',mel_ham,tgt_info)
      call set_dependency('Vring-EVAL',mel_rint,tgt_info)
      call set_rule2('Vring-EVAL',OPTIMIZE,tgt_info)
      call set_arg('Vring-EVAL',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &            val_label=(/'Vring-OPT'/))
      if (vring_mode.eq.2) then
        call set_arg('Vring-EVAL',OPTIMIZE,'LABELS_IN',2,tgt_info,
     &            val_label=(/'Vring_CABS ','Vring2_CABS'/))
      else
        call set_arg('Vring-EVAL',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &            val_label=(/'Vring_CABS'/))
      end if
      call set_rule2('Vring-EVAL',EVAL,tgt_info)
      call set_arg('Vring-EVAL',EVAL,'FORM',1,tgt_info,
     &            val_label=(/'Vring-OPT'/))

      ! set X
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_r12_xcabs
      labels(2) = form_r12_xcabs
      ncat = 1
      nint = 0
c      call add_target(fopt_r12_xcabs,ttype_frm,.true.,tgt_info)
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
      if (active_orbs)
     &    call set_dependency(fopt_r12_bcabs,'Favg-INT',tgt_info)
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
      if (active_orbs)
     &    call set_dependency(fopt_r12_bhcabs,'Favg-INT',tgt_info)
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
      if (active_orbs)
     &    call set_dependency(fopt_r12_ccabs,'Favg-INT',tgt_info)
      call set_rule(fopt_r12_ccabs,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      call add_target2('C1-EVAL',.false.,tgt_info)
      call set_dependency('C1-EVAL','C1_CABS',tgt_info)
      call set_dependency('C1-EVAL','DEF-C1INT',tgt_info)
      call set_dependency('C1-EVAL',mel_ham,tgt_info)
      call set_dependency('C1-EVAL',mel_rint,tgt_info)
      call set_rule2('C1-EVAL',OPTIMIZE,tgt_info)
      call set_arg('C1-EVAL',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &            val_label=(/'C1-OPT'/))
      call set_arg('C1-EVAL',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &            val_label=(/'C1_CABS'/))
      call set_rule2('C1-EVAL',EVAL,tgt_info)
      call set_arg('C1-EVAL',EVAL,'FORM',1,tgt_info,
     &            val_label=(/'C1-OPT'/))

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
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

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

      ! set Z2 (direct evaluation)
      call add_target('Z2INT_R12_DIR',ttype_frm,.false.,tgt_info)
      call set_dependency('Z2INT_R12_DIR','Z2INT_R12',tgt_info)
      call set_dependency('Z2INT_R12_DIR',mel_ham,tgt_info)
      call set_dependency('Z2INT_R12_DIR',mel_rint,tgt_info)      

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z2INT_R12_DIR'
      labels(2) = 'Z2INT_R12'
      labels(3) = op_r12
      labels(4) = op_rint
      labels(5) = op_r12//'^+'
      labels(6) = op_rint//'^+'
      nint = 2
      call form_parameters(-1,
     &     parameters,2,'Z2 direct',nint,'---')
      call set_rule('Z2INT_R12_DIR',ttype_frm,REPLACE,
     &     labels,2+nint*2,1,
     &     parameters,2,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z2INT_R12_DIR'
      labels(2) = 'Z2INT_R12_DIR'
      ncat = 1
      nint = 0
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('Z2INT_R12_DIR',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z2INT_R12_DIR'
      call set_rule('Z2INT_R12_DIR',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      ! set Z2 (reformulated evaluation)
      call add_target('Z2INT_R12_REF',ttype_frm,.false.,tgt_info)
      call set_dependency('Z2INT_R12_REF','Z2INT_R12',tgt_info)
      call set_dependency('Z2INT_R12_REF','Z2-INT-CABS',tgt_info)
      call set_dependency('Z2INT_R12_REF',mel_ham,tgt_info)
      call set_dependency('Z2INT_R12_REF',mel_rint,tgt_info)
      call set_dependency('Z2INT_R12_REF','R.R-INTX',tgt_info)
      call set_dependency('Z2INT_R12_REF',mel_gintz,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z2INT_R12_REF'
      labels(2) = 'Z2-INT-CABS'
      ncat = 1
      nint = 0
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('Z2INT_R12_REF',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z2INT_R12_REF'
      call set_rule('Z2INT_R12_REF',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      call add_target('Z2INT_R12_EVAL',ttype_frm,.false.,tgt_info)
      call set_dependency('Z2INT_R12_EVAL','DEF-Z2LIST',tgt_info)
      if (approx2(14:17).eq.'DRCT') then
        call set_dependency('Z2INT_R12_EVAL','Z2INT_R12_DIR',tgt_info)
      else
        call set_dependency('Z2INT_R12_EVAL','Z2INT_R12_REF',tgt_info)
      end if

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

      ! old special two-electron integral list
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

      ! for MR: import the averaged FOCK (as used for integrals)
      call add_target2('Favg-INT',.false.,tgt_info)
      call set_dependency('Favg-INT','Favg',tgt_info)
      ! (a) define
      call set_rule2('Favg-INT',DEF_ME_LIST,tgt_info)
      call set_arg('Favg-INT',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'Favg-INT'/))
      call set_arg('Favg-INT',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'Favg'/))
      call set_arg('Favg-INT',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('Favg-INT',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('Favg-INT',DEF_ME_LIST,'ABSYM',1,tgt_info,
     &     val_int=(/msc/))
      ! (b) import
      call set_rule2('Favg-INT',IMPORT,tgt_info)
      call set_arg('Favg-INT',IMPORT,'LIST',1,tgt_info,
     &     val_label=(/'Favg-INT'/))
      call set_arg('Favg-INT',IMPORT,'TYPE',1,tgt_info,
     &     val_str='F_INT')
      call set_arg('Favg-INT',IMPORT,'ENV',1,tgt_info,
     &     val_str=env_type)

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

      ! special R12^2 integrals
      call add_target('FF-X-INT',ttype_opme,.false.,tgt_info)
      call set_dependency('FF-X-INT','FF-X',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FF-X-INT'
      labels(2) = 'FF-X'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('FF-X-INT',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FF-X-INT'
      call import_parameters(-1,parameters,'FF_INT',env_type)
      call set_rule('FF-X-INT',ttype_opme,IMPORT,
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

      call add_target('G.R-X-INT',ttype_opme,.false.,tgt_info)
      call set_dependency('G.R-X-INT','G.R-X',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'G.R-X-INT'
      labels(2) = 'G.R-X'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('G.R-X-INT',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'G.R-X-INT'
      call import_parameters(-1,parameters,'FG_INT',env_type)
      call set_rule('G.R-X-INT',ttype_opme,IMPORT,
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

      ! V'-list
      call add_target(mel_vp_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_vp_def,op_vp_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_vp_inter
      labels(2) = op_vp_inter
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_vp_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      if(.not.active_orbs) then
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'VRnonsym'
      labels(2) = 'VR2'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_vp_def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      end if

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
      
      ! C1-list
      call add_target('DEF-C1INT',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF-C1INT','C1',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'C1INT'
      labels(2) = 'C1'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF-C1INT',ttype_opme,DEF_ME_LIST,
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


      call add_target('Vpx-INTER',ttype_opme,.false.,tgt_info)
      call set_dependency('Vpx-INTER','Vpx',tgt_info)
      call set_dependency('Vpx-INTER','Vpx_CABS',tgt_info)
      call set_dependency('Vpx-INTER','G.R-X-INT',tgt_info)
      call set_dependency('Vpx-INTER',mel_gintx,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Vpx-INTER'
      labels(2) = 'Vpx'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('Vpx-INTER',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1) = 'Vpx_OPT'
      labels(2) = 'Vpx_CABS'
      call opt_parameters(-1,parameters,1,0)
      call set_rule('Vpx-INTER',ttype_frm,OPTIMIZE,
     &     labels,2,1,
     &     parameters,1,tgt_info)
      labels(1) = 'Vpx_OPT'
      call set_rule('Vpx-INTER',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      call add_target(mel_p_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_p_def,op_p_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_p_inter
      labels(2) = op_p_inter
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_p_def,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
c      ! work-around for frozen core
c      if (.not.pz_eval) then
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = mel_p_def
c        call import_parameters(-1,parameters,'P_LIST',env_type)
c        call set_rule(mel_p_def,ttype_opme,IMPORT,
c     &                labels,1,1,
c     &                parameters,1,tgt_info)
c      end if

      ! Z-list
      call add_target(mel_z_def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_z_def,op_z_inter,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_z_inter
      labels(2) = op_z_inter
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_z_def,ttype_opme,DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)
c      ! work-around for frozen core
c      if (.not.pz_eval) then
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = mel_z_def
c        call import_parameters(-1,parameters,'Z_LIST',env_type)
c        call set_rule(mel_z_def,ttype_opme,IMPORT,
c     &                labels,1,1,
c     &                parameters,1,tgt_info)
c      end if
      

      ! Z2-list
      call add_target('DEF-Z2LIST',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF-Z2LIST','Z2-INT',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'Z2LIST'
      labels(2) = 'Z2-INT'
      call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
      call set_rule('DEF-Z2LIST',ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      ! test
      call add_target(eval_r12_inter,ttype_gen,.false.,tgt_info)
      call set_dependency(eval_r12_inter,mel_ham,tgt_info)
      call set_dependency(eval_r12_inter,mel_rint,tgt_info)
c      call set_dependency(eval_r12_inter,mel_gintx,tgt_info)
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

      call add_target('EVAL_PZ',ttype_gen,pz_eval,tgt_info)
      call set_dependency('EVAL_PZ',mel_p_def,tgt_info)
      call set_dependency('EVAL_PZ',mel_z_def,tgt_info)
c dbg
      if (.true.) then
c dbg
c      if (pz_eval) then
        call set_dependency('EVAL_PZ',mel_rint,tgt_info)
        call set_dependency('EVAL_PZ','R.R-INTX',tgt_info)
        call set_dependency('EVAL_PZ',mel_ffg,tgt_info)
        call set_dependency('EVAL_PZ',mel_gr,tgt_info)
        call set_dependency('EVAL_PZ',fopt_r12_pcabs,tgt_info)
        call set_dependency('EVAL_PZ',eval_r12_inter,tgt_info)
        if (frozen) then
          call set_dependency('EVAL_PZ','V-Ccore-OPT',tgt_info)

          labels(1) = 'V-Ccore-OPT'
          call set_rule('EVAL_PZ',ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)
        end if
        if (approx(14:17).eq.'DRCT') then
          call set_dependency('EVAL_PZ','ZINT_R12_DIR',tgt_info)
        else
          call set_dependency('EVAL_PZ',fopt_r12_zcabs,tgt_info)
          call set_dependency('EVAL_PZ',mel_gintz,tgt_info)
        end if

        labels(1) = fopt_r12_pcabs
        call set_rule('EVAL_PZ',ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)
          
        if (approx(14:17).eq.'DRCT') then
          labels(1) = 'ZINT_R12_DIR'
        else
          labels(1) = fopt_r12_zcabs
        end if
        call set_rule('EVAL_PZ',ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)

c      end if

c dbg
      else
      if (pz_eval) then
        call set_dependency('EVAL_PZ',mel_rint,tgt_info)
        call set_dependency('EVAL_PZ',mel_ff,tgt_info)
        if (approx(14:17).eq.'DRCT') then
          call set_dependency('EVAL_PZ','ZINT_R12_DIR',tgt_info)
          call set_dependency('EVAL_PZ',mel_ham,tgt_info)
        else
          call set_dependency('EVAL_PZ',fopt_r12_zcabs,tgt_info)
          call set_dependency('EVAL_PZ',mel_gintz,tgt_info)
        end if

        if (approx(14:17).eq.'DRCT') then
          labels(1) = 'ZINT_R12_DIR'
        else
          labels(1) = fopt_r12_zcabs
        end if
        call set_rule('EVAL_PZ',ttype_opme,EVAL,
     &         labels,1,0,
     &         parameters,0,tgt_info)

      end if

      end if
c dbg

      return

      contains

      subroutine set_gxx(occ_def)

      implicit none

      integer, intent(out) ::
     &     occ_def(ngastp,2,*)

      integer ::
     &     idxc, idxa, idx_end, occ_c(ngastp), occ_a(ngastp)

      occ_def(1:ngastp,1:2,1:4*4) = 0
      ndef = 0
      idx_end = 4
      if (active_orbs) idx_end = 8
      do idxc = 1, idx_end
        occ_c = 0
        if (idxc.eq.1) occ_c(IHOLE) = 2
        if (idxc.eq.2) occ_c(IHOLE) = 1
        if (idxc.eq.2) occ_c(IPART) = 1
        if (idxc.eq.3) occ_c(IHOLE) = 1
        if (idxc.eq.3) occ_c(IEXTR) = 1
        if (idxc.eq.4) occ_c(IPART) = 2
        if (idxc.eq.5) occ_c(IHOLE) = 1
        if (idxc.eq.5) occ_c(IVALE) = 1
        if (idxc.eq.6) occ_c(IVALE) = 2
        if (idxc.eq.7) occ_c(IVALE) = 1
        if (idxc.eq.7) occ_c(IPART) = 1
        if (idxc.eq.8) occ_c(IVALE) = 1
        if (idxc.eq.8) occ_c(IEXTR) = 1
        do idxa = 1, idx_end
          occ_a = 0
          if (idxa.eq.1) occ_a(IHOLE) = 2
          if (idxa.eq.2) occ_a(IHOLE) = 1
          if (idxa.eq.2) occ_a(IPART) = 1
          if (idxa.eq.3) occ_a(IHOLE) = 1
          if (idxa.eq.3) occ_a(IEXTR) = 1
          if (idxa.eq.4) occ_a(IPART) = 2
          if (idxa.eq.5) occ_a(IHOLE) = 1
          if (idxa.eq.5) occ_a(IVALE) = 1
          if (idxa.eq.6) occ_a(IVALE) = 2
          if (idxa.eq.7) occ_a(IVALE) = 1
          if (idxa.eq.7) occ_a(IPART) = 1
          if (idxa.eq.8) occ_a(IVALE) = 1
          if (idxa.eq.8) occ_a(IEXTR) = 1

          ndef = ndef+1
          occ_def(1:ngastp,1,ndef) = occ_c(1:ngastp)
          occ_def(1:ngastp,2,ndef) = occ_a(1:ngastp)
          
        end do
      end do

      end subroutine set_gxx

      subroutine set_g_z(ndef,occ_def)

      implicit none

      integer, intent(out) ::
     &     ndef,occ_def(ngastp,2,*)

      ndef = 25
      occ_def(1:ngastp,1:2,1:2*ndef) = 0
      ! set (H,[HPX])(H,H) part
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
      
      ! set (H,[HPX])(H,P) part
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

      ! set (H,[HPX])(P,H) part
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

      ! set (P,[HPX])(H,H)
      ! 10
      occ_def(IPART,1,19) = 1
      occ_def(IHOLE,2,19) = 1
      occ_def(IHOLE,1,20) = 1
      occ_def(IHOLE,2,20) = 1
      ! 11
      occ_def(IPART,1,21) = 1
      occ_def(IPART,2,21) = 1
      occ_def(IHOLE,1,22) = 1
      occ_def(IHOLE,2,22) = 1
      ! 12
      occ_def(IPART,1,23) = 1
      occ_def(IEXTR,2,23) = 1
      occ_def(IHOLE,1,24) = 1
      occ_def(IHOLE,2,24) = 1
      ! set (X,[HP])(H,H)
      ! 13
      occ_def(IEXTR,1,25) = 1
      occ_def(IHOLE,2,25) = 1
      occ_def(IHOLE,1,26) = 1
      occ_def(IHOLE,2,26) = 1
      ! 14
      occ_def(IEXTR,1,27) = 1
      occ_def(IPART,2,27) = 1
      occ_def(IHOLE,1,28) = 1
      occ_def(IHOLE,2,28) = 1

      ! set (P,[HPX])(H,P)
      ! 15
      occ_def(IPART,1,29) = 1
      occ_def(IHOLE,2,29) = 1
      occ_def(IHOLE,1,30) = 1
      occ_def(IPART,2,30) = 1
      ! 16
      occ_def(IPART,1,31) = 1
      occ_def(IPART,2,31) = 1
      occ_def(IHOLE,1,32) = 1
      occ_def(IPART,2,32) = 1
      ! 17
      occ_def(IPART,1,33) = 1
      occ_def(IEXTR,2,33) = 1
      occ_def(IHOLE,1,34) = 1
      occ_def(IPART,2,34) = 1
      ! set (X,[HP])(H,P)
      ! 18
      occ_def(IEXTR,1,35) = 1
      occ_def(IHOLE,2,35) = 1
      occ_def(IHOLE,1,36) = 1
      occ_def(IPART,2,36) = 1
      ! 19
      occ_def(IEXTR,1,37) = 1
      occ_def(IPART,2,37) = 1
      occ_def(IHOLE,1,38) = 1
      occ_def(IPART,2,38) = 1

      ! set (P,[HP])(H,X)
      ! 20
      occ_def(IPART,1,39) = 1
      occ_def(IHOLE,2,39) = 1
      occ_def(IHOLE,1,40) = 1
      occ_def(IEXTR,2,40) = 1
      ! 21
      occ_def(IPART,1,41) = 1
      occ_def(IPART,2,41) = 1
      occ_def(IHOLE,1,42) = 1
      occ_def(IEXTR,2,42) = 1

      ! set (X,X)(H,H)
      ! 22
      occ_def(IEXTR,1,43) = 1
      occ_def(IEXTR,2,43) = 1
      occ_def(IHOLE,1,44) = 1
      occ_def(IHOLE,2,44) = 1

      ! set (X,X)(H,P)
      ! 23
      occ_def(IEXTR,1,45) = 1
      occ_def(IEXTR,2,45) = 1
      occ_def(IHOLE,1,46) = 1
      occ_def(IPART,2,46) = 1

      ! set (X,[HP])(H,X)
      ! 24
      occ_def(IEXTR,1,47) = 1
      occ_def(IHOLE,2,47) = 1
      occ_def(IHOLE,1,48) = 1
      occ_def(IEXTR,2,48) = 1
      ! 25
      occ_def(IEXTR,1,49) = 1
      occ_def(IPART,2,49) = 1
      occ_def(IHOLE,1,50) = 1
      occ_def(IEXTR,2,50) = 1

      end subroutine

      subroutine set_g_z_old(ndef,occ_def)

      implicit none

      integer, intent(out) ::
     &     ndef,occ_def(ngastp,2,*)

      ndef = 26
        occ_def(1:ngastp,1:2,1:ndef) = 0
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

      end subroutine

      end


