*----------------------------------------------------------------------*
      subroutine set_cc_pt_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
*     set targets needed specifically in CC calculations
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
      character(*), intent(in) ::
     &     env_type

      integer ::
     &     max_rank, min_rank_pt, max_rank_pt, max_extern,
     &     min_rank_tp, min_rank_tpp, ntp_min, ntp_max,
     &     ntpp_min, ntpp_max,
     &     ansatz,
     &     isim, ncat, nint, icnt, r12op,
     &     isym, ms, msc, sym_arr(8),
     &     occ_def(ngastp,2,60), ndef
      logical ::
     &     needed, explicit, r12fix, set_tp, set_tpp, gbc4pt,
     &     R2_coupling, set_R2_R2, set_Rhxhh
      character(len=8) ::
     &     mode
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(2)
      character(12) ::
     &     approx

      msc = +1 ! assuming closed-shell

      icnt = is_keyword_set('method.CCPT')
      if (icnt.eq.0) return

      explicit = is_keyword_set('method.R12').gt.0
      set_tp = .false.
      set_tpp = .false.

      if (iprlvl.gt.0.and..not.explicit)
     &     write(luout,*) 'setting special targets for CCPT ...'
      if (iprlvl.gt.0.and.explicit)
     &     write(luout,*) 'setting special targets for CCPT (R12) ...'

      call get_argument_value('method.CC','maxexc',ival=max_rank)
      min_rank_pt = max_rank+1
      call get_argument_value('method.CCPT','maxexc',ival=max_rank_pt)      
      call get_argument_value('method.CCPT','extern',ival=max_extern)      
      call get_argument_value('method.CCPT','GBC',lval=gbc4pt)
      call get_argument_value('method.CCPT','R2_coupling',
     &     lval=R2_coupling)
      call get_argument_value('method.CCPT','R2R2',
     &     lval=set_R2_R2)
      call get_argument_value('method.CCPT','hh_scatter',
     &     lval=set_Rhxhh)

      if (explicit) then
        approx(1:12) = ' '

        call get_argument_value('method.R12','ansatz',ival=ansatz)      
        call get_argument_value('method.R12','approx',str=approx)
        call get_argument_value('method.R12','fixed',lval=r12fix)      
        call get_argument_value('method.R12','r12op',ival=r12op)      
        call get_argument_value('method.R12','min_tp',ival=min_rank_tp)
        call get_argument_value('method.R12','min_tpp',
     &                                               ival=min_rank_tpp)
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
          ntp_min=max(min_rank_pt-1,min_rank_tp)
          ntp_max=max_rank_pt-1
c          n_pp=1
        case(2)
          ! T'' operators (doubly p-connected to R12)
          set_tpp = .true.
          ntpp_min=max(min_rank_pt,min_rank_tpp)
          ntpp_max=max_rank_pt
c          n_pp=2
        case(3,4)
          ! T' + T'' operators
          set_tp = .true.
          ntp_min=max(min_rank_pt-1,min_rank_tp)
          ntp_max=max_rank_pt-1
c          n_pp=1
          set_tpp = .true.
          ntpp_min=max(min_rank_pt,min_rank_tpp)
          ntpp_max=max_rank_pt
c          n_pp=2
        end select

      else
        ansatz = 0
      end if
*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! Lagrange functional 
      call add_target(op_dlpt,ttype_op,.false.,tgt_info)
      call set_rule(op_dlpt,ttype_op,DEF_SCALAR,
     &              op_dlpt,1,1,
     &              parameters,0,tgt_info)
      
      ! Energy 
      call add_target(op_dept,ttype_op,.false.,tgt_info)
      call set_rule(op_dept,ttype_op,DEF_SCALAR,
     &              op_dept,1,1,
     &              parameters,0,tgt_info)

c      ! formal R12(ic|ab) block for testing
c      call add_target('R12VV',ttype_op,.false.,tgt_info)
c      occ_def = 0
c      ndef = 1
c      occ_def(IPART,1,1) = 2
c      occ_def(IHOLE,2,1) = 1
c      occ_def(IPART,2,1) = 1
c      call op_from_occ_parameters(-1,parameters,2,
c     &                         occ_def,ndef,1,(/.true.,.true./),ndef)
c      call set_rule('R12VV',ttype_op,DEF_OP_FROM_OCC,
c     &     'R12VV',1,1,
c     &     parameters,2,tgt_info)

      ! R12<i beta|kl> block for testing
      call add_target('R12hxhh',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 1
      occ_def(IHOLE,1,1) = 1
      occ_def(IEXTR,1,1) = 1
      occ_def(IHOLE,2,1) = 2
      call op_from_occ_parameters(-1,parameters,2,
     &                         occ_def,ndef,1,(/.true.,.true./),ndef)
      call set_rule('R12hxhh',ttype_op,DEF_OP_FROM_OCC,
     &     'R12hxhh',1,1,
     &     parameters,2,tgt_info)

      ! T(pt) operator
      call add_target(op_tpt,ttype_op,.false.,tgt_info)
      if (max_extern.eq.0) then
        call xop_parameters(-1,parameters,
     &       .false.,min_rank_pt,max_rank_pt,0,1)
        call set_rule(op_tpt,ttype_op,DEF_EXCITATION,
     &       op_tpt,1,1,
     &       parameters,1,tgt_info)
      else ! try n-externals
        if (max_rank_pt-min_rank_pt.ne.0)
     &       call quit(1,'set_cc_pt_targets','not allowed!')
        occ_def = 0
        ! 0-ext
        occ_def(IPART,1,1) = max_rank_pt
        occ_def(IHOLE,2,1) = max_rank_pt
        ! 1-ext
        occ_def(IPART,1,2) = max_rank_pt-1
        occ_def(IEXTR,1,2) = 1
        occ_def(IHOLE,2,2) = max_rank_pt
        ndef = 2
        if (max_extern.ge.2) then
          ! 2-ext
          occ_def(IPART,1,3) = max_rank_pt-2
          occ_def(IEXTR,1,3) = 2
          occ_def(IHOLE,2,3) = max_rank_pt
          ndef = 3
        end if
        if (max_extern.gt.2)
     &       call quit(1,'set_cc_pt_targets','extern>2 not allowed!')
        call op_from_occ_parameters(-1,parameters,2,
     &                         occ_def,ndef,1,(/.true.,.true./),ndef)
        call set_rule(op_tpt,ttype_op,DEF_OP_FROM_OCC,
     &                op_tpt,1,1,
     &                parameters,2,tgt_info)
      end if
      ! T'(pt) operator
      if (set_tp) then
        call add_target(op_cex_pt,ttype_op,.false.,tgt_info)
        call xop_parameters(-1,parameters,
     &       .false.,ntp_min,ntp_max,0,ntp_max+2)
        call set_rule(op_cex_pt,ttype_op,DEF_EXCITATION,
     &                op_cex_pt,1,1,
     &                parameters,1,tgt_info)
      end if
      ! T''(pt) operator
      if (set_tpp) then
        call add_target(op_cexx_pt,ttype_op,.false.,tgt_info)
        call xop_parameters(-1,parameters,
     &       .false.,ntpp_min,ntpp_max,0,ntp_max+2)
        call set_rule(op_cexx_pt,ttype_op,DEF_EXCITATION,
     &                op_cexx_pt,1,1,
     &                parameters,1,tgt_info)
      end if
      ! Diagonal
      call add_target(op_dia_pt,ttype_op,.false.,tgt_info)
      call set_dependency(op_dia_pt,op_tpt,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tpt,.false.)
      call set_rule(op_dia_pt,ttype_op,CLONE_OP,
     &              op_dia_pt,1,1,
     &              parameters,1,tgt_info)
      ! H0.T(pt)
      call add_target(op_h0_r,ttype_op,.false.,tgt_info)
      call set_dependency(op_h0_r,op_tpt,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tpt,.false.)
      call set_rule(op_h0_r,ttype_op,CLONE_OP,
     &              op_h0_r,1,1,
     &              parameters,1,tgt_info)
      ! RHS
      call add_target(op_etapt,ttype_op,.false.,tgt_info)
      call set_dependency(op_etapt,op_tpt,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tpt,.false.)
      call set_rule(op_etapt,ttype_op,CLONE_OP,
     &              op_etapt,1,1,
     &              parameters,1,tgt_info)

      ! V-X operator
      call add_target('Vpx',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 2
      occ_def(IEXTR,1,1) = 1
      occ_def(IPART,2,1) = 1
      occ_def(IHOLE,1,2) = 1
      occ_def(IEXTR,1,2) = 1
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/.true.,.true./),ndef)
      call set_rule('Vpx',ttype_op,DEF_OP_FROM_OCC,
     &              'Vpx',1,1,
     &              parameters,2,tgt_info)
      
      call add_target('G.R-X',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 1
      occ_def(IHOLE,1,1) = 1
      occ_def(IEXTR,1,1) = 1
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,2,(/.true.,.true./),ndef)
      call set_rule('G.R-X',ttype_op,DEF_OP_FROM_OCC,
     &              'G.R-X',1,1,
     &              parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_ptdl0
      labels(2) = op_dlpt
      labels(3) = op_ham
      labels(4) = op_top
      labels(5) = op_tpt
      labels(6) = op_top
      labels(7) = op_tpt
      nint = 7
      call add_target(form_ptdl0,ttype_frm,.false.,tgt_info)
      if (explicit) then
        labels(8) = op_r12
        nint = 8
        if (.not.r12fix) then
          labels(9) = op_c12
          labels(10) = op_c12
          nint = 10
        else if (r12op.gt.0) then
          if (r12op.ne.2) then
            call set_dependency(form_ptdl0,op_cex,tgt_info)
            call set_dependency(form_ptdl0,op_cex_pt,tgt_info)
            labels(nint+1) = op_cex
            labels(nint+2) = op_cex_pt
            nint = nint+2
          end if
          if (r12op.ge.2) then
            call set_dependency(form_ptdl0,op_cexx,tgt_info)
            call set_dependency(form_ptdl0,op_cexx_pt,tgt_info)
            labels(nint+1) = op_cexx
            labels(nint+2) = op_cexx_pt
            nint = nint+2
          end if
        end if
      end if
      call set_dependency(form_ptdl0,op_dlpt,tgt_info)
      call set_dependency(form_ptdl0,op_ham,tgt_info)
c      call set_dependency(form_ptdl0,op_tbar,tgt_info)
      call set_dependency(form_ptdl0,op_top,tgt_info)
      call set_dependency(form_ptdl0,op_tpt,tgt_info)
      mode = '--------'
      if (max_extern.gt.0) mode(1:4) = 'EXT '
      if (gbc4pt)          mode(4:4) =    '0'
      if (.not.R2_coupling)  mode(5:8) = 'NOR2'

      if (r12op.gt.0.and.set_Rhxhh) then 
        mode(1:3) = 'hhs'
        call set_dependency(form_ptdl0,'R12hxhh',tgt_info)
        labels(nint+1) = 'R12hxhh'
        nint = nint+1
      end if
      call form_parameters(-1,
     &     parameters,2,title_ptdl0,ansatz,mode)
      call set_rule(form_ptdl0,ttype_frm,DEF_CCPT_LAGRANGIAN,
     &              labels,nint,1,
     &              parameters,2,tgt_info)
      if (explicit) then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_ptdl0 ! output formula (itself)
        labels(2) = form_ptdl0 ! input formula
        if (r12fix.or.r12op.gt.0) then
          call set_dependency(form_ptdl0,form_r12_bhint,tgt_info)
          call set_dependency(form_ptdl0,form_r12_xhint,tgt_info)
          labels(3) = form_r12_vint ! the intermediates to be factored
          labels(4) = form_r12_vint//'^+'
          labels(5) = form_r12_bint
          labels(6) = form_r12_bhint
          labels(7) = form_r12_xint
          nint = 5
          if (r12op.gt.0) then
            call set_dependency(form_ptdl0,'Vpx_formal',tgt_info)
            labels(8) = 'Vpx_formal'
            labels(9) = 'Vpx_formal^+'
            labels(10) = form_r12_xhint
            nint = 8
          end if
        else
          labels(3) = form_r12_vint ! the intermediates to be factored
          labels(4) = form_r12_vint//'^+'
          labels(5) = form_r12_bint
          labels(6) = form_r12_xint
          nint = 4
        end if
        call set_dependency(form_ptdl0,form_r12_vint,tgt_info)
        call set_dependency(form_ptdl0,form_r12_xint,tgt_info)
        call set_dependency(form_ptdl0,form_r12_bint,tgt_info)
        if (ansatz.ne.1) then
          labels(2+nint+1) = form_r12_cint
          labels(2+nint+2) = trim(form_r12_cint)//'^+'
          call set_dependency(form_ptdl0,form_r12_cint,tgt_info)
          nint = nint+2
        end if
        call form_parameters(-1,
     &       parameters,2,title_ptdl0,nint,'---')
        call set_rule(form_ptdl0,ttype_frm,FACTOR_OUT,
     &       labels,nint+2,1,
     &       parameters,2,tgt_info)
        ! there remain a few unprocessed R12 contributions
        ! for ansatz > 1
        ! as a first resort we replace r12 by the actual integrals
        if (ansatz.gt.1) then
          call set_dependency(form_ptdl0,op_rint,tgt_info)
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = form_ptdl0
          labels(2) = form_ptdl0
          labels(3) = op_r12
          labels(4) = op_rint
          nint = 1
c          if (r12op.gt.0.and.set_R2_R2) then
c            labels(5) = 'R12VV'
c            labels(6) = op_rint
c            nint = 2
c          end if
          call form_parameters(-1,
     &         parameters,2,title_ptdl0,nint,'---')
          call set_rule(form_ptdl0,ttype_frm,REPLACE,
     &         labels,2+nint*2,1,
     &         parameters,2,tgt_info)
        end if
c dbg
        call form_parameters(-2,parameters,2,
     &       'stdout',0,'---')
        call set_rule(form_ptdl0,ttype_frm,PRINT_FORMULA,
     &       labels,1,0,
     &       parameters,2,tgt_info)
c dbg

      end if

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_ptde0
      labels(2) = form_ptdl0
      labels(3) = op_dept
      labels(4) = trim(op_tpt)//'^+'
      call add_target(form_ptde0,ttype_frm,.false.,tgt_info)
      call set_dependency(form_ptde0,form_ptdl0,tgt_info)
      call set_dependency(form_ptde0,op_dept,tgt_info)
      call set_rule(form_ptde0,ttype_frm,INVARIANT,
     &              labels,4,1,
     &              title_ptde0,1,tgt_info)

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_h0_tpt
      labels(2) = form_ptdl0
      labels(3) = op_h0_r
      labels(4) = trim(op_tpt)//'^+'
      labels(5) = ' '
      call add_target(form_h0_tpt,ttype_frm,.false.,tgt_info)
      call add_target(form_etapt,ttype_frm,.false.,tgt_info)
      call set_joined_targets(form_h0_tpt,form_etapt,tgt_info)
      call set_dependency(form_h0_tpt,form_ptdl0,tgt_info)
      call set_dependency(form_h0_tpt,op_h0_r,tgt_info)
      call set_dependency(form_h0_tpt,op_etapt,tgt_info)
      call set_rule(form_h0_tpt,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_h0_tpt,1,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_h0_tpt
      labels(2) = form_etapt
      labels(3) = form_h0_tpt
      labels(4) = op_h0_r
      labels(5) = op_etapt
      labels(6) = op_tpt
      call set_rule(form_h0_tpt,ttype_frm,LEQ_SPLIT,
     &              labels,6,2,
     &              title_h0_tpt,1,tgt_info)

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

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*
      !call get_argument_value('calculate.routes','simtraf',ival=isim)
      isim=0

      ! equations:
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = fopt_h0_tpt
      labels(2) = form_etapt
      labels(3) = form_h0_tpt
      ncat = 2
      nint = 0
      call add_target(fopt_h0_tpt,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_h0_tpt,form_ptde0,tgt_info)
      call set_dependency(fopt_h0_tpt,form_h0_tpt,tgt_info)
      call set_dependency(fopt_h0_tpt,mel_omgdef,tgt_info)
      call set_dependency(fopt_h0_tpt,mel_topdef,tgt_info)
      call set_dependency(fopt_h0_tpt,mel_ham,tgt_info)
      call set_dependency(fopt_h0_tpt,mel_dept0def,tgt_info)      
      call set_dependency(fopt_h0_tpt,mel_etaptdef,tgt_info)      
      call set_dependency(fopt_h0_tpt,mel_h0_tptdef,tgt_info)      
      call set_dependency(fopt_h0_tpt,mel_tptdef,tgt_info)      
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_h0_tpt,form_cchhat,tgt_info)
        call set_dependency(fopt_h0_tpt,mel_hhatdef,tgt_info)
        labels(4) = form_cchhat
      else if (isim.eq.2) then
        nint = 1
        call set_dependency(fopt_h0_tpt,form_cchbar,tgt_info)
        call set_dependency(fopt_h0_tpt,meldef_hbar,tgt_info)
        labels(4) = form_cchbar
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_h0_tpt,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      isim = 0
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = fopt_ptde0
      labels(2) = form_ptde0
      ncat = 1
      nint = 0
      call add_target(fopt_ptde0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_ptde0,form_ptde0,tgt_info)
      call set_dependency(fopt_ptde0,mel_topdef,tgt_info)
      call set_dependency(fopt_ptde0,mel_dept0def,tgt_info)
      call set_dependency(fopt_ptde0,mel_ham,tgt_info)
      call set_dependency(fopt_ptde0,mel_tptdef,tgt_info)      
      call set_dependency(fopt_ptde0,mel_etaptdef,tgt_info)      
      call set_dependency(fopt_ptde0,'Vpx-INTER',tgt_info)
      if (set_Rhxhh)
     &     call set_dependency(fopt_ptde0,'R12hxhh-INT',tgt_info)
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_ptde0,form_cchhat,tgt_info)
        call set_dependency(fopt_ptde0,mel_hhatdef,tgt_info)
        labels(4) = form_cchhat
      else if (isim.eq.2) then
        nint = 1
        call set_dependency(fopt_ptde0,form_cchbar,tgt_info)
        call set_dependency(fopt_ptde0,meldef_hbar,tgt_info)
        labels(4) = form_cchbar
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_ptde0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! L0/E0: 
      call add_target(mel_dlpt0,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_dlpt0,op_dlpt,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_dlpt0
      labels(2) = op_dlpt
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(mel_dlpt0,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      call add_target(mel_dept0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_dept0def,op_dept,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_dept0
      labels(2) = op_dept
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(mel_dept0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! T-list definition
      call add_target(mel_tptdef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_tptdef,op_tpt,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_tpt
      labels(2) = op_tpt
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(mel_tptdef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! H0.T-list definition 
      call add_target(mel_h0_tptdef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_h0_tptdef,op_tpt,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_h0_tpt
      labels(2) = op_h0_r
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(mel_h0_tptdef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! RHS-list definition
      call add_target(mel_etaptdef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_etaptdef,op_tpt,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_etapt
      labels(2) = op_etapt
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(mel_etaptdef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! DIA-list:
        call me_list_label(me_label,mel_dia_pt,1,0,0,0,.false.)
        call add_target(me_label,ttype_opme,.false.,tgt_info)
        call set_dependency(me_label,mel_ham,tgt_info)
        call set_dependency(me_label,op_dia_pt,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = me_label
        labels(2) = op_dia_pt
        call me_list_parameters(-1,parameters,
     &       0,0,1,0,0,.false.)
        call set_rule(me_label,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
        labels(1) = me_label
        labels(2) = mel_ham
        call set_rule(me_label,ttype_opme,PRECONDITIONER,
     &              labels,2,1,
     &              parameters,0,tgt_info)

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

      ! quick-fix:
c      call add_target('R12hxhh-INT',ttype_opme,.false.,tgt_info)
      call add_target('R12hxhh-INT',ttype_opme,set_Rhxhh,tgt_info)
      call set_dependency('R12hxhh-INT','R12hxhh',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'R12hxhh-INT'
      labels(2) = 'R12hxhh'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('R12hxhh-INT',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'R12hxhh-INT'
      call import_parameters(-1,parameters,'F12_INT',env_type)
      call set_rule('R12hxhh-INT',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)
      
      call add_target('Vpx-INTER',ttype_opme,.false.,tgt_info)
      call set_dependency('Vpx-INTER','Vpx',tgt_info)
      call set_dependency('Vpx-INTER','Vpx_CABS',tgt_info)
      call set_dependency('Vpx-INTER','G.R-X-INT',tgt_info)
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

*----------------------------------------------------------------------*
*     "phony" targets 
*----------------------------------------------------------------------*
      ! totally symmetric dia for use below:
      call me_list_label(mel_dia1,mel_dia_pt,1,0,0,0,.false.)


        ! solve GS equations
        call add_target(solve_cc_pt,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_cc_pt,mel_dia1,tgt_info)
        call set_dependency(solve_cc_pt,fopt_h0_tpt,tgt_info)
        call set_dependency(solve_cc_pt,mel_tptdef,tgt_info)
        call set_dependency(solve_cc_pt,fopt_ptde0,tgt_info)
        call solve_parameters(-1,parameters,2, 1,1,'DIA')
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_tpt
        labels(2) = mel_dia1
        labels(3) = op_h0_r
        labels(4) = op_tpt  ! no metric
        labels(5) = op_etapt
        labels(6) = fopt_h0_tpt
        call set_rule(solve_cc_pt,ttype_opme,SOLVELEQ,
     &       labels,6,1,
     &       parameters,2,tgt_info)
c      else if (.not.set_tpp) then
c        call quit(1,'set_cc_pt_targets','trap 2')
c        call add_target(solve_cc_pt,ttype_gen,.true.,tgt_info)
c        call set_dependency(solve_cc_pt,mel_dia1,tgt_info)
c        call set_dependency(solve_cc_pt,fopt_h0_tpt,tgt_info)
c        call set_dependency(solve_cc_pt,mel_tptdef,tgt_info)
c        call set_dependency(solve_cc_pt,fopt_ptde0,tgt_info)
cc        call set_dependency(solve_ccr12_gs,me_bprc,tgt_info)
cc        call set_dependency(solve_ccr12_gs,me_xprc,tgt_info)
c        call solve_parameters(-1,parameters,2, 1,1,'DIA/BLK')
c        ! baustelle: ...
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = mel_tpt
cc        labels(1) = mel_cex_pt
c        labels(2) = mel_dia1
c        labels(3) = op_h0_r
c        labels(4) = op_etapt
c        labels(5) = fopt_h0_tpt
c        call set_rule(solve_cc_pt,ttype_opme,SOLVELEQ,
c     &       labels,5,1,
c     &       parameters,2,tgt_info)        
c      else
c        call quit(1,'set_cc_pt_targets','trap 1')
c      end if
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = fopt_ptde0
      call set_rule(solve_cc_pt,ttype_opme,EVAL,
     &       labels,1,0,
     &       parameters,0,tgt_info)
        
      return
      end
