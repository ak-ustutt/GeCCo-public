*----------------------------------------------------------------------*
      subroutine set_ccr12f_special_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed specifically in CC-R12 calculations with
*     fixed geminals
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'
      include 'opdim.h'

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
     &     min_rank, max_rank, ansatz, nop,
     &     isim, ncat, nint, icnt, ndef, extend,
     &     isym, ms, msc, sym_arr(8), nlabel,
     &     ninproj, navoid, nconnect,
     &     connect(20), idx_sv(20), iblkmin(20),
     &     iblkmax(20),
     &     occ_def(ngastp,2,20)
      logical ::
     &     needed, r12fix
      character(8) ::
     &     approx
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)

      character, parameter ::
     &     form_r_t*8 = 'FORM_R_T',
     &     op_r_t*6   = 'OP_R_T',
     &     me_r_t*6   = 'ME_R_T',
     &     medef_r_t*10  = 'DEF-ME_R_T'
      
      if (iprlvl.gt.0)
     &     write(luout,*) 'setting special targets for CC-R12 ...'

      approx = '        '
      ! read keyword values
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','extend',ival=extend)

      call get_argument_value('calculate.routes','simtraf',ival=isim)

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

*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*
      call add_target(form_ccr12lg0,ttype_frm,.false.,tgt_info)
      ! (a) set formal Lagrangian (in 'complete' basis)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12lg0
      labels(2) = op_ccr12lg
      labels(3) = op_ham
      labels(4) = op_r12
      labels(5) = op_r12
      labels(6) = op_tbar
      labels(7) = op_top
      nlabel = 7
      call set_dependency(form_ccr12lg0,op_ccr12lg,tgt_info)
      call set_dependency(form_ccr12lg0,op_ham,tgt_info)
      call set_dependency(form_ccr12lg0,op_r12,tgt_info)
c      call set_dependency(form_ccr12lg0,op_rba,tgt_info)
      call set_dependency(form_ccr12lg0,op_tbar,tgt_info)
      call set_dependency(form_ccr12lg0,op_top,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_ccr12lg0,ansatz,'---')
      call set_rule(form_ccr12lg0,ttype_frm,DEF_CCR12_LAGRANGIAN,
     &              labels,nlabel,1,
     &              parameters,2,tgt_info)
      ! (b) Factor out the R12 intermediates 
      ! (effectively removing all reference to the complete basis)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12lg0 ! output formula (itself)
      labels(2) = form_ccr12lg0 ! input formula
      labels(3) = form_r12_vint    ! the intermediates to be factored
      labels(4) = form_r12_vint//'^+'
      labels(5) = form_r12_bint
      labels(6) = form_r12_bhint
      labels(7) = form_r12_xint
      nint = 5
      call set_dependency(form_ccr12lg0,form_r12_vint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_xint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_bint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_bhint,tgt_info)
      if (ansatz.ne.1) then
        labels(8) = form_r12_cint
        labels(9) = trim(form_r12_cint)//'^+'
        call set_dependency(form_ccr12lg0,form_r12_cint,tgt_info)
        nint = 7
      end if
      call form_parameters(-1,
     &     parameters,2,title_ccr12lg0,nint,'---')
      call set_rule(form_ccr12lg0,ttype_frm,FACTOR_OUT,
     &              labels,nint+2,1,
     &              parameters,2,tgt_info)
      ! (c) post-processing: remove terms which do not contribute for
      !     the given R12-approximation
      ! .... to come
      ! fix for A
      if (trim(approx).eq.'A') then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_ccr12lg0
        labels(2) = form_ccr12lg0
        labels(3) = op_ccr12lg
        labels(4) = op_x_inter
        call set_rule(form_ccr12lg0,ttype_frm,INVARIANT,
     &              labels,4,1,
     &              title_ccr12lg0,1,tgt_info)
      end if

      ! there remain a few unprocessed R12 contributions
      ! for ansatz > 1
      ! as a first resort we replace r12 by the actual integrals
      if (ansatz.gt.1) then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_ccr12lg0
        labels(2) = form_ccr12lg0
        labels(3) = op_r12
        labels(4) = op_rint
        call set_dependency(form_ccr12lg0,op_rint,tgt_info)
        call form_parameters(-1,
     &       parameters,2,title_ccr12lg0,1,'---')
        call set_rule(form_ccr12lg0,ttype_frm,REPLACE,
     &              labels,4,1,
     &              parameters,2,tgt_info)
      end if
      
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12en0
      labels(2) = form_ccr12lg0
      labels(3) = op_ccr12en
      labels(4) = op_tbar
      nlabel = 4
      call add_target(form_ccr12en0,ttype_frm,.true.,tgt_info)
      call set_dependency(form_ccr12en0,form_ccr12lg0,tgt_info)
      call set_dependency(form_ccr12en0,op_ccr12en,tgt_info)
      call set_rule(form_ccr12en0,ttype_frm,INVARIANT,
     &              labels,nlabel,1,
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

      call add_target(op_r_t,ttype_frm,.false.,tgt_info)
      ndef = 1
      call r12gem_parameters(-1,parameters,
     &     0,2,max_rank,ansatz)
      call set_rule(op_r_t,ttype_op,DEF_R12GEMINAL,
     &              op_r_t,1,1,
     &              parameters,1,tgt_info)
c      occ_def(IPART,1,1) = 1
c      occ_def(IEXTR,1,1) = 1
c      occ_def(IHOLE,2,1) = 2
c      call op_from_occ_parameters(-1,parameters,2,
c     &     occ_def,ndef,1,ndef)
c      call set_rule(op_r_t,ttype_op,DEF_OP_FROM_OCC,
c     &              op_r_t,1,1,
c     &              parameters,2,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_r_t
      labels(2) = op_r_t
      labels(3) = op_r_t
      labels(4) = op_rint
      labels(5) = op_top
      labels(6) = op_r_t
      idx_sv(1:4) = (/1,2,3,1/)
      iblkmin(1:4) = (/1,10,1,1/)
      iblkmax(1:4) = (/0,10,0,0/)
      nconnect = 1
      connect(1:2) = (/2,3/)
      navoid = 0
      ninproj = 0
      call add_target(form_r_t,ttype_frm,.true.,tgt_info)
      call set_dependency(form_r_t,op_rint,tgt_info)
      call set_dependency(form_r_t,op_top,tgt_info)
      call set_dependency(form_r_t,op_r_t,tgt_info)
      call expand_parameters(-1,
     &       parameters,3,
     &       'XXX',4,idx_sv,iblkmin,iblkmax,
     &       connect,nconnect,
     &       0,navoid,
     &       0,ninproj)
      call set_rule(form_r_t,ttype_frm,EXPAND_OP_PRODUCT,
     &              labels,6,1,
     &              parameters,3,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*

      ! CC ground state:
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = fopt_ccr12_0
      labels(2) = form_ccr12en0
      labels(3) = form_ccr12rs_t
      ncat = 2
      nint = 0
      call add_target(fopt_ccr12_0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_ccr12_0,form_ccr12en0,tgt_info)
      call set_dependency(fopt_ccr12_0,form_ccr12rs_t,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_omgdef,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_topdef,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_ham,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_rint,tgt_info)      
      call set_dependency(fopt_ccr12_0,mel_v_def,tgt_info)      
      call set_dependency(fopt_ccr12_0,mel_b_def,tgt_info)      
      call set_dependency(fopt_ccr12_0,mel_bh_def,tgt_info)      
      call set_dependency(fopt_ccr12_0,mel_x_def,tgt_info)      
      call set_dependency(fopt_ccr12_0,mel_c_def,tgt_info)      
      call set_dependency(fopt_ccr12_0,mel_ccr12en0def,tgt_info)      
      if (extend.gt.0) then
        nint = 1
        call set_dependency(fopt_ccr12_0,form_r_t,tgt_info)
        call set_dependency(fopt_ccr12_0,medef_r_t,tgt_info)
        labels(4) = form_r_t
      end if
      if (isim.eq.1) then
        nint = nint+ 1
        call set_dependency(fopt_ccr12_0,form_cchhat,tgt_info)
        call set_dependency(fopt_ccr12_0,mel_hhatdef,tgt_info)
        labels(ncat+nint+1) = form_cchhat
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_ccr12_0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*
      ! L0/E0:
      call add_target(mel_ccr12lg0,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ccr12lg0,op_ccr12lg,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_ccr12lg0
      labels(2) = op_ccr12lg
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_ccr12lg0,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call add_target(mel_ccr12en0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ccr12en0def,op_ccr12en,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_ccr12en0
      labels(2) = op_ccr12en
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_ccr12en0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      call add_target(medef_r_t,ttype_opme,.false.,tgt_info)
      call set_dependency(medef_r_t,op_r_t,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = me_r_t
      labels(2) = op_r_t
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(medef_r_t,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      ! totally symmetric dia for use below:
      call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

      call add_target(solve_ccr12_gs,ttype_gen,.true.,tgt_info)
      call set_dependency(solve_ccr12_gs,mel_dia1,tgt_info)
      call set_dependency(solve_ccr12_gs,fopt_ccr12_0,tgt_info)
      call set_dependency(solve_ccr12_gs,eval_r12_inter,tgt_info)
      call solve_parameters(-1,parameters,2, 1,1,'DIA')
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_top
      labels(2) = mel_omg
      labels(3) = mel_dia1      ! dummy
c        labels(6) = mel_b_dia
      labels(4) = mel_ccr12en0
      labels(5) = fopt_ccr12_0
      labels(6) = mel_ham
      call set_rule(solve_ccr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,6,2,
     &       parameters,2,tgt_info)

      return
      end
