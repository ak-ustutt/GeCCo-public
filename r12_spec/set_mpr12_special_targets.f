*----------------------------------------------------------------------*
      subroutine set_mpr12_special_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed specifically in MP-R12 calculations
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
     &     min_rank, max_rank, level,
     &     isim, ncat, nint, icnt, ansatz,
     &     isym, ms, msc, sym_arr(8)
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(2)
      character(12) ::
     &     approx

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting special targets for MP-R12 ...'

      ! get keyword values
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.MP','level',ival=level)
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      approx(1:12) = ' '
      call get_argument_value('method.R12','approx',str=approx)

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! Lagrange functional
      call add_target(op_mpr12lg,ttype_op,.false.,tgt_info)
      call set_rule(op_mpr12lg,ttype_op,DEF_SCALAR,
     &              op_mpr12lg,1,1,
     &              parameters,0,tgt_info)
      
      ! Energy
      call add_target(op_mpr12en,ttype_op,.false.,tgt_info)
      call set_rule(op_mpr12en,ttype_op,DEF_SCALAR,
     &              op_mpr12en,1,1,
     &              parameters,0,tgt_info)

      ! residual
      call add_target(op_omgr12,ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &     .false.,min_rank,max_rank,0,2)
      call set_rule(op_omgr12,ttype_op,DEF_R12INTERM,
     &              op_omgr12,1,1,
     &              parameters,1,tgt_info)

c      ! diagonal
c      call add_target(op_diar12,ttype_op,.false.,tgt_info)
c      call set_dependency(op_diar12,op_omgr12,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_omgr12,.false.) ! <- dagger=.false.
c      call set_rule(op_diar12,ttype_op,CLONE_OP,
c     &              op_diar12,1,1,
c     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*

      call add_target(form_mpr12lg0,ttype_frm,.false.,tgt_info)
      ! (a) set formal Lagrangian (in 'complete' basis)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12lg0
      labels(2) = op_mpr12lg
      labels(3) = op_ham
      labels(4) = op_r12
      labels(5) = op_r12
      labels(6) = op_tbar
      labels(7) = op_cba
      labels(8) = op_top
      labels(9) = op_c12
      call set_dependency(form_mpr12lg0,op_mpr12lg,tgt_info)
      call set_dependency(form_mpr12lg0,op_ham,tgt_info)
      call set_dependency(form_mpr12lg0,op_r12,tgt_info)
c      call set_dependency(form_mpr12lg0,op_rba,tgt_info)
      call set_dependency(form_mpr12lg0,op_tbar,tgt_info)
      call set_dependency(form_mpr12lg0,op_top,tgt_info)
      call set_dependency(form_mpr12lg0,op_cba,tgt_info)
      call set_dependency(form_mpr12lg0,op_c12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_mpr12lg0,level,'---')
      call set_rule(form_mpr12lg0,ttype_frm,DEF_MPR12_LAGRANGIAN,
     &              labels,9,1,
     &              parameters,2,tgt_info)
      ! (b) Factor out the R12 intermediates 
      ! (effectively removing all reference to the complete basis)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12lg0 ! output formula (itself)
      labels(2) = form_mpr12lg0 ! input formula
      labels(3) = form_r12_vint    ! the intermediates to be factored
      labels(4) = trim(form_r12_vint)//'^+'
c      labels(4) = form_r12_vbint
      labels(5) = form_r12_xint
      labels(6) = form_r12_bint
      nint = 4
      call set_dependency(form_mpr12lg0,form_r12_vint,tgt_info)
c      call set_dependency(form_mpr12lg0,form_r12_vbint,tgt_info)
      call set_dependency(form_mpr12lg0,form_r12_xint,tgt_info)
      call set_dependency(form_mpr12lg0,form_r12_bint,tgt_info)
      if (ansatz.ne.1) then
        labels(7) = form_r12_cint
        labels(8) = trim(form_r12_cint)//'^+'
c        labels(8) = form_r12_cbint
        call set_dependency(form_mpr12lg0,form_r12_cint,tgt_info)
c        call set_dependency(form_mpr12lg0,form_r12_cbint,tgt_info)
        nint = 6
      end if
      call form_parameters(-1,
     &     parameters,2,title_mpr12lg0,nint,'---')
      call set_rule(form_mpr12lg0,ttype_frm,FACTOR_OUT,
     &              labels,nint+2,1,
     &              parameters,2,tgt_info)
      ! (c) post-processing: remove terms which do not contribute for
      !     the given R12-approximation
      ! .... to come
      ! fix for A
      if (trim(approx).eq.'A') then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_mpr12lg0
        labels(2) = form_mpr12lg0
        labels(3) = op_mpr12lg
        labels(4) = op_x_inter
        call set_rule(form_mpr12lg0,ttype_frm,INVARIANT,
     &              labels,4,1,
     &              title_mpr12lg0,1,tgt_info)
      end if

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12en0
      labels(2) = form_mpr12lg0
      labels(3) = op_mpr12en
      labels(4) = op_tbar
      labels(5) = op_cba
      call add_target(form_mpr12en0,ttype_frm,.true.,tgt_info)
      call set_dependency(form_mpr12en0,form_mpr12lg0,tgt_info)
      call set_dependency(form_mpr12en0,op_mpr12en,tgt_info)
      call set_rule(form_mpr12en0,ttype_frm,INVARIANT,
     &              labels,5,1,
     &              title_mpr12en0,1,tgt_info)

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12rs_t
      labels(2) = form_mpr12lg0
      labels(3) = op_omg
      labels(4) = op_tbar
      labels(5) = ' '
      call add_target(form_mpr12rs_t,ttype_frm,.true.,tgt_info)
      call set_dependency(form_mpr12rs_t,form_mpr12lg0,tgt_info)
      call set_dependency(form_mpr12rs_t,op_omg,tgt_info)
      call set_rule(form_mpr12rs_t,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_mpr12rs_t,1,tgt_info)

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12rs_c
      labels(2) = form_mpr12lg0
      labels(3) = op_omgr12
      labels(4) = op_cba
      labels(5) = ' '
      call add_target(form_mpr12rs_c,ttype_frm,.true.,tgt_info)
      call set_dependency(form_mpr12rs_c,form_mpr12lg0,tgt_info)
      call set_dependency(form_mpr12rs_c,op_omgr12,tgt_info)
      call set_rule(form_mpr12rs_c,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_mpr12rs_c,1,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*
      ! MP ground state:
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = fopt_mpr12_0
      labels(2) = form_mpr12en0
      labels(3) = form_mpr12rs_t
      labels(4) = form_mpr12rs_c
      ncat = 3
      nint = 0
      call add_target(fopt_mpr12_0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_mpr12_0,form_mpr12en0,tgt_info)
      call set_dependency(fopt_mpr12_0,form_mpr12rs_t,tgt_info)
      call set_dependency(fopt_mpr12_0,form_mpr12rs_c,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_omgdef,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_topdef,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_omgr12def,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_c12def,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_ham,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_mpr12en0def,tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_mpr12_0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)


*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*
      ! L0/E0:
      call add_target(mel_mpr12lg0,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_mpr12lg0,op_mpr12lg,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_mpr12lg0
      labels(2) = op_mpr12lg
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_mpr12lg0,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call add_target(mel_mpr12en0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_mpr12en0def,op_mpr12en,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_mpr12en0
      labels(2) = op_mpr12en
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_mpr12en0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! CBAR list definition
      call add_target(mel_cbardef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_cbardef,op_cba,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_cbar
      labels(2) = op_cba
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_cbardef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! C12  list definition
      call add_target(mel_c12def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_c12def,op_c12,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_c12
      labels(2) = op_c12
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_c12def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! OMG-R12 list definition
      call add_target(mel_omgr12def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_omgr12def,op_omgr12,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_omgr12
      labels(2) = op_omgr12
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_omgr12def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      ! totally symmetric dia for use below:
      call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

      call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
      call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
      call set_dependency(solve_mpr12_gs,mel_b_inv,tgt_info)
c      call set_dependency(solve_mpr12_gs,mel_b_dia,tgt_info)
c      call set_dependency(solve_mpr12_gs,mel_x_inv,tgt_info)
      call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
      call solve_parameters(-1,parameters,2, 2,1,'DIA/BLK')
c      call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_top
      labels(2) = mel_c12
      labels(3) = mel_omg
      labels(4) = mel_omgr12
      labels(5) = mel_dia1
      labels(6) = mel_dia1 ! dummy
c      labels(6) = mel_b_dia
      labels(7) = mel_mpr12en0
      labels(8) = fopt_mpr12_0
      if(trim(approx).eq.'A')then
        labels(9) = mel_b_inv   ! or mel_b_inter
      else
        labels(9) = mel_b_inter
      endif
      labels(10) = mel_x_inter
      labels(11) = mel_ham
      call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
     &     labels,11,4,
     &     parameters,2,tgt_info)

      return
      end
