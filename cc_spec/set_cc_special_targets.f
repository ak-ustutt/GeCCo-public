*----------------------------------------------------------------------*
      subroutine set_cc_special_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed specifically in CC calculations
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
     &     min_rank, max_rank,
     &     isim, ncat, nint, icnt,
     &     isym, ms, msc, sym_arr(8)
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters(2)

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting special targets for CC ...'

      msc = +1 ! assuming closed shell

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
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

*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cclg0
      labels(2) = op_cclg
      labels(3) = op_ham
      labels(4) = op_tbar
      labels(5) = op_top
      call add_target(form_cclg0,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cclg0,op_cclg,tgt_info)
      call set_dependency(form_cclg0,op_ham,tgt_info)
      call set_dependency(form_cclg0,op_tbar,tgt_info)
      call set_dependency(form_cclg0,op_top,tgt_info)
      call set_rule(form_cclg0,ttype_frm,DEF_CC_LAGRANGIAN,
     &              labels,5,1,
     &              title_cclg0,1,tgt_info)
      call set_rule(form_cclg0,ttype_frm,TEX_FORMULA,
     &              labels,5,1,
     &              'cc_lag.tex',1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccen0
      labels(2) = form_cclg0
      labels(3) = op_ccen
      labels(4) = op_tbar
      call add_target(form_ccen0,ttype_frm,.false.,tgt_info)
      call set_dependency(form_ccen0,form_cclg0,tgt_info)
      call set_dependency(form_ccen0,op_ccen,tgt_info)
      call set_rule(form_ccen0,ttype_frm,INVARIANT,
     &              labels,4,1,
     &              title_ccen0,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccrs0
      labels(2) = form_cclg0
      labels(3) = op_omg
      labels(4) = op_tbar
      labels(5) = ' '
      call add_target(form_ccrs0,ttype_frm,.false.,tgt_info)
      call set_dependency(form_ccrs0,form_cclg0,tgt_info)
      call set_dependency(form_ccrs0,op_omg,tgt_info)
      call set_rule(form_ccrs0,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_ccrs0,1,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*
      call get_argument_value('calculate.routes','simtraf',ival=isim)

      ! CC ground state:
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_ccrs0
      labels(2) = form_ccen0
      labels(3) = form_ccrs0
      ncat = 2
      nint = 0
      call add_target(fopt_ccrs0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_ccrs0,form_ccen0,tgt_info)
      call set_dependency(fopt_ccrs0,form_ccrs0,tgt_info)
      call set_dependency(fopt_ccrs0,mel_omgdef,tgt_info)
      call set_dependency(fopt_ccrs0,mel_topdef,tgt_info)
      call set_dependency(fopt_ccrs0,mel_ham,tgt_info)
      call set_dependency(fopt_ccrs0,mel_ccen0def,tgt_info)      
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_ccrs0,form_cchhat,tgt_info)
        call set_dependency(fopt_ccrs0,mel_hhatdef,tgt_info)
        labels(4) = form_cchhat
      else if (isim.eq.2) then
        nint = 1
        call set_dependency(fopt_ccrs0,form_cchbar,tgt_info)
        call set_dependency(fopt_ccrs0,meldef_hbar,tgt_info)
        labels(4) = form_cchbar
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_ccrs0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! L0/E0: 
      call add_target(mel_cclg0,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_cclg0,op_cclg,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_cclg0
      labels(2) = op_cclg
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_cclg0,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call add_target(mel_ccen0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ccen0def,op_ccen,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ccen0
      labels(2) = op_ccen
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_ccen0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

c      ! Hbar definition
c      call add_target(meldef_hbar,ttype_opme,.false.,tgt_info)
c      call set_dependency(meldef_hbar,op_hbar,tgt_info)
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_hbar
c      labels(2) = op_hbar
c      call me_list_parameters(-1,parameters,
c     &     0,0,1,0,0)
c      call set_rule(meldef_hbar,ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)

      
*----------------------------------------------------------------------*
*     "phony" targets 
*----------------------------------------------------------------------*
      ! totally symmetric dia for use below:
      call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

      ! solve GS equations
      call add_target(solve_cc_gs,ttype_gen,.true.,tgt_info)
      call set_dependency(solve_cc_gs,mel_dia1,tgt_info)
      call set_dependency(solve_cc_gs,fopt_ccrs0,tgt_info)
      call solve_parameters(-1,parameters,2, 1,1,'DIA')
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_top
      labels(2) = mel_omg
      labels(3) = mel_dia1
      labels(4) = mel_ccen0
      labels(5) = fopt_ccrs0
      call set_rule(solve_cc_gs,ttype_opme,SOLVENLEQ,
     &     labels,5,2,
     &     parameters,2,tgt_info)

      return
      end
