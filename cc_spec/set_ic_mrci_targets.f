*----------------------------------------------------------------------*
      subroutine set_ic_mrci_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
*     preliminary routine for multireference targets.
*     installed in order to get some complicated factorizations
*     (both FACTOR_OUT and OPTIMIZE) into the test suite
*
*     matthias, end of 2009
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
     &     min_rank, max_rank, ndef, occ_def(ngastp,2,60),
     &     isim, ncat, nint, icnt,
     &     isym, ms, msc, sym_arr(8), maxexc, ip, ih, ivv, iv,
     &     cminh, cmaxh, cminp, cmaxp, cmaxexc, minh, maxh,
     &     minp, maxp, maxv, maxvv, minexc, cbarc(2),
     &     check_S_roots
      logical ::
     &     use_jacobi, use_cbarc, use_tr, use_dens
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)

      ! first set targets for CASSCF or uncontracted CI wave function
        call set_unc_mrci_targets(tgt_info,orb_info)

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting multireference targets #2...'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1 ! assuming closed shell

      ! get minimum and maximum numbers of excitations, holes, particles,
      ! valence-valence excitations
      call get_argument_value('calculate.multiref','minh',
     &     ival=minh)
      call get_argument_value('calculate.multiref','maxh',
     &     ival=maxh)
      call get_argument_value('calculate.multiref','minp',
     &     ival=minp)
      call get_argument_value('calculate.multiref','maxp',
     &     ival=maxp)
      call get_argument_value('calculate.multiref','maxv',
     &     ival=maxv)
      call get_argument_value('calculate.multiref','maxvv',
     &     ival=maxvv)
      call get_argument_value('calculate.multiref','minexc',
     &     ival=minexc)
      call get_argument_value('calculate.multiref','maxexc',
     &     ival=maxexc)
      if (maxh.lt.0) maxh = maxexc
      if (maxp.lt.0) maxp = maxexc
      if (maxv.lt.0) maxv = 2*maxexc
      if (maxvv.lt.0) maxvv = maxexc

      ! if maxexc = 0: return because call of unc_mrci is sufficient
      if (maxexc.eq.0) return

c dbg
      print *,'minh    = ',minh
      print *,'maxh    = ',maxh
      print *,'minp    = ',minp
      print *,'maxp    = ',maxp
      print *,'maxv    = ',maxv
      print *,'maxvv   = ',maxvv
      print *,'minexc  = ',minexc
      print *,'maxexc  = ',maxexc
c dbgend

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define particle conserving excitation operator C (for icMRCI)
      call add_target('C',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (ivv.gt.0.and.ip.eq.0.and.ih.eq.0) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).gt.0.and.ip.lt.minp) cycle
            if (max(ip,ih).gt.0.and.ih.lt.minh) cycle
            if (max(ip,ih).gt.0.and.max(ip,ih)+ivv.lt.minexc) cycle
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('C',ttype_op,DEF_OP_FROM_OCC,
     &              'C',1,1,
     &              parameters,2,tgt_info)

      ! define unit operator suited to connect C^+ and C
      call add_target('1',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 0, maxexc
        ndef = ndef + 1
        occ_def(IVALE,1,ndef) = iv
        occ_def(IVALE,2,ndef) = iv
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('1',ttype_op,DEF_OP_FROM_OCC,
     &              '1',1,1,
     &              parameters,2,tgt_info)
      call opt_parameters(-1,parameters,+1,0)
      call set_rule('1',ttype_op,SET_HERMIT,
     &              '1',1,1,
     &              parameters,1,tgt_info)

      ! define scalar norm
      call add_target('NORM',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('NORM',ttype_op,DEF_HAMILTONIAN,'NORM',
     &              1,1,parameters,1,tgt_info)

      ! define density matrix
      call add_target('D',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (ip.eq.2.and.ih.eq.2) cycle
            if (ip.eq.0.and.ih.eq.0.and.ivv.gt.0) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).gt.0.and.ip.lt.minp) cycle
            if (max(ip,ih).gt.0.and.ih.lt.minh) cycle
            if (max(ip,ih).gt.0.and.max(ip,ih)+ivv.lt.2) cycle ! only doubles
            ndef = ndef + 1
            occ_def(IVALE,1,ndef*3-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef*3-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,1,ndef*3) = max(ip-ih,0) + ivv
            occ_def(IVALE,2,ndef*3-2) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      use_dens = ndef.gt.1
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,3,(/0,0,0,0,0,0/),ndef)
      call set_rule('D',ttype_op,DEF_OP_FROM_OCC,
     &              'D',1,1,
     &              parameters,2,tgt_info)

      ! define product Jacobian times C
      call add_target('A_C',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (ivv.gt.0.and.ip.eq.0.and.ih.eq.0) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).gt.0.and.ip.lt.minp) cycle
            if (max(ip,ih).gt.0.and.ih.lt.minh) cycle
            if (max(ip,ih).gt.0.and.max(ip,ih)+ivv.lt.minexc) cycle
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef*2) = ih
            occ_def(IPART,1,ndef*2) = ip
            occ_def(IVALE,1,ndef*2) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef*2-1) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('A_C',ttype_op,DEF_OP_FROM_OCC,
     &              'A_C',1,1,
     &              parameters,2,tgt_info)

      ! Metric times icCI vector product
      call add_target('SC',ttype_op,.false.,
     &                tgt_info)
      call set_dependency('SC','A_C',tgt_info)
      call cloneop_parameters(-1,parameters,'A_C',.false.)
      call set_rule('SC',ttype_op,CLONE_OP,'SC',1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*

      ! expression for the norm
c dbg
      print *,'define norm'
c dbgend
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_NORM'
      labels(2) = 'NORM'
      labels(3) = 'C0^+'
      labels(4) = 'C^+'
      labels(5) = '1'
      labels(6) = 'C'
      labels(7) = 'C0'
      call add_target('F_NORM',ttype_frm,.false.,tgt_info)
      call set_dependency('F_NORM','NORM',tgt_info)
      call set_dependency('F_NORM','C0',tgt_info)
      call set_dependency('F_NORM','C',tgt_info)
      call set_dependency('F_NORM','1',tgt_info)
      call expand_parameters(-1,
     &     parameters,3,
     &     'multireference energy expression',5,
     &     (/2,3,4,5,6/),
     &     (/-1,-1,-1,-1,-1/),
     &     (/0,0,3,0,0/),
     &     0,0,
     &     (/1,3,3,5/),2,
     &     0,0)
c dbg
      print *,'give label',trim(labels(1))
c dbgend
      call set_rule('F_NORM',ttype_frm,EXPAND_OP_PRODUCT,
     &              labels,7,1,
     &              parameters,3,tgt_info)
      ! delete terms in which C^+ and C are contracted via valence lines
      labels(2:10)(1:len_target_name) = ' '
      labels(2) = 'F_NORM'
      labels(3) = 'NORM'
      labels(4) = 'C^+'
      labels(5) = 'C'
      call form_parameters(-1,
     &       parameters,2,'norm expression',3,'delete')
      call set_rule('F_NORM',ttype_frm,SELECT_LINE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      labels(3:10)(1:len_target_name) = ' '
      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
      call set_rule('F_NORM',ttype_frm,PRINT_FORMULA,
     &                labels,2,1,parameters,2,tgt_info)

      ! factor out reduced density matrices
      if (use_dens) then
        call add_target('F_NORM_fact',ttype_frm,.false.,tgt_info)
        call set_dependency('F_NORM_fact','F_NORM',tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = 'F_NORM_fact'
        labels(2) = 'F_NORM'
        labels(3) = 'F_D'
        call set_dependency('F_NORM_fact','F_D',tgt_info)
        call form_parameters(-1,
     &       parameters,2,'norm using density matrices',1,'---')
        call set_rule('F_NORM_fact',ttype_frm,FACTOR_OUT,
     &                labels,3,1,
     &                parameters,2,tgt_info)
        call form_parameters(-1,parameters,2,'stdout',1,'stdout')
        call set_rule('F_NORM_fact',ttype_frm,PRINT_FORMULA,
     &                  labels,2,1,parameters,2,tgt_info)
      end if

      ! metric times icCI coefficient vector
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'F_SC'
      labels(2) = 'F_NORM'
      labels(3) = 'SC'
      labels(4) = 'C^+'
      labels(5) = ' '
      call add_target('F_SC',ttype_frm,.false.,tgt_info)
      if (use_dens) then
        call set_dependency('F_SC','F_NORM_fact',tgt_info)
        labels(2) = 'F_NORM_fact'
      else
        call set_dependency('F_SC','F_NORM',tgt_info)
      end if
      call set_dependency('F_SC','SC',tgt_info)
      call set_dependency('F_SC','C',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'metric times right hand vector',1,'---')
      call set_rule('F_SC',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
      call set_rule('F_SC',ttype_frm,PRINT_FORMULA,
     &                labels,2,1,parameters,2,tgt_info)

      ! reduced density matrices
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_SC0'
      labels(2) = 'F_NORM'
      labels(3) = 'SC'
      labels(4) = 'C^+'
      labels(5) = ' '
      call add_target('F_SC0',ttype_frm,.false.,tgt_info)
      call set_dependency('F_SC0','F_NORM',tgt_info)
      call set_dependency('F_SC0','SC',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'precursor for density matrices',1,'---')
      call set_rule('F_SC0',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
      call set_rule('F_SC0',ttype_frm,PRINT_FORMULA,
     &                labels,2,1,parameters,2,tgt_info)

      ! reduced density matrices
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_D'
      labels(2) = 'F_SC0'
      labels(3) = 'D'
      labels(4) = 'C'
      labels(5) = ' '
      call add_target('F_D',ttype_frm,.false.,tgt_info)
      call set_dependency('F_D','F_SC0',tgt_info)
      call set_dependency('F_D','D',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'Reduced density matrices',1,'---')
      call set_rule('F_D',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      ! delete terms which do not have any open lines (valence space)
      ! => scalar occupation class is excluded in the definition
      labels(2:20)(1:len_target_name) = ' '
      labels(2) = 'F_D'
      labels(3) = 'D'
      labels(4) = '1'
      labels(5) = 'C0'
      call form_parameters(-1,
     &       parameters,2,'D formula',3,'ext')
      call set_rule('F_D',ttype_frm,SELECT_LINE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
      call set_rule('F_D',ttype_frm,PRINT_FORMULA,
     &                labels,2,1,parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! density matrix
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_D'
      labels(2) = 'F_D'
      call add_target('FOPT_D',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_D','F_D',tgt_info)
      call set_dependency('FOPT_D','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_D','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_D','DEF_ME_D',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_D',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! SC
      call add_target('FOPT_SC',ttype_frm,.true.,tgt_info)
      call set_dependency('FOPT_SC','F_SC',tgt_info)
      call set_dependency('FOPT_SC','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_SC','DEF_ME_C',tgt_info)
      call set_dependency('FOPT_SC','DEF_ME_SC',tgt_info)
      call set_dependency('FOPT_SC','DEF_ME_1',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_SC'
      labels(2) = 'SC'
      call set_rule('FOPT_SC',ttype_opme,ASSIGN_ME2OP,
     &     labels,2,1,
     &     parameters,0,tgt_info)
      if (use_dens)
     &        call set_dependency('FOPT_SC','DEF_ME_D',tgt_info)
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_SC'
      labels(2) = 'F_SC'
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_SC',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! ME_C
      call add_target('DEF_ME_C',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_C','C',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_C'
      labels(2) = 'C'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_C',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! ME_1
      call add_target('DEF_ME_1',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_1','1',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_1'
      labels(2) = '1'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_1',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call dens_parameters(-1,parameters,0,0,0)
      call set_rule('DEF_ME_1',ttype_opme,UNITY,
     &     labels,1,1,
     &     parameters,1,tgt_info)

      ! ME_SC
      call add_target('DEF_ME_SC',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_SC','SC',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_SC'
      labels(2) = 'SC'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_SC',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! ME_D
      call add_target('DEF_ME_D',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_D','D',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_D'
      labels(2) = 'D'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_D',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! Evaluate density matrix
      call add_target('EVAL_D',ttype_gen,.true.,tgt_info)
      call set_dependency('EVAL_D','FOPT_D',tgt_info)
      call set_dependency('EVAL_D','SOLVE_REF',tgt_info)
      call set_rule('EVAL_D',ttype_opme,EVAL,
     &     'FOPT_D',1,0,
     &     parameters,0,tgt_info)

      return
      end
