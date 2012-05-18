*----------------------------------------------------------------------*
      subroutine set_gno_targets(tgt_info,orb_info,maxtop)
*----------------------------------------------------------------------*
*     set targets for use of generalized normal ordering
*
*     matthias summer 2010
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_target_info.h'
      include 'ifc_targets.h'
      include 'def_orbinf.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'

      integer, parameter ::
     &     ntest = 100

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     maxtop

      integer ::
     &     ndef, occ_def(ngastp,2,60), msc, maxexc, iv, maxcum,
     &     ioff, gno, nv, cum_appr_mode, idef
      logical ::
     &     pure_vv
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character ::
     &     op_dint*7, f_dint*9, defme_dint*14, me_dint*10
      real(8) ::
     &     factor

      ! first set targets for CASSCF or uncontracted CI wave function
      ! (if not done already)
      if (.not.is_keyword_set('method.MR').gt.0)
     &      call quit(1,'set_gno_targets',
     &      'generalized normal order requires MR wave function')

      if (iprlvl.gt.0)
     &     write(luout,*) 
     &     'setting targets for generalized normal ordering...'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1
      if (orb_info%ims.ne.0) msc = 0

      ! get minimum and maximum numbers of excitations, holes, particles,
      ! valence-valence excitations
      call get_argument_value('method.MR','maxexc',
     &     ival=maxexc)
      call get_argument_value('method.MR','maxcum',
     &     ival=maxcum)
      call get_argument_value('method.MR','cum_appr_mode',
     &     ival=cum_appr_mode)
      call get_argument_value('method.MR','GNO',
     &     ival=gno)
      call get_argument_value('method.MR','pure_vv',
     &     lval=pure_vv)
      if (cum_appr_mode.lt.0.or.cum_appr_mode.gt.3)
     &   call quit(1,'set_gno_targets','cum_appr_mode must be 0,1,2,3')
      ioff = 0
      if (gno.eq.0) then
        ioff = 2 !3
        if (pure_vv) ioff = 3 !2
      end if

      if (ntest.ge.100) then
        write(luout,*) 'maxcum       = ',maxcum
        if (maxtop.gt.0)
     &     write(luout,*) 'cum_appr_mode= ',cum_appr_mode
      end if

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define (simple) reduced densities
c      call add_target('DENS',ttype_op,.false.,tgt_info)
      call add_target2('DENS',.false.,tgt_info)
      occ_def = 0
      ndef = 0
c      do iv = 1, 2 + maxexc*(maxtop+1) - ioff
      if (maxtop.le.0) then
        nv = min(maxcum,orb_info%nactel)
      else
        nv = min(ioff + (maxexc-1)*(maxtop+1),orb_info%nactel)
      end if
      do iv = 1, nv
            ndef = ndef + 1
            occ_def(IVALE,1,ndef*2) = iv
            occ_def(IVALE,2,ndef*2-1) = iv
      end do
c      call op_from_occ_parameters(-1,parameters,2,
c     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
c      call set_rule('DENS',ttype_op,DEF_OP_FROM_OCC,
c     &              'DENS',1,1,
c     &              parameters,2,tgt_info)
      call set_rule2('DENS',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('DENS',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'DENS'/))
      call set_arg('DENS',DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
     &     val_int=(/ndef/))
      call set_arg('DENS',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/2/))
      call set_arg('DENS',DEF_OP_FROM_OCC,'OCC',ndef*2,tgt_info,
     &     val_occ=occ_def(1:ngastp,1:2,1:ndef*2))
      call set_arg('DENS',DEF_OP_FROM_OCC,'FORMAL',1,tgt_info,
c     &     val_int=(/orb_info%nactel+1/))
     &     val_int=(/max(maxcum,ioff+2*maxexc-3)+1/))

      ! define operator adjoint to reduced densities
      call add_target('DENS_dag',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
c      do iv = 1, 2 + maxexc*(maxtop+1) - ioff
      if (maxtop.le.0) then
        nv = maxcum
      else
        nv = min(ioff + (maxexc-1)*(maxtop+1),orb_info%nactel)
      end if
      do iv = 1, nv
            ndef = ndef + 1
            occ_def(IVALE,1,ndef) = iv
            occ_def(IVALE,2,ndef) = iv
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('DENS_dag',ttype_op,DEF_OP_FROM_OCC,
     &              'DENS_dag',1,1,
     &              parameters,2,tgt_info)

      ! define precursor for expression for reduced densities
      call add_target('preDENS',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('preDENS',ttype_op,DEF_HAMILTONIAN,'preDENS',
     &              1,1,parameters,1,tgt_info)

      ! define cumulants
      call add_target('CUM',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 1, maxcum
            ndef = ndef + 1
            occ_def(IVALE,1,ndef*2) = iv
            occ_def(IVALE,2,ndef*2-1) = iv
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('CUM',ttype_op,DEF_OP_FROM_OCC,
     &              'CUM',1,1,
     &              parameters,2,tgt_info)

      ! define hole density operator
      call add_target('HOLE',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 1, 2
        ndef = ndef + 1
        occ_def(IVALE,1,ndef) = iv
        occ_def(IVALE,2,ndef) = iv
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('HOLE',ttype_op,DEF_OP_FROM_OCC,
     &              'HOLE',1,1,
     &              parameters,2,tgt_info)

      ! define central reduced densities
      call add_target('CENT',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 2, maxcum
            ndef = ndef + 1
            occ_def(IVALE,1,ndef*2) = iv
            occ_def(IVALE,2,ndef*2-1) = iv
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('CENT',ttype_op,DEF_OP_FROM_OCC,
     &              'CENT',1,1,
     &              parameters,2,tgt_info)

      ! Intermediates for cumulant approximation
      if (maxcum.gt.0.and.cum_appr_mode.ge.2) then
        call add_target2('D_INT01',.false.,tgt_info)
        occ_def = 0
        ndef = 1
        occ_def(IVALE,1,ndef*2) = 2
        occ_def(IVALE,2,ndef*2-1) = 2
        call set_rule2('D_INT01',DEF_OP_FROM_OCC,tgt_info)
        call set_arg('D_INT01',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &       val_label=(/'D_INT01'/))
        call set_arg('D_INT01',DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
     &       val_int=(/ndef/))
        call set_arg('D_INT01',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &       val_int=(/2/))
        call set_arg('D_INT01',DEF_OP_FROM_OCC,'OCC',ndef*2,tgt_info,
     &       val_occ=occ_def(1:ngastp,1:2,1:ndef*2))
        op_dint = 'D_INTXX'
        do idef = 2, 12
          write(op_dint(6:7),'(i2.2)') idef
          call add_target2(op_dint,.false.,tgt_info)
          call set_dependency(op_dint,'D_INT01',tgt_info)
          call set_rule2(op_dint,CLONE_OP,tgt_info)
          call set_arg(op_dint,CLONE_OP,'LABEL',1,tgt_info,
     &         val_label=(/op_dint/))
          call set_arg(op_dint,CLONE_OP,'TEMPLATE',1,tgt_info,
     &         val_label=(/'D_INT01'/))
        end do
        call add_target2('D_INT13',.false.,tgt_info)
        occ_def = 0
        ndef = 1
        occ_def(IVALE,1,ndef*2) = 4
        occ_def(IVALE,2,ndef*2-1) = 4
        call set_rule2('D_INT13',DEF_OP_FROM_OCC,tgt_info)
        call set_arg('D_INT13',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &       val_label=(/'D_INT13'/))
        call set_arg('D_INT13',DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
     &       val_int=(/ndef/))
        call set_arg('D_INT13',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &       val_int=(/2/))
        call set_arg('D_INT13',DEF_OP_FROM_OCC,'OCC',ndef*2,tgt_info,
     &       val_occ=occ_def(1:ngastp,1:2,1:ndef*2))
      end if
*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*

      ! precursor for reduced densities
      call add_target2('F_preDENS0',.false.,tgt_info)
      call set_dependency('F_preDENS0','preDENS',tgt_info)
      call set_dependency('F_preDENS0','DENS_dag',tgt_info)
      call set_dependency('F_preDENS0','C0',tgt_info)
      call set_rule2('F_preDENS0',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_preDENS0',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_preDENS0'/))
      call set_arg('F_preDENS0',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'preDENS'/))
      call set_arg('F_preDENS0',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &     tgt_info,
     &     val_label=(/'C0^+    ','DENS_dag','C0      '/))
      call set_arg('F_preDENS0',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &     val_int=(/2,3,4/))

      ! reduced densities
      call add_target2('F_DENS0',.false.,tgt_info)
      call set_dependency('F_DENS0','F_preDENS0',tgt_info)
      call set_dependency('F_DENS0','DENS',tgt_info)
      call set_rule2('F_DENS0',DERIVATIVE,tgt_info)
      call set_arg('F_DENS0',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_DENS0'/))
      call set_arg('F_DENS0',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_preDENS0'/))
      call set_arg('F_DENS0',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'DENS'/))
      call set_arg('F_DENS0',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'DENS_dag'/))
c dbg
c      call set_rule2('F_DENS0',PRINT_FORMULA,tgt_info)
c      call set_arg('F_DENS0',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_DENS0'/))
c dbgend

      ! precursor for density expression in terms of cumulants
      call add_target2('F_preDENS',.false.,tgt_info)
      call set_dependency('F_preDENS','preDENS',tgt_info)
      call set_dependency('F_preDENS','DENS_dag',tgt_info)
      call set_dependency('F_preDENS','CUM',tgt_info)
      call set_rule2('F_preDENS',DEF_CUMULANTS,tgt_info)
      call set_arg('F_preDENS',DEF_CUMULANTS,'LABEL',1,tgt_info,
     &     val_label=(/'F_preDENS'/))
      call set_arg('F_preDENS',DEF_CUMULANTS,'OP_RES',1,tgt_info,
     &     val_label=(/'preDENS'/))
      call set_arg('F_preDENS',DEF_CUMULANTS,'OPERATORS',2,tgt_info,
     &     val_label=(/'DENS_dag','CUM     '/))
c      call set_arg('F_preDENS',DEF_CUMULANTS,'LEVEL',1,tgt_info,
c     &     val_int=(/4+2*maxexc/))
      if (maxtop.gt.0.and.cum_appr_mode.gt.0) then
        ! Express high-rank RDMs through low-rank RDMs
        ! Substract low-rank RDMs (no substitution needed)
        call set_rule2('F_preDENS',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_preDENS',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_preDENS'/))
        call set_arg('F_preDENS',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'preDENS'/))
        call set_arg('F_preDENS',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'DENS    ','DENS_dag','DENS    '/))
        call set_arg('F_preDENS',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &       val_int=(/2,3,2/))
        call set_arg('F_preDENS',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
     &       val_int=(/maxcum,-1,maxcum/))
        call set_arg('F_preDENS',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &       val_int=(/1/))
        call set_arg('F_preDENS',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &       val_int=(/1,3/))
        call set_arg('F_preDENS',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &       val_rl8=(/-1d0/))
        call set_arg('F_preDENS',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
        ! Resubstitute reduced density matrices
        call set_rule2('F_preDENS',EXPAND,tgt_info)
        call set_dependency('F_preDENS','F_CUM',tgt_info)
        call set_arg('F_preDENS',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_preDENS'/))
        call set_arg('F_preDENS',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_preDENS'/))
        call set_arg('F_preDENS',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'F_CUM'/))
        ! delete terms with factor zero
        call set_rule2('F_preDENS',SELECT_SPECIAL,tgt_info)
        call set_arg('F_preDENS',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_preDENS'/))
        call set_arg('F_preDENS',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_preDENS'/))
        call set_arg('F_preDENS',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='NONZERO')
      end if

      ! density expression in terms of cumulants
      call add_target2('F_DENS',.false.,tgt_info)
      call set_dependency('F_DENS','F_preDENS',tgt_info)
      call set_dependency('F_DENS','DENS',tgt_info)
      call set_rule2('F_DENS',DERIVATIVE,tgt_info)
      call set_arg('F_DENS',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_DENS'/))
      call set_arg('F_DENS',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_preDENS'/))
      call set_arg('F_DENS',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'DENS'/))
      call set_arg('F_DENS',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'DENS_dag'/))
c dbg
c      call set_rule2('F_DENS',KEEP_TERMS,tgt_info)
c      call set_arg('F_DENS',KEEP_TERMS,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_DENS'/))
c      call set_arg('F_DENS',KEEP_TERMS,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_DENS'/))
c      call set_arg('F_DENS',KEEP_TERMS,'TERMS',1,tgt_info,
c     &     val_int=(/12/))
c dbgend
      call set_rule2('F_DENS',PRINT_FORMULA,tgt_info)
      call set_arg('F_DENS',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_DENS'/))

      ! precursor for cumulant expression in terms of densities
      call add_target2('F_preCUM',.false.,tgt_info)
      call set_dependency('F_preCUM','preDENS',tgt_info)
      call set_dependency('F_preCUM','DENS_dag',tgt_info)
      call set_rule2('F_preCUM',DEF_CUMULANTS,tgt_info)
      call set_arg('F_preCUM',DEF_CUMULANTS,'LABEL',1,tgt_info,
     &     val_label=(/'F_preCUM'/))
      call set_arg('F_preCUM',DEF_CUMULANTS,'OP_RES',1,tgt_info,
     &     val_label=(/'preDENS'/))
      if (maxtop.le.0) then
        call set_dependency('F_preCUM','CENT',tgt_info)
        call set_arg('F_preCUM',DEF_CUMULANTS,'OPERATORS',2,tgt_info,
     &       val_label=(/'DENS_dag','CENT    '/))
        call set_arg('F_preCUM',DEF_CUMULANTS,'MODE',1,tgt_info,
     &       val_str='CUMULANT(CENT)')
        ! exception: one-particle cumulant = one-particle RD
        call set_rule2('F_preCUM',EXPAND_OP_PRODUCT,tgt_info)
        call set_dependency('F_preCUM','DENS',tgt_info)
        call set_arg('F_preCUM',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_preCUM'/))
        call set_arg('F_preCUM',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'preDENS'/))
        call set_arg('F_preCUM',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'DENS    ','DENS_dag','DENS    '/))
        call set_arg('F_preCUM',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &       val_int=(/1,2,1/))
        call set_arg('F_preCUM',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &       val_int=(/1/))
        call set_arg('F_preCUM',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &       val_int=(/1,3/))
        call set_arg('F_preCUM',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
     &       val_int=(/1,1,1/))
        call set_arg('F_preCUM',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
      else
        call set_dependency('F_preCUM','DENS',tgt_info)
        call set_arg('F_preCUM',DEF_CUMULANTS,'OPERATORS',2,tgt_info,
     &       val_label=(/'DENS_dag','DENS    '/))
        call set_arg('F_preCUM',DEF_CUMULANTS,'MODE',1,tgt_info,
     &       val_str='CUMULANT')
      end if

      ! cumulant expression in terms of densities
      call add_target2('F_CUM',.false.,tgt_info)
      call set_dependency('F_CUM','F_preCUM',tgt_info)
      call set_dependency('F_CUM','CUM',tgt_info)
      call set_rule2('F_CUM',DERIVATIVE,tgt_info)
      call set_arg('F_CUM',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_CUM'/))
      call set_arg('F_CUM',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_preCUM'/))
      call set_arg('F_CUM',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'CUM'/))
      call set_arg('F_CUM',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'DENS_dag'/))
      call set_rule2('F_CUM',PRINT_FORMULA,tgt_info)
      call set_arg('F_CUM',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_CUM'/))

c      ! define one-particle hole density
c      call add_target2('F_HOLE_1',.false.,tgt_info)
c      call set_dependency('F_HOLE_1','HOLE',tgt_info)
c      call set_dependency('F_HOLE_1','CUM',tgt_info)
c      call set_dependency('F_HOLE_1','1',tgt_info)
c      call set_rule2('F_HOLE_1',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE_1'/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'OPERATORS',3,
c     &     tgt_info,
c     &     val_label=(/'HOLE','1','HOLE'/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
c     &     val_int=(/1,2,1/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
c     &     val_int=(/1,-1,1/))
c      call set_rule2('F_HOLE_1',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE_1'/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'OPERATORS',4,
c     &     tgt_info,
c     &     val_label=(/'HOLE','CUM','CUM','HOLE'/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
c     &     val_int=(/1,2,2,1/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
c     &     val_int=(/2,3/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
c     &     val_int=(/2,1,1,2/))
c      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      call set_rule2('F_HOLE_1',PRINT_FORMULA,tgt_info)
c      call set_arg('F_HOLE_1',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE_1'/))
c
c      ! define two-particle hole density
c      call add_target2('F_HOLE_2',.false.,tgt_info)
c      call set_dependency('F_HOLE_2','HOLE',tgt_info)
c      call set_dependency('F_HOLE_2','CUM',tgt_info)
c      call set_dependency('F_HOLE_2','1',tgt_info)
c      call set_rule2('F_HOLE_2',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE_2'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OPERATORS',3,
c     &     tgt_info,
c     &     val_label=(/'HOLE','1','HOLE'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
c     &     val_int=(/1,2,1/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
c     &     val_int=(/2,1,2/))
c      call set_rule2('F_HOLE_2',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE_2'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OPERATORS',5,
c     &     tgt_info,
c     &     val_label=(/'HOLE','CUM','1','CUM','HOLE'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
c     &     val_int=(/1,2,3,2,1/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
c     &     val_int=(/3/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'AVOID',6,tgt_info,
c     &     val_int=(/2,3,2,4,3,4/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'BLK_MIN',5,tgt_info,
c     &     val_int=(/2,1,2,1,2/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      call set_rule2('F_HOLE_2',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE_2'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OPERATORS',6,
c     &     tgt_info,
c     &     val_label=(/'HOLE','CUM','CUM','CUM','CUM','HOLE'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'IDX_SV',6,tgt_info,
c     &     val_int=(/1,2,3,3,2,1/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
c     &     val_int=(/4/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'AVOID',8,tgt_info,
c     &     val_int=(/2,4,2,5,3,4,3,5/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
c     &     val_rl8=(/2d0/))
c      call set_rule2('F_HOLE_2',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE_2'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OPERATORS',4,
c     &     tgt_info,
c     &     val_label=(/'HOLE','CUM','CUM','HOLE'/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
c     &     val_int=(/1,2,2,1/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
c     &     val_int=(/2,3/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'BLK_MIN',4,tgt_info,
c     &     val_int=(/2,2,2,2/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
c     &     val_int=(/2,2,2,2/))
c      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      call set_rule2('F_HOLE_2',PRINT_FORMULA,tgt_info)
c      call set_arg('F_HOLE_2',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE_2'/))
c
c      ! define hole density
c      call add_target2('F_HOLE',.false.,tgt_info)
c      call set_dependency('F_HOLE','HOLE',tgt_info)
c      call set_dependency('F_HOLE','CUM',tgt_info)
c      call set_dependency('F_HOLE','1',tgt_info)
c      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',3,
c     &     tgt_info,
c     &     val_label=(/'HOLE','1','HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
c     &     val_int=(/1,2,1/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
c     &     val_int=(/1,-1,1/))
c      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',4,
c     &     tgt_info,
c     &     val_label=(/'HOLE','CUM','CUM','HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
c     &     val_int=(/1,2,2,1/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
c     &     val_int=(/2,3/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
c     &     val_int=(/2,1,1,2/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',3,
c     &     tgt_info,
c     &     val_label=(/'HOLE','1','HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
c     &     val_int=(/1,2,1/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
c     &     val_int=(/2,1,2/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',5,
c     &     tgt_info,
c     &     val_label=(/'HOLE','CUM','1','CUM','HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
c     &     val_int=(/1,2,3,2,1/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
c     &     val_int=(/3/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'AVOID',6,tgt_info,
c     &     val_int=(/2,3,2,4,3,4/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MIN',5,tgt_info,
c     &     val_int=(/2,1,2,1,2/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',6,
c     &     tgt_info,
c     &     val_label=(/'HOLE','CUM','CUM','CUM','CUM','HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',6,tgt_info,
c     &     val_int=(/1,2,3,3,2,1/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
c     &     val_int=(/4/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'AVOID',8,tgt_info,
c     &     val_int=(/2,4,2,5,3,4,3,5/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
c     &     val_rl8=(/2d0/))
c      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',4,
c     &     tgt_info,
c     &     val_label=(/'HOLE','CUM','CUM','HOLE'/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
c     &     val_int=(/1,2,2,1/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
c     &     val_int=(/2,3/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MIN',4,tgt_info,
c     &     val_int=(/2,2,2,2/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
c     &     val_int=(/2,2,2,2/))
c      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      call set_rule2('F_HOLE',PRINT_FORMULA,tgt_info)
c      call set_arg('F_HOLE',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HOLE'/))

      ! precursor for central density expression in terms of densities
      call add_target2('F_preCENT',.false.,tgt_info)
      call set_dependency('F_preCENT','preDENS',tgt_info)
      call set_dependency('F_preCENT','DENS_dag',tgt_info)
      call set_dependency('F_preCENT','DENS',tgt_info)
      call set_rule2('F_preCENT',DEF_CUMULANTS,tgt_info)
      call set_arg('F_preCENT',DEF_CUMULANTS,'LABEL',1,tgt_info,
     &     val_label=(/'F_preCENT'/))
      call set_arg('F_preCENT',DEF_CUMULANTS,'OP_RES',1,tgt_info,
     &     val_label=(/'preDENS'/))
      call set_arg('F_preCENT',DEF_CUMULANTS,'OPERATORS',2,tgt_info,
     &     val_label=(/'DENS_dag','DENS    '/))
      call set_arg('F_preCENT',DEF_CUMULANTS,'MODE',1,tgt_info,
     &     val_str='CENTRAL')

      ! cumulant expression in terms of densities
      call add_target2('F_CENT',.false.,tgt_info)
      call set_dependency('F_CENT','F_preCENT',tgt_info)
      call set_dependency('F_CENT','CENT',tgt_info)
      call set_rule2('F_CENT',DERIVATIVE,tgt_info)
      call set_arg('F_CENT',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_CENT'/))
      call set_arg('F_CENT',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_preCENT'/))
      call set_arg('F_CENT',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'CENT'/))
      call set_arg('F_CENT',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'DENS_dag'/))
      call set_rule2('F_CENT',PRINT_FORMULA,tgt_info)
      call set_arg('F_CENT',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_CENT'/))

      ! Intermediates for cumulant approximation
      if (maxcum.gt.0.and.cum_appr_mode.ge.2) then
        op_dint = 'D_INTXX'
        f_dint = 'F_D_INTXX'
        do idef = 1, 13
          write(op_dint(6:7),'(i2.2)') idef
          write(f_dint(8:9),'(i2.2)') idef
          select case(idef)
          case(1,13)
            factor = -1d0
          case(2)
            factor = -(sqrt(5d0)+1d0)/(2d0*sqrt(5d0))
          case(3)
            factor = -(sqrt(5d0)-1d0)/(2d0*sqrt(5d0))
          case(4)
            factor = -(sqrt(2d0)+1d0)/sqrt(2d0)
          case(5)
            factor = -(sqrt(2d0)-1d0)/sqrt(2d0)
          case(6)
            factor = +1d0/sqrt(10d0)
          case(7)
            factor = -1d0/sqrt(10d0)
          case(8)
            factor = -1d0/3d0
          case(9)
            factor = +1d0/sqrt(6d0)
          case(10)
            factor = -1d0/sqrt(6d0)
          case(11)
            factor = -(sqrt(10d0)+1d0)/(3d0*sqrt(10d0))
          case(12)
            factor = -(sqrt(10d0)-1d0)/(3d0*sqrt(10d0))
          case default
            call quit(1,'set_gno_targets','unknown D_INT')
          end select
          call add_target2(f_dint,.false.,tgt_info)
          call set_dependency(f_dint,'DENS',tgt_info)
          call set_dependency(f_dint,op_dint,tgt_info)
          call set_rule2(f_dint,DEF_MRCC_INTM,tgt_info)
          call set_arg(f_dint,DEF_MRCC_INTM,'LABEL',1,tgt_info,
     &         val_label=(/f_dint/))
          call set_arg(f_dint,DEF_MRCC_INTM,'INTERM',1,tgt_info,
     &         val_label=(/op_dint/))
          if (idef.lt.13) then
            call set_arg(f_dint,DEF_MRCC_INTM,'OPERATORS',1,tgt_info,
     &           val_label=(/'DENS'/))
            call set_arg(f_dint,DEF_MRCC_INTM,'MODE',1,tgt_info,
     &           val_str='D_INT')
          else
            call set_dependency(f_dint,'D_INT02',tgt_info)
            call set_dependency(f_dint,'D_INT03',tgt_info)
            call set_arg(f_dint,DEF_MRCC_INTM,'OPERATORS',3,tgt_info,
     &           val_label=(/'DENS','D_INT02','D_INT03'/))
            call set_arg(f_dint,DEF_MRCC_INTM,'MODE',1,tgt_info,
     &           val_str='D_INT2')
          end if
          call set_arg(f_dint,DEF_MRCC_INTM,'FAC',1,tgt_info,
     &         val_rl8=(/factor/))
          call set_arg(f_dint,DEF_MRCC_INTM,'TITLE',1,tgt_info,
     &         val_str='Intermediate for cumulant approximation')
c dbg
c          call set_rule2(f_dint,PRINT_FORMULA,tgt_info)
c          call set_arg(f_dint,PRINT_FORMULA,'LABEL',1,tgt_info,
c     &         val_label=(/f_dint/))
c dbgend
        end do

        ! Formula for the approximated reduced density matrices
        call add_target2('F_DENS_appr',.false.,tgt_info)
        call set_dependency('F_DENS_appr','DENS',tgt_info)
        call set_rule2('F_DENS_appr',DEF_MRCC_INTM,tgt_info)
        call set_arg('F_DENS_appr',DEF_MRCC_INTM,'LABEL',1,tgt_info,
     &       val_label=(/'F_DENS_appr'/))
        call set_arg('F_DENS_appr',DEF_MRCC_INTM,'INTERM',1,tgt_info,
     &       val_label=(/'DENS'/))
        select case(maxcum)
        case(2)
          call set_dependency('F_DENS_appr','D_INT08',tgt_info)
          call set_dependency('F_DENS_appr','D_INT09',tgt_info)
          call set_dependency('F_DENS_appr','D_INT10',tgt_info)
          call set_dependency('F_DENS_appr','D_INT11',tgt_info)
          call set_dependency('F_DENS_appr','D_INT12',tgt_info)
          call set_arg('F_DENS_appr',DEF_MRCC_INTM,'OPERATORS',5,
     &         tgt_info,
     &         val_label=(/'D_INT08','D_INT09','D_INT10',
     &                     'D_INT11','D_INT12'/))
        case(3)
          call set_dependency('F_DENS_appr','D_INT04',tgt_info)
          call set_dependency('F_DENS_appr','D_INT05',tgt_info)
          call set_dependency('F_DENS_appr','D_INT06',tgt_info)
          call set_dependency('F_DENS_appr','D_INT07',tgt_info)
          call set_arg('F_DENS_appr',DEF_MRCC_INTM,'OPERATORS',4,
     &         tgt_info,
     &         val_label=(/'D_INT04','D_INT05','D_INT06','D_INT07'/))
        case(4)
          call set_dependency('F_DENS_appr','D_INT01',tgt_info)
          if (cum_appr_mode.eq.2) then
            call set_dependency('F_DENS_appr','D_INT02',tgt_info)
            call set_dependency('F_DENS_appr','D_INT03',tgt_info)
            call set_arg('F_DENS_appr',DEF_MRCC_INTM,'OPERATORS',3,
     &           tgt_info,
     &           val_label=(/'D_INT01','D_INT02','D_INT03'/))
          else
            call set_dependency('F_DENS_appr','D_INT13',tgt_info)
            call set_arg('F_DENS_appr',DEF_MRCC_INTM,'OPERATORS',3,
     &           tgt_info,
     &           val_label=(/'D_INT01','D_INT13','-------'/))
          end if
        case default
        end select
        call set_arg('F_DENS_appr',DEF_MRCC_INTM,'MAXCOM',1,
     &       tgt_info,val_int=(/maxcum/)) ! We misuse MAXCOM for maxcum
        call set_arg('F_DENS_appr',DEF_MRCC_INTM,'MODE',1,tgt_info,
     &       val_str='D_FORM')
        call set_arg('F_DENS_appr',DEF_MRCC_INTM,'TITLE',1,tgt_info,
     &       val_str='Approximated higher reduced density matrices')
c dbg
        call set_rule2('F_DENS_appr',PRINT_FORMULA,tgt_info)
        call set_arg('F_DENS_appr',PRINT_FORMULA,'LABEL',1,tgt_info,
     &       val_label=(/'F_DENS_appr'/))
c dbgend
      end if
*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! reduced densities
      call add_target2('FOPT_DENS0',.false.,tgt_info)
      call set_dependency('FOPT_DENS0','F_DENS0',tgt_info)
      call set_dependency('FOPT_DENS0','DEF_ME_DENS',tgt_info)
      call set_dependency('FOPT_DENS0','DEF_ME_C0',tgt_info)
      call set_rule2('FOPT_DENS0',OPTIMIZE,tgt_info)
      call set_arg('FOPT_DENS0',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_DENS0'/))
      call set_arg('FOPT_DENS0',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_DENS0'/))

      ! reduced densities in terms of cumulants
      call add_target2('FOPT_DENS',.false.,tgt_info)
      call set_dependency('FOPT_DENS','F_DENS',tgt_info)
      call set_dependency('FOPT_DENS','DEF_ME_DENS',tgt_info)
      call set_dependency('FOPT_DENS','DEF_ME_CUM',tgt_info)
      call set_rule2('FOPT_DENS',OPTIMIZE,tgt_info)
      call set_arg('FOPT_DENS',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_DENS'/))
      call set_arg('FOPT_DENS',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_DENS'/))

      ! cumulants in terms of reduced densities
      call add_target2('FOPT_CUM',.false.,tgt_info)
      call set_dependency('FOPT_CUM','F_CUM',tgt_info)
      call set_dependency('FOPT_CUM','DEF_ME_CENT',tgt_info)
      call set_dependency('FOPT_CUM','DEF_ME_CUM',tgt_info)
      call set_rule2('FOPT_CUM',OPTIMIZE,tgt_info)
      call set_arg('FOPT_CUM',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_CUM'/))
      call set_arg('FOPT_CUM',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_CUM'/))

c      ! hole densities
c      call add_target2('FOPT_HOLE',.false.,tgt_info)
c      call set_dependency('FOPT_HOLE','F_HOLE',tgt_info)
c      call set_dependency('FOPT_HOLE','DEF_ME_1',tgt_info)
c      call set_dependency('FOPT_HOLE','DEF_ME_CUM',tgt_info)
c      call set_dependency('FOPT_HOLE','DEF_ME_HOLE',tgt_info)
c      call set_rule2('FOPT_HOLE',OPTIMIZE,tgt_info)
c      call set_arg('FOPT_HOLE',OPTIMIZE,'LABEL_OPT',1,tgt_info,
c     &             val_label=(/'FOPT_HOLE'/))
c      call set_arg('FOPT_HOLE',OPTIMIZE,'LABELS_IN',1,tgt_info,
c     &             val_label=(/'F_HOLE'/))

      ! central reduced densities in terms of reduced densities
      call add_target2('FOPT_CENT',.false.,tgt_info)
      call set_dependency('FOPT_CENT','F_CENT',tgt_info)
      call set_dependency('FOPT_CENT','DEF_ME_DENS',tgt_info)
      call set_dependency('FOPT_CENT','DEF_ME_CENT',tgt_info)
      call set_rule2('FOPT_CENT',OPTIMIZE,tgt_info)
      call set_arg('FOPT_CENT',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_CENT'/))
      call set_arg('FOPT_CENT',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_CENT'/))

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! ME for reduced densities
      call add_target2('DEF_ME_DENS',.false.,tgt_info)
      call set_dependency('DEF_ME_DENS','DENS',tgt_info)
      call set_rule2('DEF_ME_DENS',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_DENS',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_DENS'/))
      call set_arg('DEF_ME_DENS',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'DENS'/))
      call set_arg('DEF_ME_DENS',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_DENS',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_DENS',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for cumulants
      call add_target2('DEF_ME_CUM',.false.,tgt_info)
      call set_dependency('DEF_ME_CUM','CUM',tgt_info)
      call set_rule2('DEF_ME_CUM',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_CUM',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_CUM'/))
      call set_arg('DEF_ME_CUM',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'CUM'/))
      call set_arg('DEF_ME_CUM',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_CUM',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_CUM',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for hole densities
      call add_target2('DEF_ME_HOLE',.false.,tgt_info)
      call set_dependency('DEF_ME_HOLE','HOLE',tgt_info)
      call set_rule2('DEF_ME_HOLE',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_HOLE',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_HOLE'/))
      call set_arg('DEF_ME_HOLE',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'HOLE'/))
      call set_arg('DEF_ME_HOLE',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_HOLE',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_HOLE',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for central reduced densities
      call add_target2('DEF_ME_CENT',.false.,tgt_info)
      call set_dependency('DEF_ME_CENT','CENT',tgt_info)
      call set_rule2('DEF_ME_CENT',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_CENT',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_CENT'/))
      call set_arg('DEF_ME_CENT',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'CENT'/))
      call set_arg('DEF_ME_CENT',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_CENT',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_CENT',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! Intermediates for cumulant approximation
      if (maxcum.gt.0.and.cum_appr_mode.ge.2) then
        op_dint = 'D_INTXX'
        me_dint = 'ME_D_INTXX'
        defme_dint = 'DEF_ME_D_INTXX'
        do idef = 1, 13
          write(op_dint(6:7),'(i2.2)') idef
          write(me_dint(9:10),'(i2.2)') idef
          write(defme_dint(13:14),'(i2.2)') idef
          call add_target2(defme_dint,.false.,tgt_info)
          call set_dependency(defme_dint,'DENS',tgt_info)
          call set_rule2(defme_dint,DEF_ME_LIST,tgt_info)
          call set_arg(defme_dint,DEF_ME_LIST,'LIST',1,tgt_info,
     &                 val_label=(/me_dint/))
          call set_arg(defme_dint,DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &                 val_label=(/op_dint/))
          call set_arg(defme_dint,DEF_ME_LIST,'MS',1,tgt_info,
     &                 val_int=(/0/))
          call set_arg(defme_dint,DEF_ME_LIST,'IRREP',1,tgt_info,
     &                 val_int=(/1/))
          call set_arg(defme_dint,DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &                 val_int=(/msc/))
        end do
      end if
*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! evaluate reduced densities
      call add_target2('EVAL_DENS0',.false.,tgt_info)
      call set_dependency('EVAL_DENS0','FOPT_DENS0',tgt_info)
      call set_dependency('EVAL_DENS0','EVAL_REF_S(S+1)',tgt_info)
      call set_rule2('EVAL_DENS0',EVAL,tgt_info)
      call set_arg('EVAL_DENS0',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_DENS0'/))
c dbg
c      call set_rule2('EVAL_DENS0',PRINT_MEL,tgt_info)
c      call set_arg('EVAL_DENS0',PRINT_MEL,'LIST',1,tgt_info,
c     &             val_label=(/'ME_DENS'/))
c dbgend

      ! evaluate central reduced densities
      call add_target2('EVAL_CENT',.false.,tgt_info)
      call set_dependency('EVAL_CENT','FOPT_CENT',tgt_info)
      call set_dependency('EVAL_CENT','EVAL_DENS0',tgt_info)
      call set_rule2('EVAL_CENT',EVAL,tgt_info)
      call set_arg('EVAL_CENT',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_CENT'/))
c dbg
c      call set_rule2('EVAL_CENT',PRINT_MEL,tgt_info)
c      call set_arg('EVAL_CENT',PRINT_MEL,'LIST',1,tgt_info,
c     &             val_label=(/'ME_CENT'/))
c dbgend

      ! evaluate and print out cumulants
      call add_target2('EVAL_CUM',maxtop.le.0,tgt_info)
      call set_dependency('EVAL_CUM','FOPT_CUM',tgt_info)
      call set_dependency('EVAL_CUM','EVAL_CENT',tgt_info)
      call set_rule2('EVAL_CUM',EVAL,tgt_info)
      call set_arg('EVAL_CUM',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_CUM'/))
      call set_rule2('EVAL_CUM',PRINT_MEL,tgt_info)
      call set_arg('EVAL_CUM',PRINT_MEL,'LIST',1,tgt_info,
     &             val_label=(/'ME_CUM'/))

c      ! evaluate hole densities
c      call add_target2('EVAL_HOLE',.false.,tgt_info)
c      call set_dependency('EVAL_HOLE','FOPT_HOLE',tgt_info)
c      call set_dependency('EVAL_HOLE','EVAL_CUM',tgt_info)
c      call set_rule2('EVAL_HOLE',EVAL,tgt_info)
c      call set_arg('EVAL_HOLE',EVAL,'FORM',1,tgt_info,
c     &             val_label=(/'FOPT_HOLE'/))
cc dbg
c      call set_rule2('EVAL_HOLE',PRINT_MEL,tgt_info)
c      call set_arg('EVAL_HOLE',PRINT_MEL,'LIST',1,tgt_info,
c     &             val_label=(/'ME_HOLE'/))
cc dbgend

      ! evaluate reduced densities
      call add_target2('EVAL_DENS',.false.,tgt_info)
      call set_dependency('EVAL_DENS','FOPT_DENS',tgt_info)
      call set_dependency('EVAL_DENS','EVAL_CUM',tgt_info)
      call set_rule2('EVAL_DENS',EVAL,tgt_info)
      call set_arg('EVAL_DENS',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_DENS'/))
c dbg
      call set_rule2('EVAL_DENS',PRINT_MEL,tgt_info)
      call set_arg('EVAL_DENS',PRINT_MEL,'LIST',1,tgt_info,
     &             val_label=(/'ME_DENS'/))
c dbgend

      ! set up normal ordered hamiltonian
      ! by replacing blocks of Fock operator with normal ordered FREF
      call add_target2('H_GNO',.false.,tgt_info)
      call set_dependency('H_GNO',mel_ham,tgt_info)
      call set_dependency('H_GNO','EVAL_FREF',tgt_info)
c dbg
c      call set_rule2('H_GNO',PRINT_MEL,tgt_info)
c      call set_arg('H_GNO',PRINT_MEL,'LIST',1,tgt_info,
c     &             val_label=(/mel_ham/))
c      call set_arg('H_GNO',PRINT_MEL,'COMMENT',1,tgt_info,
c     &             val_str='standard Hamiltonian:')
c dbgend
      call set_rule2('H_GNO',ADD,tgt_info)
      call set_arg('H_GNO',ADD,'LIST_SUM',1,tgt_info,
     &             val_label=(/mel_ham/))
      call set_arg('H_GNO',ADD,'LISTS',1,tgt_info,
     &             val_label=(/'ME_FREF'/))
      call set_arg('H_GNO',ADD,'FAC',1,tgt_info,
     &             val_rl8=(/1d0/))
      call set_arg('H_GNO',ADD,'REPLACE',1,tgt_info,
     &             val_log=(/.true./))
c dbg
      call set_rule2('H_GNO',PRINT_MEL,tgt_info)
      call set_arg('H_GNO',PRINT_MEL,'LIST',1,tgt_info,
     &             val_label=(/mel_ham/))
      call set_arg('H_GNO',PRINT_MEL,'COMMENT',1,tgt_info,
     &             val_str='Normal ordered Hamiltonian:')
c dbgend

      return
      end
