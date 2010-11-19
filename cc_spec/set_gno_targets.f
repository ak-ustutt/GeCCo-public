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

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     maxtop

      integer ::
     &     ndef, occ_def(ngastp,2,60), msc, maxexc, iv, maxcum,
     &     ioff, gno
      logical ::
     &     pure_vv
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)

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
      msc = +1 ! assuming closed shell

      ! get minimum and maximum numbers of excitations, holes, particles,
      ! valence-valence excitations
      call get_argument_value('method.MR','maxexc',
     &     ival=maxexc)
      call get_argument_value('method.MR','maxcum',
     &     ival=maxcum)
      call get_argument_value('method.MR','GNO',
     &     ival=gno)
      call get_argument_value('method.MR','pure_vv',
     &     lval=pure_vv)
      ioff = 0
      if (gno.eq.0) then
        ioff = 3
        if (pure_vv) ioff = 2
      end if

c dbg
      print *,'maxcum  = ',maxcum
c dbgend

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define (simple) reduced densities
c      call add_target('DENS',ttype_op,.false.,tgt_info)
      call add_target2('DENS',.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 1, 2 + maxexc*(maxtop+1) - ioff
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
     &     val_int=(/orb_info%nactel+1/))

      ! define operator adjoint to reduced densities
      call add_target('DENS_dag',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 1, 2 + maxexc*(maxtop+1) - ioff
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
     &     val_label=(/'C0^+','DENS_dag','C0'/))
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
      call set_rule2('F_DENS0',PRINT_FORMULA,tgt_info)
      call set_arg('F_DENS0',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_DENS0'/))

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
     &     val_label=(/'DENS_dag','CUM'/))
      call set_arg('F_preDENS',DEF_CUMULANTS,'LEVEL',1,tgt_info,
     &     val_int=(/4+2*maxexc/))

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
      call set_dependency('F_preCUM','DENS',tgt_info)
      call set_rule2('F_preCUM',DEF_CUMULANTS,tgt_info)
      call set_arg('F_preCUM',DEF_CUMULANTS,'LABEL',1,tgt_info,
     &     val_label=(/'F_preCUM'/))
      call set_arg('F_preCUM',DEF_CUMULANTS,'OP_RES',1,tgt_info,
     &     val_label=(/'preDENS'/))
      call set_arg('F_preCUM',DEF_CUMULANTS,'OPERATORS',2,tgt_info,
     &     val_label=(/'DENS_dag','DENS'/))
      call set_arg('F_preCUM',DEF_CUMULANTS,'MODE',1,tgt_info,
     &     val_str='CUMULANT')

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

      ! define one-particle hole density
      call add_target2('F_HOLE_1',.false.,tgt_info)
      call set_dependency('F_HOLE_1','HOLE',tgt_info)
      call set_dependency('F_HOLE_1','CUM',tgt_info)
      call set_dependency('F_HOLE_1','1',tgt_info)
      call set_rule2('F_HOLE_1',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE_1'/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &     tgt_info,
     &     val_label=(/'HOLE','1','HOLE'/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &     val_int=(/1,2,1/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
     &     val_int=(/1,-1,1/))
      call set_rule2('F_HOLE_1',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE_1'/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'HOLE','CUM','CUM','HOLE'/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,2,1/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &     val_int=(/2,3/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
     &     val_int=(/2,1,1,2/))
      call set_arg('F_HOLE_1',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_HOLE_1',PRINT_FORMULA,tgt_info)
      call set_arg('F_HOLE_1',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE_1'/))

      ! define two-particle hole density
      call add_target2('F_HOLE_2',.false.,tgt_info)
      call set_dependency('F_HOLE_2','HOLE',tgt_info)
      call set_dependency('F_HOLE_2','CUM',tgt_info)
      call set_dependency('F_HOLE_2','1',tgt_info)
      call set_rule2('F_HOLE_2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE_2'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &     tgt_info,
     &     val_label=(/'HOLE','1','HOLE'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &     val_int=(/1,2,1/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
     &     val_int=(/2,1,2/))
      call set_rule2('F_HOLE_2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE_2'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &     tgt_info,
     &     val_label=(/'HOLE','CUM','1','CUM','HOLE'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &     val_int=(/1,2,3,2,1/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/3/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'AVOID',6,tgt_info,
     &     val_int=(/2,3,2,4,3,4/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'BLK_MIN',5,tgt_info,
     &     val_int=(/2,1,2,1,2/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_HOLE_2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE_2'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OPERATORS',6,
     &     tgt_info,
     &     val_label=(/'HOLE','CUM','CUM','CUM','CUM','HOLE'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'IDX_SV',6,tgt_info,
     &     val_int=(/1,2,3,3,2,1/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/4/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'AVOID',8,tgt_info,
     &     val_int=(/2,4,2,5,3,4,3,5/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/2d0/))
      call set_rule2('F_HOLE_2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE_2'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'HOLE','CUM','CUM','HOLE'/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,2,1/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &     val_int=(/2,3/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'BLK_MIN',4,tgt_info,
     &     val_int=(/2,2,2,2/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
     &     val_int=(/2,2,2,2/))
      call set_arg('F_HOLE_2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_HOLE_2',PRINT_FORMULA,tgt_info)
      call set_arg('F_HOLE_2',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE_2'/))

      ! define hole density
      call add_target2('F_HOLE',.false.,tgt_info)
      call set_dependency('F_HOLE','HOLE',tgt_info)
      call set_dependency('F_HOLE','CUM',tgt_info)
      call set_dependency('F_HOLE','1',tgt_info)
      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &     tgt_info,
     &     val_label=(/'HOLE','1','HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &     val_int=(/1,2,1/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
     &     val_int=(/1,-1,1/))
      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'HOLE','CUM','CUM','HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,2,1/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &     val_int=(/2,3/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
     &     val_int=(/2,1,1,2/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &     tgt_info,
     &     val_label=(/'HOLE','1','HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &     val_int=(/1,2,1/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
     &     val_int=(/2,1,2/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &     tgt_info,
     &     val_label=(/'HOLE','CUM','1','CUM','HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &     val_int=(/1,2,3,2,1/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/3/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'AVOID',6,tgt_info,
     &     val_int=(/2,3,2,4,3,4/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MIN',5,tgt_info,
     &     val_int=(/2,1,2,1,2/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',6,
     &     tgt_info,
     &     val_label=(/'HOLE','CUM','CUM','CUM','CUM','HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',6,tgt_info,
     &     val_int=(/1,2,3,3,2,1/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/4/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'AVOID',8,tgt_info,
     &     val_int=(/2,4,2,5,3,4,3,5/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/2d0/))
      call set_rule2('F_HOLE',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'HOLE','CUM','CUM','HOLE'/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,2,1/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &     val_int=(/2,3/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MIN',4,tgt_info,
     &     val_int=(/2,2,2,2/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
     &     val_int=(/2,2,2,2/))
      call set_arg('F_HOLE',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_HOLE',PRINT_FORMULA,tgt_info)
      call set_arg('F_HOLE',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_HOLE'/))

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
      call set_dependency('FOPT_CUM','DEF_ME_DENS',tgt_info)
      call set_dependency('FOPT_CUM','DEF_ME_CUM',tgt_info)
      call set_rule2('FOPT_CUM',OPTIMIZE,tgt_info)
      call set_arg('FOPT_CUM',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_CUM'/))
      call set_arg('FOPT_CUM',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_CUM'/))

      ! hole densities
      call add_target2('FOPT_HOLE',.false.,tgt_info)
      call set_dependency('FOPT_HOLE','F_HOLE',tgt_info)
      call set_dependency('FOPT_HOLE','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_HOLE','DEF_ME_CUM',tgt_info)
      call set_dependency('FOPT_HOLE','DEF_ME_HOLE',tgt_info)
      call set_rule2('FOPT_HOLE',OPTIMIZE,tgt_info)
      call set_arg('FOPT_HOLE',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_HOLE'/))
      call set_arg('FOPT_HOLE',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_HOLE'/))

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

      ! evaluate cumulants
      call add_target2('EVAL_CUM',.false.,tgt_info)
      call set_dependency('EVAL_CUM','FOPT_CUM',tgt_info)
      call set_dependency('EVAL_CUM','EVAL_DENS0',tgt_info)
      call set_rule2('EVAL_CUM',EVAL,tgt_info)
      call set_arg('EVAL_CUM',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_CUM'/))
c dbg
      call set_rule2('EVAL_CUM',PRINT_MEL,tgt_info)
      call set_arg('EVAL_CUM',PRINT_MEL,'LIST',1,tgt_info,
     &             val_label=(/'ME_CUM'/))
c dbgend

      ! evaluate hole densities
      call add_target2('EVAL_HOLE',.false.,tgt_info)
      call set_dependency('EVAL_HOLE','FOPT_HOLE',tgt_info)
      call set_dependency('EVAL_HOLE','EVAL_CUM',tgt_info)
      call set_rule2('EVAL_HOLE',EVAL,tgt_info)
      call set_arg('EVAL_HOLE',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_HOLE'/))
c dbg
      call set_rule2('EVAL_HOLE',PRINT_MEL,tgt_info)
      call set_arg('EVAL_HOLE',PRINT_MEL,'LIST',1,tgt_info,
     &             val_label=(/'ME_HOLE'/))
c dbgend

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
