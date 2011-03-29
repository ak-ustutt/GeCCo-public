      subroutine set_keywords()

c      use parse_input
      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'

      integer, parameter ::
     &     ntest = 00

      integer ::
     &     iprint

      iprint = max(ntest,iprlvl)

      call keyword_init()

      call keyword_add('general')
      call argument_add('title',context='general',type=vtyp_str,
     &     len=256)
      call argument_add('memmax',context='general',type=vtyp_int,
     &     len=1,idef=(/50 000 000/))
      call argument_add('print',context='general',type=vtyp_int,
     &     len=1,idef=(/3/))
      call argument_add('form_test','general',type=vtyp_log,
     &     ldef=(/.false./))
      call argument_add('statistics','general',type=vtyp_log,
     &     ldef=(/.true./))
      call argument_add('da_block',context='general',type=vtyp_int,
     &     len=1,idef=(/32/))

      call keyword_add('orb_space')
      call keyword_add('shell',context='orb_space')
      call argument_add('def','orb_space.shell',type=vtyp_int,len=8,
     &     idef=(/-1,-1,-1,-1,-1,-1,-1,-1/))
      call argument_add('type','orb_space.shell',type=vtyp_str,len=8)
      call argument_add('nfreeze','orb_space.shell',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('nactel','orb_space.shell',type=vtyp_int,
     &     idef=(/-1/))
      call keyword_add('open_shells',context='orb_space')
      call argument_add('treat','orb_space.open_shells',
     &     type=vtyp_str,len=4,cdef=(/'p','h',' ',' '/))

      call keyword_add('method',required=.true.)
      call keyword_add('MP',context='method')
      call argument_add('level','method.MP',type=vtyp_int,idef=(/2/))

      call keyword_add('CC',context='method')
      call argument_add('maxexc','method.CC',type=vtyp_int,idef=(/2/))
      call argument_add('minexc','method.CC',type=vtyp_int,idef=(/1/))
      call argument_add('truncate','method.CC',type=vtyp_str,
     &     len=8,cdef=(/'n','o',' ',' ',' ',' ',' ',' '/))
      call argument_add('T1ext','method.CC',type=vtyp_int,idef=(/0/))
      call argument_add('H0_T1ext','method.CC',type=vtyp_int,idef=(/0/))
      call argument_add('H0d','method.CC',type=vtyp_log,
     &                                                ldef=(/.false./))

      call keyword_add('CCPT',context='method')
      call argument_add('maxexc','method.CCPT',type=vtyp_int,idef=(/3/))
      call argument_add('extern','method.CCPT',type=vtyp_int,idef=(/0/))
      call argument_add('GBC','method.CCPT',type=vtyp_log,
     &                                                 ldef=(/.false./))
      call argument_add('R2_coupling','method.CCPT',type=vtyp_log,
     &                                                 ldef=(/.true./))
      call argument_add('R2R2','method.CCPT',type=vtyp_log,
     &                                                ldef=(/.false./))
      call argument_add('RT2T2','method.CCPT',type=vtyp_log,
     &                                                ldef=(/.false./))
      call argument_add('hh_scatter','method.CCPT',type=vtyp_log,
     &                                                ldef=(/.false./))
      call argument_add('screen','method.CCPT',type=vtyp_log,
     &                                                ldef=(/.false./))

      call keyword_add('ECC',context='method')
      call argument_add('maxexc','method.ECC',type=vtyp_int,idef=(/2/))
      call argument_add('minexc','method.ECC',type=vtyp_int,idef=(/1/))
      call argument_add('truncate','method.ECC',type=vtyp_str,
     &     len=8,cdef=(/'n','o',' ',' ',' ',' ',' ',' '/))
      call argument_add('T1ext','method.ECC',type=vtyp_int,idef=(/0/))
      call argument_add('H0_T1ext','method.ECC',
     &                                       type=vtyp_int,idef=(/0/))
      call argument_add('H0d','method.ECC',type=vtyp_log,
     &                                                ldef=(/.false./))

      call keyword_add('R12',context='method')
      call argument_add('ansatz','method.R12',type=vtyp_int,idef=(/1/))
      call argument_add('maxexc','method.R12',type=vtyp_int,idef=(/2/))
      call argument_add('minexc','method.R12',type=vtyp_int,idef=(/2/))
      call argument_add('min_tp','method.R12',type=vtyp_int,idef=(/1/))
      call argument_add('min_tpp','method.R12',type=vtyp_int,idef=(/2/))
      call argument_add('RGRcouple','method.R12',
     &                                         type=vtyp_int,idef=(/0/))
      call argument_add('T1ext','method.R12',type=vtyp_int,idef=(/0/))
      call argument_add('H0_T1ext','method.R12',
     &                                        type=vtyp_int,idef=(/-1/))
      call argument_add('H0d','method.R12',type=vtyp_log,
     &                                                ldef=(/.false./))
      call argument_add('approx','method.R12',type=vtyp_str,len=8,
     &     cdef=(/'A',' ',' ',' ',' ',' ',' ',' '/))
      call argument_add('F_appr','method.R12',type=vtyp_str,len=8,
     &     cdef=(/'n','o','n','e',' ',' ',' ',' '/))
      call argument_add('K_appr','method.R12',type=vtyp_str,len=8,
     &     cdef=(/'n','o','n','e',' ',' ',' ',' '/))
      call argument_add('Z_appr','method.R12',type=vtyp_str,len=8,
     &     cdef=(/'n','o','n','e',' ',' ',' ',' '/))
c     &     cdef=(/'J','2','K','1',' ',' ',' ',' '/))
      call argument_add('Z2_appr','method.R12',type=vtyp_str,len=8,
     &     cdef=(/'a','s','-','Z',' ',' ',' ',' '/))
c     &     cdef=(/'J','1','K','1',' ',' ',' ',' '/))
      call argument_add('fixed','method.R12',type=vtyp_log,
     &     ldef=(/.false./))
      call argument_add('fix_new','method.R12',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('extend','method.R12',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('r12op','method.R12',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('pz_eval','method.R12',type=vtyp_log,
     &     ldef=(/.false./))
      call argument_add('trunc','method.R12',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('xsp1','method.R12',type=vtyp_log,
     &     ldef=(/.false./))
      call argument_add('f12x','method.R12',type=vtyp_str,len=8,
     &     cdef=(/' ',' ',' ',' ',' ',' ',' ',' '/))
      call argument_add('screen','method.R12',type=vtyp_log,
     &     ldef=(/.false./))
      call argument_add('vring','method.R12',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('use_CS','method.R12',type=vtyp_log,
     &     ldef=(/.false./))
      call argument_add('pert','method.R12',type=vtyp_log,
     &     ldef=(/.false./))
      call argument_add('opt','method.R12',type=vtyp_log,
     &     ldef=(/.false./))

      ! special keywords for multireference wave functions
      call keyword_add('MR',context='method')
      call argument_add('cminh','method.MR',type=vtyp_int,
     &                  idef=(/0/))  ! min. number of holes for wf
      call argument_add('cmaxh','method.MR',type=vtyp_int,
     &                  idef=(/-1/)) ! max. number of holes for wf
      call argument_add('cminp','method.MR',type=vtyp_int,
     &                  idef=(/0/))  ! min. number of particles for wf
      call argument_add('cmaxp','method.MR',type=vtyp_int,
     &                  idef=(/-1/)) ! max. number of particles for wf
      call argument_add('cmaxexc','method.MR',type=vtyp_int,
     &                  idef=(/0/))  ! max. excitation for wf
      call argument_add('ciroot','method.MR',type=vtyp_int,
     &                  idef=(/1/))  ! root to be taken for wf
      call argument_add('minh','method.MR',type=vtyp_int,
     &                  idef=(/0/))  ! min. number of holes
      call argument_add('maxh','method.MR',type=vtyp_int,
     &                  idef=(/-1/)) ! max. number of holes
      call argument_add('minp','method.MR',type=vtyp_int,
     &                  idef=(/0/))  ! min. number of particles
      call argument_add('maxp','method.MR',type=vtyp_int,
     &                  idef=(/-1/)) ! max. number of particles
      call argument_add('maxv','method.MR',type=vtyp_int,
     &                  idef=(/-1/)) ! max. number of valence ops
      call argument_add('maxvv','method.MR',type=vtyp_int,
     &                  idef=(/-1/)) ! max. number of val-val exc.
      call argument_add('maxexc','method.MR',type=vtyp_int,
     &                  idef=(/2/))  ! max. excitation
      call argument_add('minexc','method.MR',type=vtyp_int,
     &                  idef=(/1/))  ! min. excitation
      call argument_add('pure_vv','method.MR',type=vtyp_log,
     &                  ldef=(/.true./)) ! pure act.-act. excitations
      call argument_add('excrestr','method.MR',type=vtyp_int,len=6,
     &                  idef=(/-1,-1,-1,-1,-1,-1/)) ! restr. exc. in T
      call argument_add('GNO','method.MR',type=vtyp_int,
     &                  idef=(/0/))  ! 1 for generalized normal order
      call argument_add('maxcum','method.MR',type=vtyp_int,
     &                  idef=(/4/))  ! max. cumulant rank
      call argument_add('oldref','method.MR',type=vtyp_log,
     &                  ldef=(/.false./)) ! use existing CASSCF coeff.
      call argument_add('cheap_prc','method.MR',type=vtyp_log,
     &                  ldef=(/.false./)) ! try cheap preconditioner
      call argument_add('project','method.MR',type=vtyp_log,
     &     ldef=(/.false./)) ! project out singles from doubles a.s.o.

      call keyword_add('MRCI',context='method')
      call argument_add('nroots','method.MRCI',type=vtyp_int,
     &                  idef=(/1/))  ! number of roots (for ic-MRCI)

      call keyword_add('MRCC',context='method')
      call argument_add('maxcom_res','method.MRCC',type=vtyp_int,
     &     idef=(/4/))
      call argument_add('maxcom_en','method.MRCC',type=vtyp_int,
     &     idef=(/4/))
      call argument_add('maxtt','method.MRCC',type=vtyp_int,
     &     idef=(/-1/))
      call argument_add('G_level','method.MRCC',type=vtyp_int,
     &     idef=(/-1/)) ! max. power in e^(-T) for sim. trans.
      call argument_add('Op_eqs','method.MRCC',type=vtyp_log,
     &     ldef=(/.false./))

      ! Truncations (obsolete)
      call keyword_add('truncate',context='method')
      call argument_add('trunc_type','method.truncate',
     &     type=vtyp_int,idef=(/0/))

      call keyword_add('calculate')
      ! internal tests
      call keyword_add('check',context='calculate')
      call keyword_add('formulae',context='calculate.check')
      call argument_add('contr_test','calculate.check',
     &     type=vtyp_int,idef=(/1/))
      call argument_add('algebra_test','calculate.check',
     &     type=vtyp_int,idef=(/1/))

      ! general
      call keyword_add('solve',context='calculate')
      call argument_add('maxiter','calculate.solve',type=vtyp_int,
     &     idef=(/30/))
      call argument_add('maxmic','calculate.solve',type=vtyp_int,
     &     idef=(/20/))
      call argument_add('maxsub','calculate.solve',type=vtyp_int,
     &     idef=(/8/))
      call argument_add('check_incore','calculate.solve',type=vtyp_int,
     &     idef=(/-1/))
      call argument_add('conv','calculate.solve',type=vtyp_rl8,
     &     xdef=(/1d-6/))

      ! specials for: non-linear, linear, eigenvalue
      call keyword_add('non_linear',context='calculate.solve')
      call argument_add('maxiter','calculate.solve.non_linear',
     &     type=vtyp_int,
     &     idef=(/30/))
      call argument_add('maxmic','calculate.solve.non_linear',
     &     type=vtyp_int,
     &     idef=(/20/))
      call argument_add('maxsub','calculate.solve.non_linear',
     &     type=vtyp_int,
     &     idef=(/8/))
      call argument_add('conv','calculate.solve.non_linear',
     &     type=vtyp_rl8,
     &     xdef=(/1d-6/))
      call argument_add('tr_ini','calculate.solve.non_linear',
     &     type=vtyp_rl8,
     &     xdef=(/1.0d0/))
      call argument_add('method','calculate.solve.non_linear',
     &     type=vtyp_str,len=8,
     &     cdef=(/'d','i','i','s',' ',' ',' ',' '/))
      call argument_add('optref','calculate.solve.non_linear',
     &     type=vtyp_int,
     &     idef=(/0/)) ! optimize reference fct.
      call argument_add('update_prc','calculate.solve.non_linear',
     &     type=vtyp_log,
     &     ldef=(/.true./)) ! update precond. when metric is updated
      call argument_add('preopt','calculate.solve.non_linear',
     &     type=vtyp_log,
     &     ldef=(/.false./)) ! first one optimization with fixed metric

      call keyword_add('linear',context='calculate.solve')
      call argument_add('maxiter','calculate.solve.linear',
     &     type=vtyp_int,
     &     idef=(/30/))
      call argument_add('maxsub','calculate.solve.linear',
     &     type=vtyp_int,
     &     idef=(/8/))
      call argument_add('conv','calculate.solve.linear',
     &     type=vtyp_rl8,
     &     xdef=(/1d-6/))
      call argument_add('method','calculate.solve.linear',
     &     type=vtyp_str,len=8,
     &     cdef=(/'s','u','b','s','p','a','c','e'/))

      call keyword_add('eigen',context='calculate.solve')
      call argument_add('maxiter','calculate.solve.eigen',
     &     type=vtyp_int,
     &     idef=(/30/))
      call argument_add('maxsub','calculate.solve.eigen',
     &     type=vtyp_int,
     &     idef=(/8/))
      call argument_add('conv','calculate.solve.eigen',
     &     type=vtyp_rl8,
     &     xdef=(/1d-6/))
      call argument_add('method','calculate.solve.eigen',
     &     type=vtyp_str,len=8,
     &     cdef=(/'d','a','v','i','d','s','o','n'/))

      call keyword_add('CC_solve_tbar',context='calculate')
      call keyword_add('CC_solve_sim',context='calculate')
      call keyword_add('properties',context='calculate')

      call keyword_add('skip_E',context='calculate')

      call keyword_add('check_S',context='calculate')
      call argument_add('sym','calculate.check_S',
     &     type=vtyp_int,len=8,
     &     idef=(/1,0,0,0,0,0,0,0/))
      call argument_add('msc','calculate.check_S',
     &     type=vtyp_int,len=1,idef=(/0/))

      call keyword_add('excitation',context='calculate')
      call argument_add('sym','calculate.excitation',
     &     type=vtyp_int,len=8,
     &     idef=(/1,0,0,0,0,0,0,0/))
      call argument_add('msc','calculate.excitation',
     &     type=vtyp_int,len=1,idef=(/0/))
      call keyword_add('normalize',context='calculate.excitation')
      call keyword_add('analyze',context='calculate.excitation')

      call keyword_add('ionization',context='calculate')
      call argument_add('sym','calculate.ionization',
     &     type=vtyp_int,len=8,
     &     idef=(/1,0,0,0,0,0,0,0/))

      call keyword_add('routes',context='calculate')
      call argument_add('schedule','calculate.routes',type=vtyp_int,
c     &     idef=(/0/))
     &     idef=(/1/))
      call argument_add('contract','calculate.routes',type=vtyp_int,
     &     idef=(/3/))
      call argument_add('str_block','calculate.routes',type=vtyp_int,
     &     idef=(/200/))
      call argument_add('cnt_block','calculate.routes',type=vtyp_int,
     &     idef=(/-1/))
      call argument_add('force_batching',
     &     'calculate.routes',type=vtyp_int,
     &     idef=(/-1/))
      call argument_add('force_ooc_sort',
     &     'calculate.routes',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('maxbranch','calculate.routes',type=vtyp_int,
     &     idef=(/-1/))
      call argument_add('use_tr','calculate.routes',type=vtyp_log,
     &     ldef=(/.true./))
      call argument_add('simtraf','calculate.routes',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('sv_thresh','calculate.routes',type=vtyp_rl8,
     &     xdef=(/1d-14/))
      call argument_add('sv_fix','calculate.routes',type=vtyp_log,
     &     ldef=(/.false./))

      ! special keywords for response theory
      call keyword_add('response',context='method')
      call argument_add('order','method.response',type=vtyp_int,
     &                  idef=(/0/)) ! perturbation order
      call argument_add('pert','method.response',
     &                  type=vtyp_str,len=1,cdef=(/'d'/)) ! pert. op.
      call argument_add('comp','method.response',
     &                  type=vtyp_str,len=1,cdef=(/'Z'/)) ! cartesian components
      call argument_add('freq','method.response',
     &                  type=vtyp_rl8,xdef=(/0d0/)) ! frequencies
      call argument_add('BX_3C','method.response',
     &     type=vtyp_log,ldef=(/.true./)) ! treat BX intermed. as in approx.3C
      call argument_add('rules','method.response',
     &     type=vtyp_log,ldef=(/.true./)) ! use 2n+1 / 2n+2 rules
      call argument_add('restart','method.response',
     &     type=vtyp_int,idef=(/0/)) ! restart calc. at given prop. order

      call keyword_add('experimental',context='calculate')
      ! set additional experimental keyword in this subroutine:
      call set_experimental_keywords()

      if (iprint.ge.50)
     &     call show_keywords(luout)

      return
      end
