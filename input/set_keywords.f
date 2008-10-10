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

      call keyword_add('orb_space')
      call keyword_add('shell',context='orb_space')
      call argument_add('def','orb_space.shell',type=vtyp_int,len=8)
      call argument_add('type','orb_space.shell',type=vtyp_str,len=8)
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

      call keyword_add('CCPT',context='method')
      call argument_add('maxexc','method.CCPT',type=vtyp_int,idef=(/3/))

      call keyword_add('ECC',context='method')
      call argument_add('maxexc','method.ECC',type=vtyp_int,idef=(/2/))
      call argument_add('minexc','method.ECC',type=vtyp_int,idef=(/1/))
      call argument_add('truncate','method.ECC',type=vtyp_str,
     &     len=8,cdef=(/'n','o',' ',' ',' ',' ',' ',' '/))

      call keyword_add('R12',context='method')
      call argument_add('ansatz','method.R12',type=vtyp_int,idef=(/1/))
      call argument_add('maxexc','method.R12',type=vtyp_int,idef=(/2/))
      call argument_add('minexc','method.R12',type=vtyp_int,idef=(/2/))
      call argument_add('min_tp','method.R12',type=vtyp_int,idef=(/1/))
      call argument_add('min_tpp','method.R12',type=vtyp_int,idef=(/2/))
      call argument_add('T1ext','method.R12',type=vtyp_int,idef=(/0/))
      call argument_add('approx','method.R12',type=vtyp_str,len=8,
     &     cdef=(/'A',' ',' ',' ',' ',' ',' ',' '/))
      call argument_add('F_appr','method.R12',type=vtyp_str,len=8,
     &     cdef=(/'n','o','n','e',' ',' ',' ',' '/))
      call argument_add('K_appr','method.R12',type=vtyp_str,len=8,
     &     cdef=(/'n','o','n','e',' ',' ',' ',' '/))
      call argument_add('fixed','method.R12',type=vtyp_log,
     &     ldef=(/.false./))
      call argument_add('extend','method.R12',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('r12op','method.R12',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('pz_eval','method.R12',type=vtyp_log,
     &     ldef=(/.false./))

      ! Truncations
      call keyword_add('truncate',context='method')
      call argument_add('trunc_type','method.truncate',
     &     type=vtyp_int,idef=(/0/))

      call keyword_add('calculate')
      ! internal tests
      call keyword_add('check',context='calculate')
      call keyword_add('formulae',context='calculate.check')

      ! general
      call keyword_add('solve',context='calculate')
      call argument_add('maxiter','calculate.solve',type=vtyp_int,
     &     idef=(/30/))
      call argument_add('maxmic','calculate.solve',type=vtyp_int,
     &     idef=(/20/))
      call argument_add('maxsub','calculate.solve',type=vtyp_int,
     &     idef=(/8/))
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

      call keyword_add('excitation',context='calculate')
      call argument_add('sym','calculate.excitation',
     &     type=vtyp_int,len=8,
     &     idef=(/1,0,0,0,0,0,0,0/))
      call argument_add('msc','calculate.excitation',
     &     type=vtyp_int,len=1,idef=(/0/))

      call keyword_add('ionization',context='calculate')
      call argument_add('sym','calculate.ionization',
     &     type=vtyp_int,len=8,
     &     idef=(/1,0,0,0,0,0,0,0/))

      call keyword_add('routes',context='calculate')
      call argument_add('schedule','calculate.routes',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('contract','calculate.routes',type=vtyp_int,
     &     idef=(/3/))
      call argument_add('str_block','calculate.routes',type=vtyp_int,
     &     idef=(/200/))
      call argument_add('use_tr','calculate.routes',type=vtyp_log,
     &     ldef=(/.true./))
      call argument_add('simtraf','calculate.routes',type=vtyp_int,
     &     idef=(/0/))

      call keyword_add('experimental',context='calculate')
      ! set additional experimental keyword in this subroutine:
      call set_experimental_keywords()

      if (iprint.ge.50)
     &     call show_keywords(luout)

      return
      end
