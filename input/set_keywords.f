      subroutine set_keywords()

      use parse_input
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 100

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

      call keyword_add('method',required=.true.)
      call keyword_add('CC',context='method')
      call argument_add('maxexc','method.CC',type=vtyp_int,idef=(/2/))
      call argument_add('minexc','method.CC',type=vtyp_int,idef=(/1/))

      call keyword_add('R12',context='method')
      call argument_add('ansatz','method.R12',type=vtyp_int,idef=(/1/))
      call argument_add('triples','method.R12',type=vtyp_int,idef=(/0/))

      call keyword_add('calculate')
      call keyword_add('solve',context='calculate')
      call argument_add('maxiter','calculate.solve',type=vtyp_int,
     &     idef=(/30/))
      call argument_add('maxmic','calculate.solve',type=vtyp_int,
     &     idef=(/20/))
      call argument_add('maxsub','calculate.solve',type=vtyp_int,
     &     idef=(/8/))
      call argument_add('conv','calculate.solve',type=vtyp_rl8,
     &     xdef=(/1d-6/))
      call argument_add('tr_ini','calculate.solve',type=vtyp_rl8,
     &     xdef=(/1.0d0/))
      call argument_add('method','calculate.solve',type=vtyp_str,len=8,
     &     cdef=(/'d','i','i','s',' ',' ',' ',' '/))

      call keyword_add('routes',context='calculate')
      call argument_add('schedule','calculate.routes',type=vtyp_int,
     &     idef=(/0/))
      call argument_add('contract','calculate.routes',type=vtyp_int,
     &     idef=(/1/))
      call argument_add('simtraf','calculate.routes',type=vtyp_int,
     &     idef=(/0/))

      if (iprint.ge.50)
     &     call keyword_list(luout,keyword_root,show_args=.true.)

      return
      end
