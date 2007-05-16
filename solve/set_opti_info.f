*----------------------------------------------------------------------*
      subroutine set_opti_info(opti_info,nopt,op_opt)
*----------------------------------------------------------------------*
      implicit none
      
      include 'def_optimize_info.h'
      include 'ifc_input.h'
      include 'ifc_memman.h'
      include 'def_operator.h'
      include 'def_operator_array.h'

      type(optimize_info), intent(inout) ::
     &     opti_info
      integer, intent(in) ::
     &     nopt
      type(operator_array), intent(in) ::
     &     op_opt(nopt)

      character ::
     &     str*16
      integer ::
     &     ifree, iopt

      opti_info%variational = .false.
      opti_info%linear      = .false.

      call get_argument_value('calculate.solve','maxiter',
     &     ival=opti_info%maxmacit)
      call get_argument_value('calculate.solve','maxmic',
     &     ival=opti_info%micifac)
      opti_info%maxmicit = opti_info%maxmacit*opti_info%micifac
      call get_argument_value('calculate.solve','maxsub',
     &     ival=opti_info%maxsbsp)

      call get_argument_value('calculate.solve','method',
     &     str=str)
      call uppcas(str)
      select case(trim(str(1:4)))
      case('PERT') 
        opti_info%mode_nleq = mode_nleq_pert
        opti_info%norder = 1
      case('DIIS') 
        opti_info%mode_nleq = mode_nleq_diis
        opti_info%norder = 1
      case('ASSJ','RLE ') 
        opti_info%mode_nleq = mode_nleq_assj
        opti_info%norder = 1
      case default
        call quit(0,'set_opti','invalid method: '//trim(str))
      end select
      
      ifree = mem_alloc_int (opti_info%nwfpar,nopt,'nwfpar')
      ifree = mem_alloc_real(opti_info%thrgrd,nopt,'thrgrd')

      call get_argument_value('calculate.solve','conv',
     &     xval=opti_info%thrgrd(1))
      if (nopt.gt.1)
     &     opti_info%thrgrd(2:nopt) = opti_info%thrgrd(1)

      call get_argument_value('calculate.solve','tr_ini',
     &     xval=opti_info%trini)

      opti_info%nopt = nopt
      do iopt = 1, nopt
        opti_info%nwfpar(iopt) = op_opt(iopt)%op%len_op
      end do

      return
      end
