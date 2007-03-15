*----------------------------------------------------------------------*
      subroutine do_calc(orb_info)
*----------------------------------------------------------------------*
      
      implicit none
      include 'stdunit.h'
      include 'ifc_memman.h'
      include 'def_orbinf.h'
      include 'def_operator.h'
      include 'def_operator_list.h'
      include 'def_filinf.h'
      include 'def_file_list.h'

      type(orbinf), intent(inout) ::
     &     orb_info

      type(operator_list), pointer ::
     &     op_list
      type(file_list), pointer ::
     &     form_list
      type(operator), pointer ::
     &     ops(:)
      integer ::
     &     ifree, nops, nform

      ifree = mem_setmark('do_calc')
      
      ! set up orbital info
      call set_orbinf(orb_info,.true.)

      ifree = mem_setmark('operator def')
      allocate(op_list)
      nullify(op_list%op)
      nullify(op_list%prev)
      nullify(op_list%next)
      nops = 0
      ! set up operators
      call set_operators(op_list,nops,orb_info)
      if (nops.eq.0)
     &     call quit(0,'do_calc','no operators defined?')

      ifree = mem_setmark('formula def')
      allocate(form_list)
      nullify(form_list%fhand)
      nullify(form_list%prev)
      nullify(form_list%next)
      nform = 0
      ! set up formulae
      call set_formulae(form_list,nform,op_list,nops)

      ! set up graphs
c      call set_graphs_for_ops()
      
      ! set up operator dimensions

      ! optimize the formulae

      ! loop over requests from current calculate block

        ! initialize operator matrix elements
      
        ! call the appropriate solver
      

      ! free memory allocated for operators

      ifree = mem_flushmark('do_calc')

      return
      end
