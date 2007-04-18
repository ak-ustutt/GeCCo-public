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
      include 'def_action.h'
      include 'def_action_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'

      type(orbinf), intent(inout) ::
     &     orb_info

      type(operator_list), pointer ::
     &     op_list
      type(file_list), pointer ::
     &     form_list
      type(action_list), pointer ::
     &     act_list
      type(operator), pointer ::
     &     ops(:)
      integer ::
     &     ifree, nops, nform, nactions
      type(strinf), pointer ::
     &     str_info

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
      ! set up (basic) formulae
      call set_formulae(form_list,nform,op_list,nops)
      if (nform.eq.0)
     &     call quit(0,'do_calc','no formulae/method defined?')

      ifree = mem_setmark('action def')
      allocate(act_list)
      nullify(act_list%act)
      nullify(act_list%prev)
      nullify(act_list%next)
      nactions = 0
      ! set up actions
      call set_actions(act_list,nactions,
     &     form_list,nform,op_list,nops)

      ! set up graphs
      ifree = mem_setmark('graph def')
      allocate(str_info)
      call set_graphs_for_ops(str_info,op_list,nops,orb_info)
      
      ! set up operator dimensions


      ! loop over requested actions

        ! initialize operator matrix elements
      
        ! optimize the formulae
        ! call the appropriate solver
      

      ! free memory allocated for operators

      ifree = mem_flushmark('do_calc')

      return
      end
