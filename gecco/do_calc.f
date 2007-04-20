*----------------------------------------------------------------------*
      subroutine do_calc(orb_info,env_type)
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
      character, intent(in) ::
     &     env_type*(*)

      type(operator_list), pointer ::
     &     op_list
      type(file_list), pointer ::
     &     form_list
      type(action_list), pointer ::
     &     act_list, current_act
      type(operator), pointer ::
     &     ops(:)
      type(filinf), pointer ::
     &     ffform(:), ffops(:)
      integer ::
     &     ifree, nops, nform, nactions
      type(strinf), pointer ::
     &     str_info

      ifree = mem_setmark('do_calc')
      
      ! set up orbital info
      call set_orbinf(orb_info,.true.)

      ifree = mem_setmark('operator_def')
      allocate(op_list)
      nullify(op_list%op)
      nullify(op_list%prev)
      nullify(op_list%next)
      nops = 0
      ! set up operators
      call set_operators(op_list,nops,orb_info)
      if (nops.eq.0)
     &     call quit(0,'do_calc','no operators defined?')

      ifree = mem_setmark('formula_def')
      allocate(form_list)
      nullify(form_list%fhand)
      nullify(form_list%prev)
      nullify(form_list%next)
      nform = 0
      ! set up (basic) formulae
      call set_formulae(form_list,nform,op_list,nops)
      if (nform.eq.0)
     &     call quit(0,'do_calc','no formulae/method defined?')

      ifree = mem_setmark('action_def')
      allocate(act_list)
      nullify(act_list%act)
      nullify(act_list%prev)
      nullify(act_list%next)
      nactions = 0
      ! set up actions
      call set_actions(act_list,nactions,
     &     form_list,nform,op_list,nops)
      if (nform.eq.0)
     &     call quit(0,'do_calc','no actions defined?')

      ! set up graphs
      ifree = mem_setmark('graph_def')
      allocate(str_info)
      call set_graphs_for_ops(str_info,op_list,nops,orb_info)
      
      ! set up operator dimensions
      ifree = mem_setmark('operator_dim')
      call set_dim_for_ops(op_list,nops,str_info,orb_info)

      ! turn linked lists into arrays
      allocate(ops(nops))
      call op_list2arr(op_list,ops,nops)
      allocate(ffform(nops))
      call file_list2arr(form_list,ffform,nform)

      ! initialize files for operator elements
      allocate(ffops(nops))
      call init_op_files(ffops,ops,nops)

      ! loop over requested actions
      current_act => act_list
      do
        if (.not.associated(current_act%act))
     &       call quit(1,'do_calc','action list is buggy')

        select case (current_act%act%action_type)
          case (iaction_import)
            ! import operator matrix elements
            call import_op_el(current_act%act%idxopdef_out(1),
     &                        current_act%act%idxopfile_out(1,1),
     &                        ffops,ops,nops,
     &                        env_type,str_info,orb_info)
          case (iaction_evaluate)
            ! evaluate a single formula expression
            call quit(1,'do_calc','action not implemented yet')
          case (iaction_solve_leq)
            ! Solve system of linear equations
            call quit(1,'do_calc','action not implemented yet')
          case (iaction_solve_nleq)
            ! Solve system of non-linear equations
            call quit(1,'do_calc','action not implemented yet')
          case (iaction_solve_evp)
            ! Solve eigenvalue problem
            call quit(1,'do_calc','action not implemented yet')
          case (iaction_solve_gevp)
            ! Solve general eigenvalue problem
            call quit(1,'do_calc','action not implemented yet')
          case default
            write(luout,*) 'action = ',current_act%act%action_type
            call quit(1,'do_calc','unknown action')
        end select

      
        ! optimize the formulae
        ! call the appropriate solver
      
      end do
        
      ! free memory allocated for operators etc.
      deallocate(ops,ffform,ffops)

      ifree = mem_flushmark('do_calc')

      return
      end
