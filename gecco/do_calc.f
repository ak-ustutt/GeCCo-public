*----------------------------------------------------------------------*
      subroutine do_calc(orb_info,env_type)
*----------------------------------------------------------------------*
      
      implicit none
      include 'stdunit.h'
      include 'par_globalmarks.h'
      include 'par_formnames_gen.h'
      include 'ifc_memman.h'
      include 'def_orbinf.h'
      include 'def_operator.h'
      include 'def_operator_list.h'
      include 'def_operator_array.h'
      include 'def_filinf.h'
      include 'def_file_list.h'
      include 'def_file_array.h'
      include 'def_operator_info.h'
      include 'def_action.h'
      include 'def_action_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'mdef_formula_info.h'
      include 'explicit.h'
      include 'ifc_input.h'

      type(orbinf), intent(inout) ::
     &     orb_info
      character, intent(in) ::
     &     env_type*(*)

      type(operator_info) ::
     &     op_info
      type(formula_info) ::
     &     form_info

c      type(operator_list), pointer ::
c     &     op_list
c      type(file_list), pointer ::
c     &     form_list
      type(action_list), pointer ::
     &     act_list, current_act
c      type(operator_array), pointer ::
c     &     op_arr(:)
c      type(file_array), pointer ::
c     &     ffform(:)
c      type(filinf), pointer ::
c     &     , ffops(:)
      type(filinf) ::
     &     ffform_opt
      integer ::
     &     ifree, nactions
      type(strinf), pointer ::
     &     str_info
      type(strmapinf), pointer ::
     &     strmap_info

      ifree = mem_setmark('do_calc')
      
      if(is_keyword_set('method.R12').gt.0)then
        explicit=.true.
      else
        explicit=.false.
      endif  

      ! set up orbital info
      call set_orbinf(orb_info,.true.)

      ifree = mem_setmark(operator_def)
      call init_operator_info(op_info)
      ! set up operators
      op_info%nops = 0
      call set_operators(op_info%op_list,op_info%nops,orb_info)

      if (op_info%nops.eq.0)
     &     call quit(0,'do_calc','no operators defined?')

      ! turn linked lists into arrays
      call update_op_arr(op_info)
      ! preliminary fix for generating unique operator IDs:
      op_info%id_cnt = op_info%nops

      ifree = mem_setmark(formula_def)
c      allocate(form_list)
c      nullify(form_list%fhand)
c      nullify(form_list%prev)
c      nullify(form_list%next)
      call init_formula_info(form_info)
      form_info%nform = 0
      ! set up (basic) formulae
      call set_formulae(form_info,op_info,orb_info)
      call update_form_arr(form_info)
      if (form_info%nform.eq.0)
     &     call quit(0,'do_calc','no formula/method defined?')

      ifree = mem_setmark(action_def)
      allocate(act_list)
      nullify(act_list%act)
      nullify(act_list%prev)
      nullify(act_list%next)
      nactions = 0
      ! set up actions
      call set_actions(act_list,nactions,
     &     form_info,op_info%op_list,op_info%nops)
      if (nactions.eq.0)
     &     call quit(0,'do_calc','no actions defined?')

      ! set up graphs
      ifree = mem_setmark(graph_def)
      allocate(str_info)
      call set_graphs_for_ops(str_info,
     &     op_info%op_list,op_info%nops,orb_info)
      ifree = mem_setmark(strmaps)
      allocate(strmap_info)
      call init_strmap(str_info,strmap_info)

      ! set up operator dimensions
      call mem_pushmark() ! push current memory section
      ifree = mem_gotomark(operator_def)
      call set_dim_for_ops(op_info%op_list,op_info%nops,
     &     str_info,orb_info)
      call mem_popmark() ! pop current memory section

      ! turn linked lists into arrays
      call update_op_arr(op_info)

      ! initialize files for operator elements
      ifree = mem_setmark(op_files)
      call init_op_files(op_info)

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
     &                        op_info,
     &                        env_type,str_info,orb_info)
          case (iaction_evaluate)
            ! evaluate a single formula expression
            call quit(1,'do_calc','action not implemented yet')
          case (iaction_setup_prc)
            call set_prc4op(current_act%act%idxopdef_out(1),
     &                      current_act%act%idxopfile_out(1,1),
     &                      current_act%act%idxopdef_in(1),
     &                      current_act%act%idxopdef_in(2),
     &                      current_act%act%idxopfile_in(1,1),
     &                      current_act%act%idxopfile_in(1,2),
     &                      op_info,
     &                      str_info,orb_info)
          case (iaction_solve_leq)
            ! Solve system of linear equations
            call quit(1,'do_calc','action not implemented yet')
          case (iaction_solve_nleq)
            ! get optimized formula file
            call file_init(ffform_opt,name_form_opt,ftyp_sq_unf,0)
            call form_opt(ffform_opt,
     &           current_act%act%nform,current_act%act%idx_formula,
     &           form_info,op_info,str_info,orb_info)
            ! Solve system of non-linear equations
            call solve_nleq(current_act%act%nop_opt,
     &                      current_act%act%nop_out,
     &                      current_act%act%idxopdef_out,
     &                      current_act%act%idxopfile_out,
     &                      current_act%act%nop_in,
     &                      current_act%act%idxopdef_in,
     &                      current_act%act%idxopfile_in,
     &                      ffform_opt,
     &                      op_info,str_info,strmap_info,orb_info
     &                     )
c            call file_delete(ffform_opt)
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

        if (.not.associated(current_act%next)) exit
        current_act => current_act%next

      end do
        
c      ! free memory allocated for operators etc.
c      deallocate(ffform)
      ! still a few deallocs missing .... !!      
      deallocate(str_info)
      ifree = mem_flushmark(op_files)
      call clean_strmap(strmap_info)
      deallocate(strmap_info)
      ifree = mem_flushmark(strmaps)
      ifree = mem_flushmark(graph_def)
      ifree = mem_flushmark(action_def)
      ifree = mem_flushmark(formula_def)
      ifree = mem_flushmark(operator_def)

      ifree = mem_flushmark('do_calc')

      return
      end
