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
      include 'def_me_list.h'
      include 'def_me_list_list.h'
      include 'def_me_list_array.h'
      include 'def_operator_info.h'
      include 'mdef_target_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'mdef_formula_info.h'
      include 'ifc_input.h'

      type(orbinf), intent(inout) ::
     &     orb_info
      character, intent(in) ::
     &     env_type*(*)

      type(target_info) ::
     &     tgt_info

      type(operator_info) ::
     &     op_info
      type(formula_info) ::
     &     form_info
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info

      integer ::
     &     ifree, idx, jdx, kdx, ldx
      type(target), pointer ::
     &     tgt
      type(action), pointer ::
     &     rule

      integer, external ::
     &     idx_next_target, idx_target
      

      ifree = mem_setmark('do_calc')
      
      ! set up orbital info
      call set_orbinf(orb_info,.false.)!.true.)

      ! initialize target list
      call init_target_info(tgt_info)
      call set_target_list(tgt_info,orb_info,env_type)

      ! initialize basis info blocks and set memory blocks:
      !  operators:
      ifree = mem_setmark(operator_def)
      call init_operator_info(op_info)
      ! ME-lists
      ifree = mem_setmark(me_list_def)

      !  formulae:
      ifree = mem_setmark(formula_def)
      call init_formula_info(form_info)
      form_info%nform = 0

      ! graphs:
      ifree = mem_setmark(graph_def)
      call init_str_info(str_info)
      ifree = mem_setmark(strmaps)
      call init_strmap(str_info,strmap_info)

      ! loop until all dependencies are fulfilled
      do
        ! get next target to process
        idx = idx_next_target(tgt_info)
        if (idx.le.0) exit

        tgt => tgt_info%array(idx)%tgt

        write(luout,*)
     &       'My next target: ',trim(tgt_info%array(idx)%tgt%name)

        if (tgt%n_rules.eq.0)
     &       call quit(1,'do_calc','no rules for target?')

        ! loop over rules for this target
        do jdx = 1, tgt%n_rules

          rule => tgt%rules(jdx)
          write(luout,*)
     &       'Rule: ',trim(rule%command)

          ! which type of target gets modified?
          select case (rule%type)
          case(ttype_op)
            ! set up operator definitions
            call process_operators(rule,
     &                              op_info,orb_info)
          case(ttype_frm)
            ! set up formula definitions
            call process_formulae(rule,
     &           form_info,op_info,str_info,orb_info)
          case(ttype_opme)
            ! import operators/evaluate formulae
            call process_me_lists(rule,
     &           form_info,op_info,str_info,strmap_info,orb_info)
          case default
            call quit(1,'do_calc','unknown target type')
          end select

          do kdx = 1, rule%n_update
            ldx = idx_target(rule%labels,tgt_info)
            if (ldx.le.0) cycle ! needs not necessarily be a def'd target
            call touch_target(ldx,.false.,tgt_info)
          end do

        end do

        call touch_target(idx,.true.,tgt_info)

      end do

      write(luout,*) '... all targets processed!'
        
      ! still a few deallocs missing .... !!      
      call clean_strmap(strmap_info)
      ifree = mem_flushmark(strmaps)
      ifree = mem_flushmark(graph_def)
      ifree = mem_flushmark(formula_def)
      ifree = mem_flushmark(me_list_def)
      ifree = mem_flushmark(operator_def)

      ifree = mem_flushmark('do_calc')

      return
      end
