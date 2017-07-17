*----------------------------------------------------------------------*
!>    Driver subroutine for target setting and evaluation on target level
!>
!>    calls set_target and process_*** with the relevant rule to be evaluated.
!>    @param[inout] orb_info
!>    @param[in] env_type 
!>    @param[in] name_infile
*----------------------------------------------------------------------*
      subroutine do_calc(orb_info,env_type,name_infile)
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
      character(*), intent(in) ::
     &     name_infile

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

      type(filinf) ::
     &     fforbinf

      integer ::
     &     ifree, idx, jdx, kdx, ldx
      type(target), pointer ::
     &     tgt
      type(action), pointer ::
     &     rule

      integer, external ::
     &     idx_next_target, idx_target

      logical ::
     &     print_tgt_graph

      real(8)::
     &     cpu0_r,sys0_r,wall0_r, ! beginning of a rule
     &     cpu0_t,sys0_t,wall0_t, ! beginning of a target
     &     cpu,sys,wall ! variables for timing information
      character(len=512)::
     &     timing_msg                   ! to write a message for timings information



      ifree = mem_setmark('do_calc')
      
      ! set up orbital info
      call set_orbinf(orb_info,env_type,.false.)!.true.)

      ! print orbital informations to be accessed by an interface
      call put_orbinfo(orb_info, fforbinf)

      ! initialize target list
      call init_target_info(tgt_info)
      call set_command_prototypes(tgt_info,env_type)
      call set_target_list(tgt_info,orb_info,env_type,name_infile,
     &     fforbinf)

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

        write(lulog,*)
     &       'My next target: ',trim(tgt_info%array(idx)%tgt%name)

          call atim_csw(cpu0_t,sys0_t,wall0_t)
c        if (tgt%n_rules.eq.0)
c     &       call quit(1,'do_calc','no rules for target?')
        ! loop over rules for this target
        do jdx = 1, tgt%n_rules

          rule => tgt%rules(jdx)
          write(lulog,*)
     &       'Rule: ',trim(rule%command)
          call atim_csw(cpu0_r,sys0_r,wall0_r)
          if (.not.rule%new) then 
            ! old route:
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

          else
            ! new route
            call process_rule(rule,tgt_info,
     &           form_info,op_info,str_info,strmap_info,orb_info)
          end if
          call atim_csw(cpu,sys,wall)
          if(iprlvl.ge.10)then
          write (timing_msg,"(x,'time for rule ',A)")rule%command
          call prtim(lulog,trim(timing_msg),
     &           cpu-cpu0_r,sys-sys0_r,wall-wall0_r)
          end if
c          do kdx = 0, rule%n_update
c            ldx = idx_target(rule%labels,tgt_info)
c            if (ldx.le.0) cycle ! needs not necessarily be a def'd target
c            call touch_target(ldx,.false.,tgt_info)
c          end do

        end do

        call touch_target(idx,.true.,tgt_info)
        call atim_csw(cpu,sys,wall)
         if(iprlvl.ge.5)then
         write (timing_msg,"(x,'time for target ',A)") 
     &                   trim(tgt_info%array(idx)%tgt%name)
         call prtim(lulog,trim(timing_msg),
     &          cpu-cpu0_t,sys-sys0_t,wall-wall0_t)
         end if

      end do

      write(lulog,*) '... all targets processed!'

      call get_argument_value('general','print_tgt_graph',
     &     lval=print_tgt_graph)
      if (print_tgt_graph) call print_target_graph(tgt_info,.true.)
     
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
