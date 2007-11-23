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
      include 'mdef_target_info.h'
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

      type(target_info) ::
     &     tgt_info

      type(operator_info) ::
     &     op_info
      type(formula_info) ::
     &     form_info

      integer ::
     &     ifree
      type(target), pointer ::
     &     tgt
      type(strinf), pointer ::
     &     str_info
      type(strmapinf), pointer ::
     &     strmap_info

      ifree = mem_setmark('do_calc')
      
      ! set up orbital info
      call set_orbinf(orb_info,.true.)

      ! initialize target list
      call init_target_info(tgt_info)
      call set_target_list(tgt_info)

      ! initialize basis info blocks and set memory blocks:
      !  operators:
      ifree = mem_setmark(operator_def)
      call init_operator_info(op_info)

      !  formulae:
      ifree = mem_setmark(formula_def)
      call init_formula_info(form_info)
      form_info%nform = 0
c dbg
      stop 'BAUSTELLE'
c dbg

      ! graphs:
      ifree = mem_setmark(graph_def)
      ifree = mem_setmark(strmaps)

      ! initialize files for operator elements
      ifree = mem_setmark(op_files)

      ! loop until all dependencies are fulfilled
      do
        ! get next target to process
c        idx = get_next_target(tgt_info)
c        tgt => tgt_info%array(idx)%tgt

        ! which type of target ?
        select case (tgt%type)
        case(ttype_op)
          ! set up operator definitions
c          call do_opdef
        case(ttype_frm,ttype_frmopt)
          ! set up formulat definitions
c          call do_form
        case(ttype_opme)
          ! import operators/evaluate formulae
c          call do_eval
        end select

      end do
        
c      ! free memory allocated for operators etc.
c      deallocate(ffform)
      ! still a few deallocs missing .... !!      
      deallocate(str_info)
      ifree = mem_flushmark(op_files)
c      call clean_strmap(strmap_info)
c      deallocate(strmap_info)
      ifree = mem_flushmark(strmaps)
      ifree = mem_flushmark(graph_def)
      ifree = mem_flushmark(formula_def)
      ifree = mem_flushmark(operator_def)

      ifree = mem_flushmark('do_calc')

      return
      end
