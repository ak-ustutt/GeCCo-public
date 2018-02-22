*----------------------------------------------------------------------*
      subroutine print_itf(lulog,mode,idx,fl_item,op_info)
*----------------------------------------------------------------------*
*     Print ITF info to lulog
*
*     mode: "shrt", "long"; short takes only effect for [CONTR] type
*        definitions (whole diagrams, not factorized)
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     lulog
      integer, intent(inout) ::
     &     idx
      character(len=4), intent(in) ::
     &     mode
      type(formula_item), intent(in), target ::
     &     fl_item
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     long,
     &     parent_inter=.false.  ! True if intermediate is constructed from another intermediate
      integer ::
     &     inter=0    ! Counter of intermediates
      type(operator) ::
     &     tensor     ! Debug delete
      integer ::
     &     nops(4,2)  ! Matrix of index info
      character, dimension(4) ::
     &     hol=(/ 'i','j','k','l' /),
     &     par=(/ 'a','b','c','d' /),
     &     val=(/ 'u','v','w','x' /)
      character ::
     &     i_array(2,2)=reshape((/'>','>','>','>'/),(/2,2/)),
     &     p1_array(2,2)=reshape((/'>','>','>','>'/),(/2,2/)),
     &     p2_array(2,2)=reshape((/'>','>','>','>'/),(/2,2/)),
     &     k_array(2,2)=reshape((/'>','>','>','>'/),(/2,2/))
      integer ::
     &     i,j      ! loop indcies
      character(len=46) ::
     &     fomt

      long = mode.eq.'long'.or.mode.eq.'LONG'

      select case(fl_item%command)
      case(command_end_of_formula)
        write(lulog,*) '[END]'
      case(command_set_target_init)
        write(lulog,*) '[INIT TARGET]',fl_item%target
      case(command_set_target_update)
        write(lulog,*) '[SET TARGET]',fl_item%target
      case(command_new_intermediate)
        write(lulog,*) '[NEW INTERMEDIATE]',fl_item%target
        write(lulog,'(2x,a)') trim(fl_item%interm%name)
        write(lulog,'(2x,"attribute parentage: ",a," ",a)')
     &                        trim(fl_item%parent1),
     &                        trim(fl_item%parent2)
        write(lulog,'(2x,"incore: ",i2)') fl_item%incore
        call print_op_occ(lulog,fl_item%interm)

        inter=inter+1

      case(command_del_intermediate)
        write(lulog,*) '[DELETE INTERMEDIATE]',fl_item%target
        write(lulog,'(2x,a)') trim(fl_item%label)
      case(command_add_contribution)
        idx = idx+1
        if (long)
     &       write(lulog,*) '[CONTR]',fl_item%target,'( term #',idx,')'
        if (long)
     &       call prt_contr2(lulog,fl_item%contr,op_info)
        if (.not.long)
     &       call prt_contr_short(lulog,idx,fl_item%contr,op_info)
      case(command_add_intm)
        idx = idx+1
        write(lulog,*) '[ADD]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(lulog,fl_item%bcontr)
      case(command_cp_intm)
        idx = idx+1
        write(lulog,*) '[COPY]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(lulog,fl_item%bcontr)
      case(command_add_bc)
        idx = idx+1
        write(lulog,*) '[CONTRACT][ADD]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(lulog,fl_item%bcontr)
      case(command_add_bc_reo)
        idx = idx+1
        write(lulog,*) '[CONTRACT][REORDER][ADD]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(lulog,fl_item%bcontr)
        call prt_reorder(lulog,fl_item%reo)
      case(command_bc)
        idx = idx+1
        write(lulog,*) '[CONTRACT]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(lulog,fl_item%bcontr)

        ! Assuming that this is called only after NEW INTERMEDIATE
        call index_array2(lulog,fl_item%bcontr,p1_array,p2_array,
     &                    k_array) 
        
        ! Check if intermediate is constructed from previous intermediate
        if(trim(fl_item%bcontr%label_op1).eq.'_STIN0001') then
          parent_inter=.true.
        end if

        write(lulog,*) 'TENSOR:'
        if (inter<10) then
          if (parent_inter) then
            write(lulog,'(a1,i0,a1,4a1,a4,i0,a1,
     &                   4a1,a1,a1,a1,4a1,a1)'),
     &            'I',inter,'[',
     &            k_array,']+=I',inter-1,
     &            '[',p1_array,']',fl_item%bcontr%label_op2,
     &            '[',p2_array,']'
          else
            write(lulog,'(a1,i1,a1,4a1,a3,a1,a1,4a1,a1,a1,a1,4a1,a1)'),
     &            'I',inter,'[',
     &            k_array,']+=',fl_item%bcontr%label_op1,
     &            '[',p1_array,']',fl_item%bcontr%label_op2,
     &            '[',p2_array,']'
          end if
        else if (inter>=10 .and. inter<100) then
          if (parent_inter) then
            write(lulog,'(a1,i0,a1,4a1,a4,i0,a1,a1,4a1,
     &                    a1,a1,a1,4a1,a1)'),
     &            'I',inter,'[',
     &            k_array,']+=I',inter-1,'[',
     &            p1_array,']',fl_item%bcontr%label_op2,'[',p2_array,']'
          else
            write(lulog,'(a1,i0,a1,4a1,a3,a1,a1,4a1,a1,a1,a1,4a1,a1)'),
     &            'I',inter,'[',
     &            k_array,']+=',fl_item%bcontr%label_op1,'[',
     &            p1_array,']',fl_item%bcontr%label_op2,'[',p2_array,']'
          end if 
        else if (inter>=100 .and. inter<1000) then
          if (parent_inter) then
            write(lulog,'(a1,i0,a1,4a1,a4,i0,a1,a1,
     &                    4a1,a1,a1,a1,4a1,a1)'),
     &            'I',inter,'[',
     &            k_array,']+=I',inter-1,'[',
     &            p1_array,']',fl_item%bcontr%label_op2,'[',p2_array,']'
          else
            write(lulog,'(a1,i0,a1,4a1,a3,a1,a1,4a1,a1,a1,a1,4a1,a1)'),
     &            'I',inter,'[',
     &            k_array,']+=',fl_item%bcontr%label_op1,'[',
     &            p1_array,']',fl_item%bcontr%label_op2,'[',p2_array,']'
          end if
        end if

        ! Reset array
        do i=1, 2
          do j=1, 2
            i_array(i,j)='>'
            p1_array(i,j)='>'
            p2_array(i,j)='>'
            k_array(i,j)='>'
          end do
        end do

        parent_inter=.false.

      case(command_bc_reo)
        idx = idx+1
        write(lulog,*) '[CONTRACT][REORDER]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(lulog,fl_item%bcontr)
        call prt_reorder(lulog,fl_item%reo)
      case(command_reorder)
        write(lulog,*) '[REORDER]',fl_item%target
        call prt_reorder(lulog,fl_item%reo)
      case(command_add_reo)
        idx = idx+1
        write(lulog,*) '[REORDER][ADD]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(lulog,fl_item%bcontr)
        call prt_reorder(lulog,fl_item%reo)
      case(command_symmetrise)
        write(lulog,*) '[SYMMETRISE]',fl_item%target
      case default
        write(lulog,*) 'unknown command ',fl_item%command,
     &       fl_item%target
      end select
      
      end
