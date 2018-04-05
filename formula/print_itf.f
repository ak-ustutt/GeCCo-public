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
      character(len=4)
     &     index1
      character(len=maxlen_bc_label) ::
     &     tensor1, tensor2
      integer ::
     &     i,j      ! loop indcies
      character(len=5) ::
     &     s_int,
     &     p_int
      character(len=4) ::
     &     istr1='    ', istr2='    ', istr3='    '  ! Temporary index string, can contain whitspace inbetween index
      character(len=50) ::
     &     itf_line
      integer ::
     &     last=1, actual=1

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
        call index_array(lulog,fl_item%bcontr,p1_array,p2_array,
     &                    k_array)

        ! Get index letters from arrays and copy to strings
        do i=1, 2
          do j=1, 2
            if (p1_array(j,i).ne.'>') then
              istr1(j+(2*i-2):)=p1_array(j,i)
            end if
            if (p2_array(j,i).ne.'>') then
              istr2(j+(2*i-2):)=p2_array(j,i)
            end if
            if (k_array(j,i).ne.'>') then
              istr3(j+(2*i-2):)=k_array(j,i)
            end if
          end do
        end do

        ! Remove whitespaces from wihtin strings
        call remove_whitespace(istr1)
        call remove_whitespace(istr2)
        call remove_whitespace(istr3)

        tensor1=fl_item%bcontr%label_op1
        tensor2=fl_item%bcontr%label_op2

        ! Change integral tensor name
        if (tensor1.eq.'INT_D') then
          tensor1='K    '
        else if (tensor2.eq.'INT_D') then
          tensor2='K'
        end if
          
        write(s_int,'(i0)') inter

        ! Output ITF line
        write(lulog,*) 'TENSOR:'
        if (tensor1.eq.'_STIN0001') then
          ! Previous intermediate appears on rhs
          tensor1='I'
          write(p_int,'(i0)') inter-1
          itf_line='I'//trim(s_int)//'['//trim(istr3)//']+='//
     &        trim(tensor1)//trim(p_int)//'['//trim(istr1)//']'//
     &        trim(tensor2)//'['//trim(istr2)//']'
          write(lulog,*) trim(itf_line)
        else
          itf_line='I'//trim(s_int)//'['//trim(istr3)//']+='//
     &        trim(tensor1)//'['//trim(istr1)//']'//trim(tensor2)//'['//
     &        trim(istr2)//']'
          write(lulog,*) trim(itf_line)
        end if

        ! Reset array and strings
        do i=1, 2
          do j=1, 2
            i_array(i,j)='>'
            p1_array(i,j)='>'
            p2_array(i,j)='>'
            k_array(i,j)='>'
          end do
        end do

        istr1='    '
        istr2='    '
        istr3='    '

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
