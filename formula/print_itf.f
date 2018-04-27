*----------------------------------------------------------------------*
      subroutine print_itf(itflog,fl_head,op_info,print_form,formlog)
*----------------------------------------------------------------------*
*     Print ITF info to itflog
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_itf_tensor.h'

      integer, intent(in) ::
     &     itflog,
     &     formlog
      type(formula_item), intent(in), target ::
     &     fl_head
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in) ::
     &     print_form        ! Print to optional formulae file
      logical ::
     &     long
      integer ::
     &     inter=0    ! Counter of intermediates
      type(formula_item), pointer ::
     &     fl_item,  ! Current formula_item
     &     next_item,  ! Next formula_item
     &     prev_item  ! Previous formula_item
      integer ::
     &     nops(4,2)  ! Matrix of index info
      character(len=maxlen_bc_label) ::
     &     old_res='>',     ! Name of tensors involved in the contraction
     &     contract_next='>', ! Name of tensors involved in the contraction
     &     rename_tensor
      integer ::
     &     i,j,      ! loop indcies
     &     idx=1     ! Dummy argument to call print_form_item2()
      character(len=8) ::
     &     istr1='        ',         ! Operator 1 index
     &     istr2='        ',         ! Operator 2 index
     &     istr3='        ',         ! Result index
     &     nstr1='        ',         ! Next operator 1 index
     &     nstr2='        ',         ! Next operator 2 index
     &     nstr3='        ',         ! Next result index
     &     contract_next_index='        ',         ! Next CONTRACT result index
     &     prev_str3='>       '      ! Previous result index
      type(itf_tensor) ::
     &     itf1, itf2, itf3

      ! Point to start of linked list
      fl_item=>fl_head

      ! Loop over formula_items
      do
      select case(fl_item%command)
      case(command_end_of_formula)
!        write(itflog,*) '[END]'
      case(command_set_target_init)
!        write(itflog,*) '[INIT TARGET]',fl_item%target
      case(command_set_target_update)
        write(itflog,*) '[SET TARGET]',fl_item%target
      case(command_new_intermediate)
!        write(itflog,*) '[NEW INTERMEDIATE]',fl_item%target
!        write(itflog,'(2x,a)') trim(fl_item%interm%name)
!        write(itflog,'(2x,"attribute parentage: ",a," ",a)')
!     &                        trim(fl_item%parent1),
!     &                        trim(fl_item%parent2)
!        write(itflog,'(2x,"incore: ",i2)') fl_item%incore
!        call print_op_occ(itflog,fl_item%interm)

        inter=inter+1

!        prev_item=>fl_item%prev
!        if (prev_item%command.eq.8) then
!            write(itflog,*) "PREV res: ", prev_item%bcontr%label_res
!        end if
!        if (prev_item%command.eq.0) then
!            write(itflog,*) "PREV res: ", prev_item%target
!        end if

      case(command_del_intermediate)
        write(itflog,*) '[DELETE INTERMEDIATE]',fl_item%target
        write(itflog,'(2x,a)') trim(fl_item%label)
      case(command_add_contribution)
        write(itflog,*) '[CONTR]',fl_item%target
      case(command_add_intm)
        write(itflog,*) '[ADD]',
     &       fl_item%target
        call prt_bcontr(itflog,fl_item%bcontr)
      case(command_cp_intm)
        write(itflog,*) '[COPY]',
     &       fl_item%target
        call prt_bcontr(itflog,fl_item%bcontr)
      case(command_add_bc)
!        write(itflog,*) '[CONTRACT][ADD]',
!     &       fl_item%target
!        call prt_bcontr(itflog,fl_item%bcontr)


        ! Get index for current contraction
        call assign_index(fl_item%bcontr,istr1,istr2,istr3)
        !call itf_tensor_init(next_item%bcontr,itf1,itf2,itf3)

        !write(itflog,*) "current index ", istr3
        !write(itflog,*) "previous index ", prev_str3
        !write(itflog,*) "contract next ", contract_next
        !write(itflog,*) "contract next index ", contract_next_index

        ! Check if still part of old block (ie. old result == new result)
        ! Check result name and result index
        if (trim(old_res).ne.trim(fl_item%bcontr%label_res) .or.
     &      trim(istr3).ne.trim(prev_str3)) then
            ! Check if still part of the old block (ie. if the next
            ! result from the previous CONTRACT case == new result)
            ! This checks if a new intermediate is part of the new block
            if (trim(fl_item%bcontr%label_res).ne.trim(contract_next)
     &      .or. trim(istr3).ne.trim(contract_next_index)) then
                ! Store result using previous result index string,
                ! stored when wrote alloc
                if (old_res.ne.'>') then
                    write(itflog,*) "store ",
     &              trim(rename_tensor(old_res)), "[",
     &                             trim(contract_next_index), "]"
                end if
                write(itflog,*) 
                write(itflog,*) "alloc ",
     &          trim(rename_tensor(fl_item%bcontr%label_res)),
     &                         "[", trim(istr3), "]"
                prev_str3=istr3
                ! Update info from CONTRACT ADD
                contract_next_index=prev_str3
            end if
        end if

        call print_itf_line(fl_item%bcontr%label_res,
     &                      fl_item%bcontr%label_op1,
     &                      fl_item%bcontr%label_op2, 
     &                      fl_item%bcontr%fact,
     &                      istr1, istr2, istr3, inter, itflog)
        call clear_index(istr1,istr2, istr3)

        ! Update previous results for use next time around
        old_res=fl_item%bcontr%label_res


      case(command_add_bc_reo)
        write(itflog,*) '[CONTRACT][REORDER][ADD]',
     &       fl_item%target
        call prt_bcontr(itflog,fl_item%bcontr)
        call prt_reorder(itflog,fl_item%reo)
      case(command_bc)
!        write(itflog,*) '[CONTRACT]',
!     &       fl_item%target
!        call prt_bcontr(itflog,fl_item%bcontr)

        ! Assuming that this is called only after NEW INTERMEDIATE
        call assign_index(fl_item%bcontr,istr1,istr2,istr3)

        ! If old result does not equal the next result, then the intermediate
        ! belongs to the next block.
        ! So end current block, start new block and print intermediate
        ! line.
        next_item=>fl_item%next
        if (next_item%command.eq.8) then
            ! command_add_bc

            ! Assign index of next_item
            call assign_index(next_item%bcontr,nstr1,nstr2,nstr3)
            call itf_tensor_init(next_item%bcontr,itf1,itf2,itf3,itflog)

            ! Check if name and index are same as previous result
            if (trim(next_item%bcontr%label_res).ne.trim(old_res) .or.
     &          trim(istr3).ne.trim(prev_str3))
     &      then
                ! Store the tensor, use previous result index stored
                ! when wrote alloc
                if (old_res.ne.'>') then
                    write(itflog,*) "store ",
     &              trim(rename_tensor(old_res)), "[",
     &                             trim(prev_str3), "]"
                end if
                write(itflog,*)
!                prev_str3=istr3
            end if
            ! Update varible for use in CONTRACT ADD
            contract_next=next_item%bcontr%label_res
            contract_next_index=nstr3
        end if

        ! Assume that this is called after NEW INTERMEDIATE
        ! We need to allocate memory for next result, this should be
        ! done before the memory is allocated for the intermediate.
        ! So result tensor is allocated even if ITF block starts by
        ! constructing intermediate
        if (next_item%command.eq.8) then
            ! Do this unless the result is the same as the previous result
            if (trim(next_item%bcontr%label_res).ne.trim(old_res) .or.
     &          trim(istr3).ne.trim(prev_str3)) then
               write(itflog,*) "alloc ",
     &         trim(rename_tensor(next_item%bcontr%label_res)),
     &                        "[",trim(nstr3),"]"
               prev_str3=istr3
            end if
        end if

        call print_itf_line(fl_item%bcontr%label_res,
     &                      fl_item%bcontr%label_op1,
     &                      fl_item%bcontr%label_op2,
     &                      fl_item%bcontr%fact,
     &                      istr1, istr2, istr3, inter, itflog)
!        call construct_tensor(fl_item%bcontr%label_res,
!     &                      fl_item%bcontr%label_op1,
!     &                      fl_item%bcontr%label_op2,
!     &                      istr1, istr2, istr3,
!     &                      itf1, itf2, itf3,
!     &                      fl_item%bcontr%fact, itflog)
        call clear_index(istr1,istr2, istr3)

      case(command_bc_reo)
        write(itflog,*) '[CONTRACT][REORDER]',
     &       fl_item%target
        call prt_bcontr(itflog,fl_item%bcontr)
        call prt_reorder(itflog,fl_item%reo)
      case(command_reorder)
        write(itflog,*) '[REORDER]',fl_item%target
        call prt_reorder(itflog,fl_item%reo)
      case(command_add_reo)
        write(itflog,*) '[REORDER][ADD]',
     &       fl_item%target
        call prt_bcontr(itflog,fl_item%bcontr)
        call prt_reorder(itflog,fl_item%reo)
      case(command_symmetrise)
        write(itflog,*) '[SYMMETRISE]',fl_item%target
      case default
        write(itflog,*) 'unknown command ',fl_item%command,
     &       fl_item%target
      end select

      if (print_form) then
        call print_form_item2(formlog,'LONG',idx,fl_item,op_info)
      end if

      if (.not.associated(fl_item%next)) exit
      fl_item=>fl_item%next

      end do
      
      end
