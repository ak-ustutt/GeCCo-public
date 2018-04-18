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
     &     long
      integer ::
     &     inter=0    ! Counter of intermediates
      type(formula_item), pointer ::
     &     next_item,  ! Next formula_item
     &     prev_item  ! Next formula_item
      integer ::
     &     nops(4,2)  ! Matrix of index info
      character(len=maxlen_bc_label) ::
     &     old_res='>',     ! Name of tensors involved in the contraction
     &     contract_next='>'     ! Name of tensors involved in the contraction
      integer ::
     &     i,j      ! loop indcies
      character(len=8) ::
     &     istr1='        ',         ! Operator 1 index
     &     istr2='        ',         ! Operator 2 index
     &     istr3='        ',         ! Result index
     &     nstr1='        ',         ! Next operator 1 index
     &     nstr2='        ',         ! Next operator 2 index
     &     nstr3='        ',         ! Next result index
     &     contract_next_index='        ',         ! Next CONTRACT result index
     &     prev_str3='>       '      ! Previous result index

      long = mode.eq.'long'.or.mode.eq.'LONG'

      select case(fl_item%command)
      case(command_end_of_formula)
!        write(lulog,*) '[END]'
      case(command_set_target_init)
        write(lulog,*) '[INIT TARGET]',fl_item%target
      case(command_set_target_update)
        write(lulog,*) '[SET TARGET]',fl_item%target
      case(command_new_intermediate)
!        write(lulog,*) '[NEW INTERMEDIATE]',fl_item%target
!        write(lulog,'(2x,a)') trim(fl_item%interm%name)
!        write(lulog,'(2x,"attribute parentage: ",a," ",a)')
!     &                        trim(fl_item%parent1),
!     &                        trim(fl_item%parent2)
!        write(lulog,'(2x,"incore: ",i2)') fl_item%incore
!        call print_op_occ(lulog,fl_item%interm)

        inter=inter+1

!        prev_item=>fl_item%prev
!        if (prev_item%command.eq.8) then
!            write(lulog,*) "PREV res: ", prev_item%bcontr%label_res
!        end if
!        if (prev_item%command.eq.0) then
!            write(lulog,*) "PREV res: ", prev_item%target
!        end if

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
!        write(lulog,*) '[CONTRACT][ADD]',
!     &       fl_item%target,'( term #',idx,')'
!        call prt_bcontr(lulog,fl_item%bcontr)


        ! Get index for current contraction
        call assign_index(fl_item%bcontr,istr1,istr2,istr3)

        !write(lulog,*) "current index ", istr3
        !write(lulog,*) "previous index ", prev_str3
        !write(lulog,*) "contract next ", contract_next
        !write(lulog,*) "contract next index ", contract_next_index

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
                    write(lulog,*) "store ", trim(old_res), "[",
     &                             trim(contract_next_index), "]"
                end if
                write(lulog,*) 
                write(lulog,*) "alloc ", trim(fl_item%bcontr%label_res),
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
     &                      istr1, istr2, istr3, inter, lulog)
        call clear_index(istr1,istr2, istr3)

        ! Update previous results for use next time around
        old_res=fl_item%bcontr%label_res


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
        call assign_index(fl_item%bcontr,istr1,istr2,istr3)

        ! Need to compare previous result from before the intermediate
        ! and next result after the intermediate. If these aren't the
        ! same, then store previous result and allocate for the next
        ! result and the intermediate
        next_item=>fl_item%next

        ! If old result does not equal the next result, then the intermediate
        ! belongs to the next block.
        ! So end current block, start new block and print intermediate
        ! line.
        next_item=>fl_item%next
        if (next_item%command.eq.8) then
            ! command_add_bc

            ! Assign index of next_item
            call assign_index(next_item%bcontr,nstr1,nstr2,nstr3)

            ! Check if name and index are same as previous result
            if (trim(next_item%bcontr%label_res).ne.trim(old_res) .or.
     &          trim(istr3).ne.trim(prev_str3))
     &      then
                ! Store the tensor, use previous result index stored
                ! when wrote alloc
                if (old_res.ne.'>') then
                    write(lulog,*) "store ", trim(old_res), "[",
     &                             trim(prev_str3), "]"
                end if
                write(lulog,*)
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
               write(lulog,*) "alloc ",trim(next_item%bcontr%label_res),
     &                        "[",trim(nstr3),"]"
               prev_str3=istr3
            end if
        end if

        call print_itf_line(fl_item%bcontr%label_res,
     &                      fl_item%bcontr%label_op1,
     &                      fl_item%bcontr%label_op2,
     &                      fl_item%bcontr%fact,
     &                      istr1, istr2, istr3, inter, lulog)
        call clear_index(istr1,istr2, istr3)

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
