*----------------------------------------------------------------------*
      subroutine op_adv_state(label_op,nop,max_state,op_info,use_1,
     &     last_state,new_state)
*----------------------------------------------------------------------*
*     advance the state of the operators in label_op
*     if the state goes beyond max_state, restore to the first state
*     case that last_state is set to T. If new_state is present,
*     set the state to that particular state.
*
*     yuri 2014
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'mdef_target_info.h'

      integer, parameter ::
     &     ntest = 10

      integer, intent(in) ::
     &     max_state, nop
      character(*), intent(in) ::
     &     label_op(nop)
      type(operator_info), intent(inout) ::
     &     op_info
      logical, intent(in) ::
     &     use_1
      logical, intent(out), optional ::
     &     last_state
      integer, intent(in), optional ::
     &     new_state

      integer ::
     &     idx, i_state, i_label, size_name, ios_var

      character(len=3) ::
     &     c_st

      character(len_target_name), external ::
     &     state_label
      integer, external ::
     &     idx_oplist2

      character(mxlen_melabel) ::
     &     current_mel, new_mel


      do i_label = 1,nop
       
       if(ntest.GE.10) write(lulog,*) 'Advancing state in operator ',
     &      trim(label_op(i_label))
       
       idx = idx_oplist2(label_op(i_label),op_info)
       if (idx.le.0)
     &      call quit(1,'op_adv_state',
     &      'operator not found: "'//trim(label_op(i_label))//'"')

       current_mel = trim(op_info%op_arr(idx)%op%assoc_list)

       if(ntest.GE.100) write(lulog,*)
     &      "current ME: ", trim(current_mel)

       size_name = index(trim(current_mel),"_",.true.)

       if(size_name.GT.0) then
        read(current_mel(size_name+1:),fmt='(i2)',iostat=ios_var)
     &       i_state
        if(ios_var.NE.0) then
         size_name = len_trim(current_mel)+1
         i_state = 1
        end if
       else
        size_name = len_trim(current_mel)
        i_state = 1
       end if

       if (present(new_state)) then
        i_state = new_state
        if(ntest.GT.1) write(lulog,FMT='(" State set to: ",i0)') i_state
       else
        if (i_state.ge.max_state) then
         i_state = 1
         if(ntest.GE.10) write(lulog,FMT='(" All states processed.")')
         if(present(last_state)) last_state = .true.
        else
         i_state = i_state + 1
         if(present(last_state)) last_state = .false.
         if(ntest.GE.10) write(lulog,FMT='(" Next state: ",i0)') i_state
        end if
       end if
       c_st = state_label(i_state,use_1)
       new_mel = current_mel(1:size_name-1)//trim(c_st)
       if(ntest.GE.100) write(lulog,*)
     &      "new ME: ", trim(new_mel)

       call assign_me2op(trim(new_mel),
     &      trim(label_op(i_label)),op_info)

      end do

      return
      end
