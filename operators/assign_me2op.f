*----------------------------------------------------------------------*
      subroutine assign_me2op(label_mel,label_op,op_info)
*----------------------------------------------------------------------*
*     assign an existing ME-list with label "label_mel"
*     to operator "label_op" (must exist as well),
*     but NOT the operator to this list!
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'

      character(*), intent(in) ::
     &     label_mel, label_op
      type(operator_info), intent(inout) ::
     &     op_info

      integer ::
     &     idx
      type(me_list), pointer ::
     &     mel

      integer, external ::
     &     idx_mel_list, idx_oplist2

      ! the label must exist previously ...
      idx = idx_mel_list(label_mel,op_info)
      if (idx.lt.0)
     &     call quit(1,'assign_me2op',
     &     'list does not already exist: "'//trim(label_mel)//'"')
      mel => op_info%mel_arr(idx)%mel

      ! ... as well as the operator-label
      idx = idx_oplist2(label_op,op_info)
      if (idx.le.0)
     &       call quit(1,'assign_me2op',
     &     'operator not found: "'//trim(label_op)//'"')
      ! set associated list on operator:
      op_info%op_arr(idx)%op%assoc_list(1:mxlen_melabel) = ' '
      op_info%op_arr(idx)%op%assoc_list(1:mxlen_melabel) = 
     &       trim(label_mel)

      ! update op_list array in order to set up the lookup-table
      call update_op_arr(op_info)

      if (iprlvl.ge.20) then
        write(luout,'(3x,4a)')
     &       'assigned list: ',trim(mel%label),' (file ',
     &       ') to operator: ',trim(op_info%op_arr(idx)%op%name)
      end if

      return
      end
