*----------------------------------------------------------------------*
      subroutine assign_op2me(label_mel,label_op,op_info)
*----------------------------------------------------------------------*
*     assign an operator "label_op" (must exist)
*     to an existing ME-list with label "label_mel",
*     but NOT the list to the operator!
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
     &     call quit(1,'assign_op2me',
     &     'list does not already exist: "'//trim(label_mel)//'"')
      mel => op_info%mel_arr(idx)%mel

      ! ... as well as the operator-label
      idx = idx_oplist2(label_op,op_info)
      if (idx.le.0)
     &       call quit(1,'assign_op2me',
     &     'operator not found: "'//trim(label_op)//'"')
      mel%op => op_info%op_arr(idx)%op

      ! update op_list array in order to set up the lookup-table
      ! is this needed here???
      call update_op_arr(op_info)

      if (iprlvl.ge.20) then
        write(lulog,'(3x,4a)')
     &       'assigned operator: ',trim(mel%op%name),' to me-list:  ',
     &       trim(mel%label)
      end if

      return
      end
