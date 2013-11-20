*----------------------------------------------------------------------*
      subroutine diag_type2cnt_info(cnt_info,idxop,contr,op_info)
*----------------------------------------------------------------------*
*     put information about the type of diagonal operator structure
*     of contracted operators (but not for intermediates)
*     and of the result operator (if final contr.) on cnt_info
*
*     matthias, Nov. 2010
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_contraction_info.h'

      type(contraction_info), intent(inout) ::
     &     cnt_info
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     idxop(2)

      integer ::
     &     idxmel, idxres, idx

      ! final level? set diag_type of result
      if (contr%narc.eq.0) then
        idxres = contr%idx_res
        idxmel = op_info%op2list(idxres)
        if (idxmel.le.0) then
          call prt_contr2(lulog,contr,op_info)
          write(lulog,*) op_info%op2list(1:op_info%nops)
          call quit(1,'diag_type2cnt_info',
     &       'No list associated with operator '//
     &       trim(op_info%op_arr(idxres)%op%name))
        end if
        cnt_info%diag_type12 = op_info%mel_arr(idxmel)%mel%diag_type
      end if

      ! set diag_type of operators 1 and 2
      ! (not yet available for intermediates)
      do idx = 1, 2
        if (idxop(idx).gt.0) then
          idxmel = op_info%op2list(idxop(idx))
          if (idxmel.le.0) then
            call prt_contr2(lulog,contr,op_info)
            write(lulog,*) op_info%op2list(1:op_info%nops)
            call quit(1,'diag_type2cnt_info',
     &         'No list associated with operator '//
     &         trim(op_info%op_arr(idxop(idx))%op%name))
          end if
          if (idx.eq.1) then
            cnt_info%diag_type1 = op_info%mel_arr(idxmel)%mel%diag_type
          else
            cnt_info%diag_type2 = op_info%mel_arr(idxmel)%mel%diag_type
          end if
        end if
      end do

      return
      end
