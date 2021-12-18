*----------------------------------------------------------------------*
      integer function idx_xret(label,op_info,depend)
*----------------------------------------------------------------------*
*     run set_formula_dependencies() first
*     then, the function tells you on which index on xret 
*     the subroutine frm_schedX() will return the info value
*     for the ME-list labelled with "label" 
*     may also serve as a test, whether the formula analyzed in
*     set_formula_dependencies() does contain an update for the
*     required ME-list
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_operator_info.h'
      include 'def_dependency_info.h'

      character(*), intent(in) ::
     &     label
      type(operator_info), intent(in) ::
     &     op_info
      type(dependency_info), intent(in) ::
     &     depend

      integer ::
     &     idx, jdx, idxmel

      integer, external ::
     &     idx_mel_list

      idxmel = idx_mel_list(label,op_info)
      jdx = -1
      do idx = 1, depend%ntargets
        if (depend%idxlist(idx).eq.idxmel) then
          jdx = idx
          exit
        end if
      end do

      idx_xret = jdx

      return
      end
