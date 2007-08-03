*----------------------------------------------------------------------*
      subroutine set_operators(op_info,orb_info)
*----------------------------------------------------------------------*
*     driver routine for setting up operator info structures
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'ifc_input.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'explicit.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf) ::
     &     orb_info

      call set_h_operators(op_info,orb_info,explicit)
      if (is_keyword_set('method.CC').gt.0) then
        call set_cc_operators(op_info,orb_info)
      endif
      if(explicit)then
        call set_r12_operators(op_info,orb_info)
      end if

      ! further operator definitions may follow

      if (iprlvl.gt.0) then
        write(luout,*) 'List of defined operators:'
        call list_operators(luout,op_info)
      end if
      
      return

      end
