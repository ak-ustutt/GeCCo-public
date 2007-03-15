*----------------------------------------------------------------------*
      subroutine set_cc_formula(form_list,nform,ops,nops)
*----------------------------------------------------------------------*
*     hard-wired set-up of the formulae needed for CC calculations
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'ioparam.h'
      include 'ifc_input.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'def_file_list.h'
      include 'def_operator.h'
      include 'def_operator_list.h'

      character, parameter ::
     &     name_lagrange*14 = 'cclagrange.fml'

      type(file_list), intent(inout), target ::
     &     form_list
      integer, intent(inout) ::
     &     nform
      integer, intent(in) ::
     &     nops
      type(operator) ::
     &     ops(nops)

      type(file_list), pointer ::
     &     list_pnt

      integer ::
     &     idxham, idxtop, idxlag

      ! advance to end of operator list:
      list_pnt => form_list
      do while (associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      ! is last entry already in use?
      if (associated(list_pnt%fhand)) then
        allocate(list_pnt%next)
        list_pnt%next%prev = list_pnt
        list_pnt => form_list%next
      end if
      allocate (list_pnt%fhand)

c      idxham = idx_oplist()
      idxham = 1
      idxtop = 2
      idxlag = 3

      ! set up Lagrangian
      call file_init(list_pnt%fhand,name_lagrange,ftyp_sq_unf,0)
      call set_cc_lagrangian(list_pnt%fhand,nops,ops,
     &     idxham,idxlag,idxtop)

      return
      end
