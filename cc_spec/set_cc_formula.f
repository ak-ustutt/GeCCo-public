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
     &     name_lagrange*14 = 'cclagrange.fml',
     &     name_ccenergy*12 = 'ccenergy.fml',
     &     name_vectfunc*14 = 'ccvectfunc.fml',
     &     name_ccen*12 = 'CC energy',
     &     name_resi*14 = 'CC residual'

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
     &     idxham, idxtop, idxlag, idxr12, idxc12, idxrba, idxcba, 
     &     idxomg
      logical ::
     &     explicit

      ! advance to end of operator list:
      list_pnt => form_list
      do while (associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      ! is last entry already in use?
      if (associated(list_pnt%fhand)) then
        allocate(list_pnt%next)
        list_pnt%next%prev => list_pnt
        list_pnt => form_list%next
      end if
      allocate (list_pnt%fhand)

c      idxham = idx_oplist()
      idxham = 1
      idxtop = 2
      idxlag = 3
      idxomg = 4
      idxr12 = 5
      idxc12 = 6
      idxrba = 7
      idxcba = 8

      explicit=.false.
      if(is_keyword_set('method.R12'))then
        explicit=.true.
      endif

      ! set up Lagrangian
      nform = nform+1
      call file_init(list_pnt%fhand,name_lagrange,ftyp_sq_unf,0)
      call set_cc_lagrangian(list_pnt%fhand,nops,ops,
     &     idxham,idxlag,idxtop,idxr12,idxc12,idxrba,idxcba,explicit)

      ! set up CC-energy 
      ! (part of Lagragian that does not depend on TBAR)
      nform = nform+1
      allocate(list_pnt%next)
      list_pnt%next%prev => list_pnt
      list_pnt => list_pnt%next
      nullify(list_pnt%next)
      allocate(list_pnt%fhand)
      call file_init(list_pnt%fhand,name_ccenergy,ftyp_sq_unf,0)
      call form_indep(list_pnt%fhand,list_pnt%prev%fhand,name_ccen,
     &     idxlag,ops,nops)

      ! set up CC-residual (=vector function)
      nform = nform+1
      allocate(list_pnt%next)
      list_pnt%next%prev => list_pnt
      list_pnt => list_pnt%next
      nullify(list_pnt%next)
      allocate(list_pnt%fhand)
      call file_init(list_pnt%fhand,name_vectfunc,ftyp_sq_unf,0)
      call form_deriv(list_pnt%fhand,list_pnt%prev%prev%fhand,name_resi,
     &     idxlag,0,idxomg,ops,nops)

      return
      end
