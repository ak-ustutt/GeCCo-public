*----------------------------------------------------------------------*
      subroutine set_cc_formula(form_list,nform,ops,nops)
*----------------------------------------------------------------------*
*     hard-wired set-up of the formulae needed for CC calculations
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'ioparam.h'
      include 'ifc_input.h'
c      include 'ifc_operators.h'
      include 'par_opnames_gen.h'
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
     &     idxham, idxtop, idxlag, idxomg

      ! explicit interface does not work with ifort
      integer, external ::
     &     idx_oplist

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

      idxham = idx_oplist(op_ham,ops,nops)
      if (idxham.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_ham))
      idxtop = idx_oplist(op_top,ops,nops)
      if (idxtop.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_top))
      idxlag = idx_oplist(op_tbar,ops,nops)
      if (idxlag.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_tbar))
      idxomg = idx_oplist(op_omg,ops,nops)
      if (idxomg.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_omg))

      ! set up Lagrangian
      nform = nform+1
      call file_init(list_pnt%fhand,name_lagrange,ftyp_sq_unf,0)
      call set_cc_lagrangian(list_pnt%fhand,nops,ops,
     &     idxham,idxlag,idxtop)

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
     &     1,idxlag,
     &     ops,nops)

      ! set up CC-residual (=vector function)
      nform = nform+1
      allocate(list_pnt%next)
      list_pnt%next%prev => list_pnt
      list_pnt => list_pnt%next
      nullify(list_pnt%next)
      allocate(list_pnt%fhand)
      call file_init(list_pnt%fhand,name_vectfunc,ftyp_sq_unf,0)
      call form_deriv(list_pnt%fhand,list_pnt%prev%prev%fhand,name_resi,
     &     1,idxlag,0,idxomg,
     &     ops,nops)

      return
      end
