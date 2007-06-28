*----------------------------------------------------------------------*
      subroutine set_cc_formula2(form_info,op_info,orb_info)
*----------------------------------------------------------------------*
*     hard-wired set-up of the formulae needed for CC calculations
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'par_formnames_gen.h'
      include 'par_opnames_gen.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
c      include 'def_filinf.h'
c      include 'def_file_list.h'
c      include 'def_operator.h'
c      include 'def_operator_list.h'
      include 'mdef_formula_info.h'
      include 'ifc_input.h'
c      include 'ifc_operators.h'

      type(formula_info), intent(inout) ::
     &     form_info
      type(operator_info) ::
     &     op_info
      type(orbinf) ::
     &     orb_info

      type(formula_list), pointer ::
     &     list_pnt
      type(formula), pointer ::
     &     cclg_pnt

      integer ::
     &     idxham, idxtop, idxtba, idxomg, idxecc, idxhhat

      ! explicit interface does not work with ifort
      integer, external ::
     &     idx_oplist2

      ! advance to end of formula list:
      list_pnt => form_info%form_list
      do while (associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      ! is last entry already in use?
      if (associated(list_pnt%form)) then
        allocate(list_pnt%next)
        list_pnt%next%prev => list_pnt
        list_pnt => list_pnt%next
      end if
      allocate (list_pnt%form)

      idxham = idx_oplist2(op_ham,op_info)
      if (idxham.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_ham))
      idxecc = idx_oplist2(op_ccen,op_info)
      if (idxecc.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_ccen))
      idxtop = idx_oplist2(op_top,op_info)
      if (idxtop.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_top))
      idxtba = idx_oplist2(op_tbar,op_info)
      if (idxtba.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_tbar))
      idxomg = idx_oplist2(op_omg,op_info)
      if (idxomg.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_omg))

      ! set up Lagrangian
      form_info%nform = form_info%nform+1
      cclg_pnt => list_pnt%form
c dbg
c      call test_formgen2(op_info,orb_info)
c dbg

      call set_cc_lagrangian2(list_pnt%form,op_info,
     &     idxham,idxtba,idxtop,idxecc)

      ! is Hhat operator on list?
      idxhhat = idx_oplist2(op_hhat,op_info)
      if (idxhhat.gt.0) then
        form_info%nform = form_info%nform+1
        allocate(list_pnt%next)
        list_pnt%next%prev => list_pnt
        list_pnt => list_pnt%next
        nullify(list_pnt%next)
        allocate(list_pnt%form)
        call set_hhat2(list_pnt%form,op_info,
     &       idxhhat,idxham,idxtop)
      end if
c
      ! set up CC-energy 
      ! (part of Lagragian that does not depend on TBAR)
      form_info%nform = form_info%nform+1
      allocate(list_pnt%next)
      list_pnt%next%prev => list_pnt
      list_pnt => list_pnt%next
      nullify(list_pnt%next)
      allocate(list_pnt%form)
      call form_indep2(list_pnt%form,
     &     cclg_pnt,
     &     label_ccen0,title_ccen0,
     &     1,idxtba,
     &     op_info)

      ! set up CC-residual (=vector function)
      form_info%nform = form_info%nform+1
      allocate(list_pnt%next)
      list_pnt%next%prev => list_pnt
      list_pnt => list_pnt%next
      nullify(list_pnt%next)
      allocate(list_pnt%form)
      call form_deriv2(list_pnt%form,cclg_pnt,
     &     label_ccrs0,title_ccrs0,
     &     1,idxtba,0,idxomg,
     &     op_info)

      return
      end
