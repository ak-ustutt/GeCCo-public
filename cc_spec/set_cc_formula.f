*----------------------------------------------------------------------*
      subroutine set_cc_formula(form_info,ops,nops)
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
      include 'mdef_formula_info.h'
      include 'par_formnames_gen.h'
      include 'explicit.h'

      type(formula_info), intent(inout) ::
     &     form_info
      integer, intent(in) ::
     &     nops
      type(operator) ::
     &     ops(nops)

      type(formula_list), pointer ::
     &     list_pnt

      integer ::
     &     idxham, idxtop, idxtba, idxr12, idxc12, idxrba, idxcba, 
     &     idxomg, idxecc

      ! explicit interface does not work with ifort
      integer, external ::
     &     idx_oplist

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

      idxham = idx_oplist(op_ham,ops,nops)
      if (idxham.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_ham))
      idxecc = idx_oplist(op_ccen,ops,nops)
      if (idxecc.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_ccen))
      idxtop = idx_oplist(op_top,ops,nops)
      if (idxtop.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_top))
      idxtba = idx_oplist(op_tbar,ops,nops)
      if (idxtba.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_tbar))
      idxomg = idx_oplist(op_omg,ops,nops)
      if (idxomg.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_omg))

      if(explicit)then
        idxr12=idx_oplist(op_r12,ops,nops)
        if(idxr12.le.0)
     &       call quit(1,'set_cc_formula','operator not on list: '
     &       //trim(op_r12))
        idxrba=idx_oplist(op_rba,ops,nops)
        if(idxrba.le.0)
     &       call quit(1,'set_cc_formula','operator not on list: '
     &       //trim(op_rba))
        idxc12=idx_oplist(op_c12,ops,nops)
        if(idxc12.le.0)
     &       call quit(1,'set_cc_formula','operator not on list: '
     &       //trim(op_c12))
        idxcba=idx_oplist(op_cba,ops,nops)
        if(idxcba.le.0)
     &       call quit(1,'set_cc_formula','operator not on list: '
     &       //trim(op_cba))
      endif
  
      ! set up Lagrangian
      form_info%nform = form_info%nform+1
      call set_cc_lagrangian(list_pnt%form,nops,ops,
     &     idxecc,idxham,idxtba,idxtop,idxr12,idxc12,idxrba,idxcba)
      ! is Hhat operator on list?
c      idxhhat = idx_oplist(op_hhat,ops,nops)
c      if (idxhhat.gt.0) then
c        form_info%nform = form_info%nform+1
c        allocate(list_pnt%next)
c        list_pnt%next%prev => list_pnt
c        list_pnt => list_pnt%next
c        nullify(list_pnt%next)
c        allocate(list_pnt%form)
c        call set_hhat(list_pnt%form,nops,ops,
c     &       idxhhat,idxham,idxtop)
c      end if

      ! set up CC-energy 
      ! (part of Lagragian that does not depend on TBAR)
      form_info%nform = form_info%nform+1
      allocate(list_pnt%next)
      list_pnt%next%prev => list_pnt
      list_pnt => list_pnt%next
      nullify(list_pnt%next)
      allocate(list_pnt%form)
      call form_indep(list_pnt%form,
     &     list_pnt%prev%form,
     &     label_ccen0,title_ccen0,
     &     1,idxtba,
     &     ops,nops)

      ! set up CC-residual (=vector function)
      form_info%nform = form_info%nform+1
      allocate(list_pnt%next)
      list_pnt%next%prev => list_pnt
      list_pnt => list_pnt%next
      nullify(list_pnt%next)
      allocate(list_pnt%form)
c      call file_init(list_pnt%form%fhand,name_vectfunc,ftyp_sq_unf,0)
      call form_deriv(list_pnt%form,list_pnt%prev%prev%form,
     &     label_ccrs0,title_ccrs0,
     &     1,idxtba,0,idxomg,
     &     ops,nops)

      return
      end
