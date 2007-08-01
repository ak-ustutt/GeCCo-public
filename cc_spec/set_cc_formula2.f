*----------------------------------------------------------------------*
      subroutine set_cc_formula2(form_info,op_info,orb_info)
*----------------------------------------------------------------------*
*     hard-wired set-up of the formulae needed for CC calculations
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'cc_routes.h'
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

c      type(formula_list), pointer ::
c     &     list_pnt
      type(formula), pointer ::
     &     cclg_pnt, form_pnt, form_eta

      integer ::
     &     idxham, idxtop, idxtba, idxomg, idxecc, idxhhat, idx,
     &     idxtbtrf, idxlcc, idxeta

      ! explicit interface does not work with ifort
      integer, external ::
     &     idx_oplist2, idx_formlist

        call write_title(luout,wst_section,
     &     'Setting up formulae for coupled-cluster calculations')      

      idxham = idx_oplist2(op_ham,op_info)
      if (idxham.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_ham))
      idxecc = idx_oplist2(op_ccen,op_info)
      if (idxecc.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_ccen))
      idxlcc = idx_oplist2(op_cclg,op_info)
      if (idxlcc.le.0)
     &     call quit(1,'set_cc_formula','operator not on list: '
     &     //trim(op_cclg))
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

c      call test_formgen3(op_info,orb_info)

      ! set up Lagrangian
      ! new entry
      call add_formula(form_info,label_cclg0)

      ! get index of formula and point there (needed below)
      idx = idx_formlist(label_cclg0,form_info)
      cclg_pnt => form_info%form_arr(idx)%form

      if (iprlvl.gt.0)
     &     write(luout,'(2x,"* ",a)')
     &     'Setting ground state Lagrangian'
      call set_cc_lagrangian2(cclg_pnt,op_info,
     &     idxham,idxtba,idxtop,idxlcc)

      ! is Hhat operator on list?
      idxhhat = idx_oplist2(op_hhat,op_info)
      if (idxhhat.gt.0) then
        call add_formula(form_info,label_cchhat)
        idx = idx_formlist(label_cchhat,form_info)
        form_pnt => form_info%form_arr(idx)%form
        if (iprlvl.gt.0)
     &     write(luout,'(2x,"* ",a)')
     &     'Setting e^-T1 H e^T1'
        call set_hhat2(form_pnt,op_info,
     &       idxhhat,idxham,idxtop)
      end if
c
      ! set up CC-energy 
      ! (part of Lagragian that does not depend on TBAR)
      call add_formula(form_info,label_ccen0)
      idx = idx_formlist(label_ccen0,form_info)
      form_pnt => form_info%form_arr(idx)%form

      if (iprlvl.gt.0)
     &     write(luout,'(2x,"* ",a)')
     &     'Setting ground state energy and residual'
      call form_indep2(form_pnt,
     &     cclg_pnt,
     &     label_ccen0,title_ccen0,idxecc,
     &     1,idxtba,
     &     op_info)

      ! set up CC-residual (=vector function)
      call add_formula(form_info,label_ccrs0)
      idx = idx_formlist(label_ccrs0,form_info)
      form_pnt => form_info%form_arr(idx)%form

      call form_deriv2(form_pnt,cclg_pnt,
     &     label_ccrs0,title_ccrs0,
     &     1,idxtba,0,idxomg,
     &     op_info)

      if (solve_tbar) then
        if (iprlvl.gt.0)
     &     write(luout,'(2x,"* ",a)')
     &     'Setting ground state left-hand equations'
        call add_formula(form_info,label_cctbar_a)
        idx = idx_formlist(label_cctbar_a,form_info)
        idxtbtrf = idx_oplist2(op_tbar_a,op_info)
        form_pnt => form_info%form_arr(idx)%form
        ! d L / d t = tbar A + eta
        call form_deriv2(form_pnt,cclg_pnt,
     &       label_cctbar_a,title_cctbar_a,
     &       1,idxtop,0,idxtbtrf,
     &       op_info)

        ! split into tbar A and eta contribution
        call add_formula(form_info,label_cceta)
        idx = idx_formlist(label_cceta,form_info)
        form_eta => form_info%form_arr(idx)%form
        idxeta = idx_oplist2(op_eta,op_info)

        call leq_post(form_pnt,form_eta,form_pnt,
     &                idxtbtrf,idxeta,idxtba,
     &                title_cctbar_a,title_cceta,
     &                op_info)
      end if

      return
      end
