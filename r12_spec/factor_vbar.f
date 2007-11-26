      subroutine factor_vbar(form,op_info,orb_info)
*-----------------------------------------------------------------------
*     Routine which generates a formal expression for the R12 
*     V+ intermediate and then searches for terms in the R12-Lagrangian
*     which may be factorised into a product of that expression and a
*     remainder. GWR November 2007.
*-----------------------------------------------------------------------
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_contraction_list.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'par_opnames_gen.h'
      include 'mdef_formula_info.h'
      include 'def_formula_item.h'
      include 'par_formnames_gen.h'
      include 'ifc_memman.h'

      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(formula_item), target, intent(inout) ::
     &     form

      type(contraction) ::
     &     contr
      type(formula_item), target ::
     &     form_rg, form_vbar
      type(formula_item), pointer ::
     &     form_pnt, form_rg_pnt
      type(operator), pointer ::
     &     vb_pnt
      type(operator_array), pointer ::
     &     ops_array(:)
      integer ::
     &     luinput, idxopv, idxham, idxc12, idxrba, idxlcc
      logical ::
     &     reo
      integer, pointer ::
     &     ivtx_reo(:),occ_vtx(:,:,:)
      logical, pointer ::
     &     fix_vtx(:)
      character*256 ::
     &     opin, opout

      integer, external ::
     &     idx_oplist2
      logical, external ::
     &     rd_contr


      if(ntest.ge.100)then
        write(luout,*)'==============================='
        write(luout,*)' V+-intermediate factorisation '
        write(luout,*)'==============================='
      endif

      idxham = idx_oplist2(op_ham,op_info)
      idxc12 = idx_oplist2(op_c12,op_info)
      idxrba = idx_oplist2(op_rba,op_info)
      idxlcc = idx_oplist2(op_cclg,op_info)

      ! Form the formal V+-operator.
      call add_operator(op_vbar_formal,op_info)
      idxopv = idx_oplist2(op_vbar_formal,op_info)
      vb_pnt => op_info%op_arr(idxopv)%op
      allocate(ops_array(2))
      ops_array(1)%op => op_info%op_arr(idxlcc)%op
      ops_array(2)%op => op_info%op_arr(idxc12)%op
      call set_gen_intermediate(vb_pnt,op_vbar_formal,
     &     ops_array,2,orb_info)
      deallocate(ops_array)

      ! Point to the input formula.
      form_pnt => form

      ! Initialise a temporary formula.
      call init_formula(form_rg)
      form_rg_pnt => form_rg
      call new_formula_item(form_rg_pnt,command_set_target_init,idxlcc)
      form_rg_pnt => form_rg_pnt%next

      call expand_op_product(form_rg_pnt,idxlcc,
     &     1d0,3,(/idxrba,idxc12,idxham/),-1,-1,
     &     (/1,2,1,3,2,3/),3,.false.,op_info)

      ! Take the derivative wrt C. 
      ! NB form_vbar is initialised in form_deriv.
      call form_deriv3(form_vbar,form_rg,
     &     1,idxc12,0,idxopv,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'Formal V+')
        call print_form_list(luout,form_vbar,op_info)
      endif

      ! Factor the V-terms out of the input formula.
      call factor_out_subexpr(form,form_vbar,op_info)

      ! Replace the formal V-terms with the actual integrals.
      opin = op_vbar_formal
      opout = op_vbar_inter
      call form_op_replace(opin,opout,form,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'V+-factored R12 Lagrangian')
        call print_form_list(luout,form,op_info)
      endif

      call del_operator(idxopv,op_info)

      return
      end
