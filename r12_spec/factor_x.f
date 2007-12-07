      subroutine factor_x(form,op_info,orb_info)
*-----------------------------------------------------------------------
*     Routine which generates a formal expression for the R12 
*     X-intermediate and then searches for terms in the R12-Lagrangian
*     which may be factorised into a product of that expression and a
*     remainder. 
*     GWR December 2007.
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
     &     form_rr, form_x
      type(formula_item), pointer ::
     &     form_pnt, form_rr_pnt
      type(operator), pointer ::
     &     x_pnt
      type(operator_array), pointer ::
     &     ops_array(:)
      integer ::
     &     luinput, idxopx, idxc12, idxr12, idxrba, idxlcc
      logical ::
     &     reo
      logical, pointer ::
     &     fix_vtx(:)
      character*256 ::
     &     opin, opout

      integer, external ::
     &     idx_oplist2
      logical, external ::
     &     rd_contr


      if(ntest.ge.100)then
        write(luout,*)'=============================='
        write(luout,*)' X-intermediate factorisation '
        write(luout,*)'=============================='
      endif

      idxc12 = idx_oplist2(op_c12,op_info)
      idxrba = idx_oplist2(op_rba,op_info)
      idxr12 = idx_oplist2(op_r12,op_info)
      idxlcc = idx_oplist2(op_cclg,op_info)

      ! Form the formal X-operator.
      call add_operator(op_x_formal,op_info)
      idxopx = idx_oplist2(op_x_formal,op_info)
      x_pnt => op_info%op_arr(idxopx)%op
      allocate(ops_array(2))
      ops_array(1)%op => op_info%op_arr(idxlcc)%op
      ops_array(2)%op => op_info%op_arr(idxc12)%op
      call set_gen_intermediate(x_pnt,op_x_formal,
     &     ops_array,2,orb_info)
      deallocate(ops_array)

      ! Point to the input formula.
      form_pnt => form

      ! Initialise a temporary formula.
      call init_formula(form_rr)
      form_rr_pnt => form_rr
      call new_formula_item(form_rr_pnt,command_set_target_init,idxlcc)
      form_rr_pnt => form_rr_pnt%next

      call expand_op_product(form_rr_pnt,idxlcc,
     &     1d0,3,(/idxrba,idxc12,idxr12/),-1,-1,
     &     (/1,2,1,3,2,3/),3,.false.,op_info)

      ! Take the derivative wrt C. 
      ! NB form_x is initialised in form_deriv.
      call form_deriv3(form_x,form_rr,
     &     1,idxc12,0,idxopx,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'Formal X')
        call print_form_list(luout,form_x,op_info)
      endif

      ! Factor the V-terms out of the input formula.
      call factor_out_subexpr(form,form_x,op_info)

      ! Replace the formal V-terms with the actual integrals.
      opin = op_x_formal
      opout = op_x_inter
      call form_op_replace(opin,opout,form,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'X-factored R12 Lagrangian')
        call print_form_list(luout,form,op_info)
      endif

      call del_operator(idxopx,op_info)

      return
      end
