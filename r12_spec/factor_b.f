      subroutine factor_b(form,op_info,orb_info)
*-----------------------------------------------------------------------
*     Routine which generates a formal expression for the R12 
*     B-intermediate and then searches for terms in the R12-Lagrangian
*     which may be factorised into a product of that expression and a
*     remainder. GWR November 2007.
*-----------------------------------------------------------------------
      implicit none

      integer, parameter ::
     &     ntest = 000

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
     &     form_h, form_rfr, form_b
      type(formula_item), pointer ::
     &     form_h_pnt, form_pnt, form_rfr_pnt
      type(operator), pointer ::
     &     f_temp_pnt, b_pnt
      type(operator_array), pointer ::
     &     ops_array(:)
      integer ::
     &     luinput, idxopb, idxham, idxc12, idxr12, idxrba, idxlcc,
     &     idx_f_temp, ndef
      logical ::
     &     reo
      integer, allocatable ::
     &     occ_def(:,:,:)
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
        write(luout,*)'=============================='
        write(luout,*)' B-intermediate factorisation '
        write(luout,*)'=============================='
      endif

      idxham = idx_oplist2(op_ham,op_info)
      idxc12 = idx_oplist2(op_c12,op_info)
      idxr12 = idx_oplist2(op_r12,op_info)
      idxrba = idx_oplist2(op_rba,op_info)
      idxlcc = idx_oplist2(op_cclg,op_info)

      ! Form the formal B-operator.
      call add_operator(op_b_formal,op_info)
      idxopb = idx_oplist2(op_b_formal,op_info)
      b_pnt => op_info%op_arr(idxopb)%op
      allocate(ops_array(2))
      ops_array(1)%op => op_info%op_arr(idxlcc)%op
      ops_array(2)%op => op_info%op_arr(idxc12)%op
      call set_gen_intermediate(b_pnt,op_b_formal,
     &     ops_array,2,orb_info)
      deallocate(ops_array)

      ! Point to the input formula.
      form_pnt => form

      ! Initialise a temporary formula.
      call init_formula(form_rfr)
      form_rfr_pnt => form_rfr
      call new_formula_item(form_rfr_pnt,command_set_target_init,idxlcc)
      form_rfr_pnt => form_rfr_pnt%next

      ! Set up a temporary operator representing F.
      call add_operator(op_f_temp,op_info)
      idx_f_temp = idx_oplist2(op_f_temp,op_info)
      f_temp_pnt => op_info%op_arr(idx_f_temp)%op

      ndef = 1
      allocate(occ_def(ngastp,2,ndef))
      occ_def(1:ngastp,1,1) = (/0,0,0,1/)
      occ_def(1:ngastp,2,1) = (/0,0,0,1/)
      call set_uop(f_temp_pnt,op_f_temp,.false.,0,0,1,1,0,
     &     occ_def,ndef,orb_info)
      deallocate(occ_def)

      ! Substitute this temporary operator with the actual Hamiltonian.
      call init_formula(form_h)
      form_h_pnt => form_h
      call new_formula_item(form_h_pnt,command_set_target_init,
     &     idx_f_temp)
      form_h_pnt => form_h_pnt%next
      call expand_op_product(form_h_pnt,idx_f_temp,
     &     1d0,1,idxham,-1,-1,
     &     0,0,.false.,op_info)

      ! Form the necessary shape of the formal B-intermediate.
      call expand_op_product(form_rfr_pnt,idxlcc,
     &     1d0,4,(/idxrba,idx_f_temp,idxc12,idxr12/),-1,-1,
     &     (/1,2,1,3,2,4,3,4/),4,.false.,op_info)

      ! Replace formal F with the correct block of H.
      call expand_subexpr(form_rfr,form_h,.false.,op_info)

      ! Take the derivative wrt C. 
      ! NB form_v is initialised in form_deriv.
      call form_deriv3(form_b,form_rfr,
     &     1,idxc12,0,idxopb,op_info)

      if(ntest.ge.1000)then
        call write_title(luout,wst_title,'Formal B')
        call print_form_list(luout,form_b,op_info)
      endif

      ! Factor the B-terms out of the input formula.
      call factor_out_subexpr(form,form_b,op_info)

      ! Replace the formal B-terms with the actual integrals 
      ! (symmetrised).
      opin = op_b_formal
c      opout = op_b_inter
      opout = op_b_symm
      call form_op_replace(opin,opout,form,op_info)

      if(ntest.ge.100)then
        call write_title(luout,wst_title,'B-factored R12 Formula')
        call print_form_list(luout,form,op_info)
      endif

      ! Delete the temporary F-operator and its integral formula.
      call del_operator(idx_f_temp,op_info)
      call del_operator(idxopb,op_info)
      call dealloc_formula_list(form_h)

      return
      end
