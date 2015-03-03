      subroutine test_formgen2(op_info,orb_info)

      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'ifc_operators.h'
      include 'par_opnames_gen.h'

      type(operator_info) ::
     &     op_info
      type(orbinf) ::
     &     orb_info

      type(operator), pointer ::
     &     op

      integer ::
     &     idxopa, idxopb, idxopc, idxopd, idxopr, idxoph, iformal
      
      type(formula_item), target ::
     &     form_a_dc, form_bhb

      type(formula_item), pointer ::
     &     fl_a_dc_pnt, fl_bhb_pnt

      type(operator_array), pointer ::
     &     defop(:)

      integer, external ::
     &     idx_oplist2
      call quit(1,'test_formgen2','call to obsolete routine')


      write(lulog,*) '-------------------'
      write(lulog,*) ' adding A operator'
      write(lulog,*) '-------------------'
      call add_operator('A',op_info)

      iformal=3
      idxopa = op_info%nops
      op => op_info%op_arr(idxopa)%op
      call set_xop(op,'A',.false.,0,0,1,0,0,
     &     2,2,0,0,iformal,orb_info)

      write(lulog,*) '-------------------'
      write(lulog,*) ' adding B operator'
      write(lulog,*) '-------------------'
      call add_operator('B',op_info)

      idxopb = op_info%nops
      op => op_info%op_arr(idxopb)%op
      call set_xop(op,'B',.false.,0,0,1,0,0,
     &     2,2,0,iformal,orb_info)

      write(lulog,*) '-------------------'
      write(lulog,*) ' adding C operator'
      write(lulog,*) '-------------------'
      call add_operator('C',op_info)

      idxopc = op_info%nops
      op => op_info%op_arr(idxopc)%op
      call set_xop(op,'C',.false.,0,0,1,0,0,
     &     2,2,0,iformal,orb_info)

      write(lulog,*) '-------------------'
      write(lulog,*) ' adding R operator'
      write(lulog,*) '-------------------'
      call add_operator('R',op_info)

      idxopr = op_info%nops
      op => op_info%op_arr(idxopr)%op
      call set_xop(op,'R',.false.,0,0,1,0,0,
     &     2,2,0,iformal,orb_info)

      write(lulog,*) '-------------------'
      write(lulog,*) ' adding D operator'
      write(lulog,*) '-------------------'
      call add_operator('D',op_info)

      idxopd = op_info%nops
      op => op_info%op_arr(idxopd)%op
      op%name = 'D'
      op%n_occ_cls = 1
      allocate(op%ihpvca_occ(ngastp,2,1))
      allocate(op%ica_occ(2,1))
      op%ihpvca_occ(1:ngastp,1:2,1) = 0
      op%ihpvca_occ(1,1,1) = 2
      op%ihpvca_occ(1,2,1) = 2

      write(lulog,*) '-------------------'
      write(lulog,*) ' generate B = A+DC'
      write(lulog,*) '-------------------'
      ! define: B = A + DC
      call init_formula(form_a_dc)
      fl_a_dc_pnt => form_a_dc
      call new_formula_item(fl_a_dc_pnt,command_set_target_init,idxopb)
      fl_a_dc_pnt => fl_a_dc_pnt%next
      call expand_op_product(fl_a_dc_pnt,idxopb,
     &     1d0,1,idxopa,-1,-1,
     &     0,0,.false.,op_info)
      do while(associated(fl_a_dc_pnt%next))
        fl_a_dc_pnt => fl_a_dc_pnt%next
      end do
      call expand_op_product(fl_a_dc_pnt,idxopb,
     &     1d0,2,(/idxopd,idxopc/),-1,-1,
     &     (/1,2/),1,.false.,op_info)

      call print_form_list(lulog,form_a_dc,op_info)

      write(lulog,*) '-----------------------------'
      write(lulog,*) ' generate R = e^{-B} H e^{B}'
      write(lulog,*) '-----------------------------'
      ! generate: e^{-B} H e^{B}
      idxoph = idx_oplist2(op_ham,op_info)
      call init_formula(form_bhb)
      fl_bhb_pnt => form_bhb
      call new_formula_item(fl_bhb_pnt,command_set_target_init,idxopr)
      fl_bhb_pnt => fl_bhb_pnt%next
      
      call expand_op_bch(fl_bhb_pnt,4,idxopr,
     &     1d0,-1,idxoph,1d0,idxopb,-1,-1,op_info)

      call print_form_list(lulog,form_bhb,op_info)

      ! replace B -> A + DC
      write(lulog,*) '------------------------------------'
      write(lulog,*) ' insert B: R = e^{-A-DC} H e^{B+DC}'
      write(lulog,*) '------------------------------------'

      call expand_subexpr(form_bhb,form_a_dc,0,op_info)

      call print_form_list(lulog,form_bhb,op_info)

      write(lulog,*) '------------------------------------'
      write(lulog,*) ' after summing terms: '
      write(lulog,*) '------------------------------------'

      call sum_terms(form_bhb,op_info)

      call print_form_list(lulog,form_bhb,op_info)

      call quit(1,'test_formgen2','test exit')

      return
      end
