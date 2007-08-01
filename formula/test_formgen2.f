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
     &     ndef,
     &     idxopa, idxopb, idxopc, idxopd, idxopr, idxoph,
     &     idxham, idxecc, idxopi, idxopf1, idxopf2, idxopu1

      integer, pointer ::
     &     occ_def(:,:,:)
      
      type(formula_item), target ::
     &     form_a_dc, form_bhb

      type(formula_item), pointer ::
     &     fl_a_dc_pnt, fl_bhb_pnt

      type(operator_array), pointer ::
     &     defop(:)

      integer, external ::
     &     idx_oplist2


      idxecc = idx_oplist2(op_ccen,op_info)
      idxham = idx_oplist2(op_ham,op_info)
      
      write(luout,*) '-------------------'
      write(luout,*) ' adding F1 operator'
      write(luout,*) '-------------------'
      call add_operator('F1',op_info)
      idxopf1 = idx_oplist2('F1',op_info)
      op => op_info%op_arr(idxopf1)%op

      call set_hop(op,'F1',.false.,0,0,1,0,0,
     &     1,1,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding F2 operator'
      write(luout,*) '-------------------'
      call add_operator('F2',op_info)
      idxopf2 = idx_oplist2('F2',op_info)
      op => op_info%op_arr(idxopf2)%op

      call set_hop(op,'F2',.false.,0,0,1,0,0,
     &     2,2,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding U1 operator'
      write(luout,*) '-------------------'
      call add_operator('U1',op_info)
      idxopu1 = idx_oplist2('U1',op_info)
      op => op_info%op_arr(idxopu1)%op

      allocate(occ_def(ngastp,2,3))
      occ_def(1:ngastp,1,1) = (/1,1,0,0/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/1,0,0,1/)
      occ_def(1:ngastp,2,2) = (/2,0,0,0/)
      occ_def(1:ngastp,1,3) = (/0,1,0,1/)
      occ_def(1:ngastp,2,3) = (/2,0,0,0/)

      call set_uop(op,'U1',.false.,0,0,1,0,0,
     &     occ_def,3,orb_info)

      write(luout,*) '-----------------------'
      write(luout,*) ' joining F1+F2 operator'
      write(luout,*) '-----------------------'

      call join_operator(op_info%op_arr(idxopf1)%op,
     &                   op_info%op_arr(idxopf2)%op,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding I operator'
      write(luout,*) '-------------------'
      call add_operator('INT',op_info)
      idxopi = idx_oplist2('INT',op_info)
      op => op_info%op_arr(idxopi)%op

      allocate(defop(3))

      defop(1)%op => op_info%op_arr(idxecc)%op
      defop(2)%op => op_info%op_arr(idxopf1)%op
      defop(3)%op => op_info%op_arr(idxopf2)%op
 
      call set_gen_intermediate(op,'INT',
     &     defop,3,orb_info)

      call quit(1,'test_formgen','exit 0')

      write(luout,*) '-------------------'
      write(luout,*) ' adding A operator'
      write(luout,*) '-------------------'
      call add_operator('A',op_info)

      idxopa = op_info%nops
      op => op_info%op_arr(idxopa)%op
      call set_xop(op,'A',.false.,0,0,1,0,0,
     &     2,2,0,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding B operator'
      write(luout,*) '-------------------'
      call add_operator('B',op_info)

      idxopb = op_info%nops
      op => op_info%op_arr(idxopb)%op
      call set_xop(op,'B',.false.,0,0,1,0,0,
     &     2,2,0,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding C operator'
      write(luout,*) '-------------------'
      call add_operator('C',op_info)

      idxopc = op_info%nops
      op => op_info%op_arr(idxopc)%op
      call set_xop(op,'C',.false.,0,0,1,0,0,
     &     2,2,0,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding R operator'
      write(luout,*) '-------------------'
      call add_operator('R',op_info)

      idxopr = op_info%nops
      op => op_info%op_arr(idxopr)%op
      call set_xop(op,'R',.false.,0,0,1,0,0,
     &     2,2,0,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding D operator'
      write(luout,*) '-------------------'
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

      write(luout,*) '-------------------'
      write(luout,*) ' generate B = A+DC'
      write(luout,*) '-------------------'
      ! define: B = A + DC
      call init_formula(form_a_dc)
      fl_a_dc_pnt => form_a_dc
      call new_formula_item(fl_a_dc_pnt,command_set_target_init,idxopb)
      fl_a_dc_pnt => fl_a_dc_pnt%next
      call expand_op_product(fl_a_dc_pnt,idxopb,
     &     1d0,1,idxopa,-1,-1,
     &     0,0,op_info)
      do while(associated(fl_a_dc_pnt%next))
        fl_a_dc_pnt => fl_a_dc_pnt%next
      end do
      call expand_op_product(fl_a_dc_pnt,idxopb,
     &     1d0,2,(/idxopd,idxopc/),-1,-1,
     &     (/1,2/),1,op_info)

      call print_form_list(luout,form_a_dc,op_info)

      write(luout,*) '-----------------------------'
      write(luout,*) ' generate R = e^{-B} H e^{B}'
      write(luout,*) '-----------------------------'
      ! generate: e^{-B} H e^{B}
      idxoph = idx_oplist2(op_ham,op_info)
      call init_formula(form_bhb)
      fl_bhb_pnt => form_bhb
      call new_formula_item(fl_bhb_pnt,command_set_target_init,idxopr)
      fl_bhb_pnt => fl_bhb_pnt%next
      
      call expand_op_bch(fl_bhb_pnt,4,idxopr,
     &     1d0,-1,idxoph,1d0,idxopb,-1,-1,op_info)

      call print_form_list(luout,form_bhb,op_info)

      ! replace B -> A + DC
      write(luout,*) '------------------------------------'
      write(luout,*) ' insert B: R = e^{-A-DC} H e^{B+DC}'
      write(luout,*) '------------------------------------'

      call expand_subexpr(form_bhb,form_a_dc,op_info)

      call print_form_list(luout,form_bhb,op_info)

      write(luout,*) '------------------------------------'
      write(luout,*) ' after summing terms: '
      write(luout,*) '------------------------------------'

      call sum_terms(form_bhb,op_info)

      call print_form_list(luout,form_bhb,op_info)

      call quit(1,'test_formgen','test exit')

      return
      end
