      subroutine test_formgen3(op_info,orb_info)

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
     &     idxecc, idxham, idxopr, idxoprd, idxopd1, idxopb,
     &     idxopf, idxopc, idxopcd

      integer, pointer ::
     &     occ_def(:,:,:)
      
      type(formula_item), target ::
     &     form_rdr, form_b, form_test

      type(formula_item), pointer ::
     &     fl_rdr_pnt, fl_b_pnt, fl_test_pnt

      type(operator_array), pointer ::
     &     defop(:)

      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      integer, external ::
     &     idx_oplist2

      call atim_csw(cpu0,sys0,wall0)
      
      idxecc = idx_oplist2(op_ccen,op_info)
      idxham = idx_oplist2(op_ham,op_info)

      allocate(occ_def(ngastp,2,3))
      
      write(luout,*) '-------------------'
      write(luout,*) ' adding R operator '
      write(luout,*) '-------------------'
      call add_operator('R',op_info)
      idxopr = idx_oplist2('R',op_info)
      op => op_info%op_arr(idxopr)%op

      occ_def(1:ngastp,1,1) = (/0,1,0,1/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)
      occ_def(1:ngastp,1,2) = (/0,0,0,2/)
      occ_def(1:ngastp,2,2) = (/2,0,0,0/)

      call set_uop(op,'R',.false.,0,0,1,0,0,
     &     occ_def,2,orb_info)

      write(luout,*) '--------------------'
      write(luout,*) ' adding R+ operator '
      write(luout,*) '--------------------'
      call add_operator('R+',op_info)
      idxoprd = idx_oplist2('R+',op_info)
      op => op_info%op_arr(idxoprd)%op

      call clone_operator(op,op_info%op_arr(idxopr)%op,orb_info)
      op%dagger = .true.

      write(luout,*) '-------------------'
      write(luout,*) ' adding C operator'
      write(luout,*) '-------------------'
      call add_operator('C',op_info)
      idxopc = idx_oplist2('C',op_info)
      op => op_info%op_arr(idxopc)%op

      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/2,0,0,0/)

      call set_uop(op,'C',.false.,0,0,1,0,0,
     &     occ_def,1,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding C+ operator'
      write(luout,*) '-------------------'
      call add_operator('C+',op_info)
      idxopcd = idx_oplist2('C+',op_info)
      op => op_info%op_arr(idxopcd)%op

      call clone_operator(op,op_info%op_arr(idxopc)%op,orb_info)
      op%dagger = .true.

      ! a dummy operator, clone of C
      write(luout,*) '-------------------'
      write(luout,*) ' adding D1 operator'
      write(luout,*) '-------------------'
      call add_operator('D1',op_info)
      idxopd1 = idx_oplist2('D1',op_info)
      op => op_info%op_arr(idxopd1)%op

      call clone_operator(op,op_info%op_arr(idxopc)%op,orb_info)

      ! fock operator
      write(luout,*) '-------------------'
      write(luout,*) ' adding F operator'
      write(luout,*) '-------------------'
      call add_operator('F',op_info)
      idxopf = idx_oplist2('F',op_info)
      op => op_info%op_arr(idxopf)%op

      occ_def(1:ngastp,1,1) = (/0,1,0,0/)
      occ_def(1:ngastp,2,1) = (/0,0,0,1/)
      occ_def(1:ngastp,1,2) = (/0,0,0,1/)
      occ_def(1:ngastp,2,2) = (/0,1,0,0/)
      occ_def(1:ngastp,1,3) = (/0,0,0,1/)
      occ_def(1:ngastp,2,3) = (/0,0,0,1/)

      call set_uop(op,'F',.false.,0,0,1,0,0,
     &     occ_def,3,orb_info)
c      call set_hop(op,'F',.false.,0,0,1,0,0,
c     &     1,1,iformal,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding B operator'
      write(luout,*) '-------------------'
      call add_operator('B',op_info)
      idxopb = idx_oplist2('B',op_info)
      op => op_info%op_arr(idxopb)%op

      allocate(defop(2))

      defop(1)%op => op_info%op_arr(idxecc)%op
      defop(2)%op => op_info%op_arr(idxopd1)%op
 
      call set_gen_intermediate(op,'B',
     &     defop,2,orb_info)

      call write_title(luout,wst_dbg_func,'making B formula')

      ! define  B = R+ D1 R
      call init_formula(form_rdr)
      fl_rdr_pnt => form_rdr
      call new_formula_item(fl_rdr_pnt,command_set_target_init,idxecc)
      fl_rdr_pnt => form_rdr%next
      call expand_op_product(fl_rdr_pnt,idxecc,
     &     1d0,4,(/idxoprd,idxopd1,idxopf,idxopr/),-1,-1,
     &     (/1,3,3,4/),2,op_info)

      call print_form_list(luout,form_rdr,op_info)

      call form_deriv3(form_b,form_rdr,
     &     1,idxopd1,0,idxopb,
     &     op_info)

      ! so we have B
      call write_title(luout,wst_dbg_func,'B formula')
      call print_form_list(luout,form_b,op_info)

      ! make a more complicated formula
      call init_formula(form_test)
      fl_test_pnt => form_test
      call new_formula_item(fl_test_pnt,command_set_target_init,idxecc)
      fl_test_pnt => form_test%next
      call expand_op_product(fl_test_pnt,idxecc,
     &     1d0,5,(/idxoprd,idxopcd,idxopf,idxopc,idxopr/),-1,-1,
     &     (/1,2,1,3,3,5,4,5/),4,op_info)

      call write_title(luout,wst_dbg_func,'test formula')
      call print_form_list(luout,form_test,op_info)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'1st part',cpu-cpu0,sys-sys0,wall-wall0)
      
      call atim_csw(cpu0,sys0,wall0)

      call write_title(luout,wst_dbg_func,'now trying to factorize')

      call factor_out_subexpr(form_test,form_b,op_info)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'factor',cpu-cpu0,sys-sys0,wall-wall0)

      call print_form_list(luout,form_test,op_info)

      call quit(1,'test_formgen','exit 0')

      return
      end
