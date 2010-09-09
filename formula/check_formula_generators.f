      subroutine check_formula_generators(form,op_info,orb_info)

      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'def_contraction.h'
      include 'def_formula.h'
      include 'def_formula_item.h'
      include 'ifc_operators.h'
      include 'par_opnames_gen.h'

      type(formula) ::
     &     form
      type(operator_info) ::
     &     op_info
      type(orbinf) ::
     &     orb_info
      type(contraction) ::
     &     ctest, ctest2

      type(operator), pointer ::
     &     opa, opb, opc, opd, ope, opf, opg, opr, opr2

      integer, parameter ::
     &     maxdef = 20
      integer ::
     &     occ_def(ngastp,2,maxdef)

      integer ::
     &     idxopa, idxopb, idxopc, idxopd, idxope, idxopf, idxopg,
     &     idxopr, idxoph, iformal, idxopr2
      
      type(formula_item), target ::
     &     form_r, form_r2, form_c

      type(formula_item), pointer ::
     &     fl_a_dc_pnt, fl_bhb_pnt

      type(operator_array), pointer ::
     &     defop(:)
      integer, external ::
     &     idx_oplist2
      logical, external ::
     &     contr_in_contr


      call write_title(luout,wst_dbg_subr,
     &     'Checking Formula Generator Routines')

      write(luout,*) 'I will define a few operators for testing'

      write(luout,*) '-------------------'
      write(luout,*) ' adding A operator'
      write(luout,*) '-------------------'
      call add_operator('A',op_info)

      iformal=3
      idxopa = idx_oplist2('A',op_info)
      opa => op_info%op_arr(idxopa)%op
      call set_xop(opa,'A',.false.,
     &     2,2,0,0,iformal,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding B operator'
      write(luout,*) '-------------------'
      call add_operator('B',op_info)

      idxopb =  idx_oplist2('B',op_info)
      opb => op_info%op_arr(idxopb)%op
      occ_def(1:ngastp,1,1) = (/2,0,0,0/)
      occ_def(1:ngastp,2,1) = (/0,2,0,0/)
      call set_uop(opb,'B',.false.,
     &     occ_def,1,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding C operator'
      write(luout,*) '-------------------'
      call add_operator('C',op_info)

      idxopc =  idx_oplist2('C',op_info)
      opc => op_info%op_arr(idxopc)%op
      occ_def(1:ngastp,1,1) = (/0,1,0,0/)
      occ_def(1:ngastp,2,1) = (/0,1,0,0/)
      call set_uop(opc,'C',.false.,
     &     occ_def,1,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding D operator'
      write(luout,*) '-------------------'
      call add_operator('D',op_info)

      idxopd =  idx_oplist2('D',op_info)
      opd => op_info%op_arr(idxopd)%op
      occ_def(1:ngastp,1,1) = (/1,0,0,0/)
      occ_def(1:ngastp,2,1) = (/0,1,0,0/)
      call set_uop(opd,'D',.false.,
     &     occ_def,1,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding E operator'
      write(luout,*) '-------------------'
      call add_operator('E',op_info)
      idxope =  idx_oplist2('E',op_info)
      ope => op_info%op_arr(idxope)%op
      occ_def(1:ngastp,1,1) = (/0,2,0,0/)
      occ_def(1:ngastp,2,1) = (/0,2,0,0/)
      call set_uop(ope,'E',.false.,
     &     occ_def,1,orb_info)

      write(luout,*) '-------------------'
      write(luout,*) ' adding R operator'
      write(luout,*) '-------------------'
      call add_operator('R',op_info)

      idxopr =  idx_oplist2('R',op_info)
      opr => op_info%op_arr(idxopr)%op

      allocate(defop(2))
      defop(1)%op => opa
      defop(2)%op => opc

      call set_gen_intermediate(opr,'R',
     &     defop,2,orb_info)

      deallocate(defop)

      call list_operators(luout,op_info)

      write(luout,*) 'now we generate: R = B B B x A A x'

      call init_formula(form_r)

      call expand_op_product2(form_r,idxopr,
     &     1d0,9,6,
     &     (/idxopr,idxopb,idxopb,idxopr,
     &                            idxopr,idxopa,idxopa,idxopa,idxopr/),
     &     (/1     ,2     ,3     ,1     ,
     &                            1     ,4     ,5     ,6,     1     /),
     &     -1,-1,
     &     (/2,7,2,8,3,7,3,8/),4,
     &     0,0,
     &     .false.,op_info)

c      write(luout,*) 'expansion 1 (result: r)'
c      call print_form_list(luout,form_r,op_info)
c      call dealloc_formula_list(form_r)
c      call init_formula(form_c)
c
c      call expand_op_product2(form_c,idxopc,
c     &     1d0,9,8,
c     &     (/idxopc,idxopb,idxopb,idxopb,
c     &                            idxope,idxopa,idxopa,idxopa,idxopc/),
c     &     (/1     ,2     ,3     ,4     ,
c     &                            5     ,6     ,7     ,8,     1     /),
c     &     -1,-1,
c     &     (/1,6, 2,5, 3,6, 4,7, 5,6, 2,8, 3,9/),7,
c     &     0,0,
c     &     .false.,op_info)
c
c      write(luout,*) 'expansion 2 (result: c)'
c      call print_form_list(luout,form_c,op_info)
c
c      call add_operator('R2',op_info)
c
c      idxopr2 = idx_oplist2('R2',op_info)
c      opr2 => op_info%op_arr(idxopr2)%op
c
c      occ_def(1:ngastp,1,1) = (/0,0,0,0/)
c      occ_def(1:ngastp,2,1) = (/0,2,0,0/)
c      occ_def(1:ngastp,1,2) = (/0,2,0,0/)
c      occ_def(1:ngastp,2,2) = (/0,0,0,0/)
c      call set_uop2(opr2,'R2',
c     &     occ_def,1,2,(/.true.,.true./),orb_info)
c
c      call init_formula(form_r2)
cc      call expand_op_product2(form_r2,idxopr2,
cc     &     1d0,6,3,
cc     &     (/idxopr2,idxopb,idxopr2,idxopr2,idxopa,idxopr2/),
cc     &     (/1     ,2     ,1     ,1  , 3     ,1/),
cc     &     -1,-1,
cc     &     (/2,5/),1,
cc     &     0,0,
cc     &     .false.,op_info)
c      call expand_op_product2(form_r2,idxopr2,
c     &     1d0,8,5,
c     &     (/idxopr2,idxopb,idxopb,idxopr2,
c     &                         idxopr2,idxopa,idxopa,idxopr2/),
c     &     (/1     ,2     ,4,1     ,1  , 3,5     ,1/),
c     &     -1,-1,
c     &     (/2,6,2,7,3,6,3,7/),4,
c     &     0,0,
c     &     .false.,op_info)
c
c      write(luout,*) 'expansion 3 (result: r2)'
c      call print_form_list(luout,form_r2,op_info)
c
c      write(luout,*) 'trying to split'
c
c      if (.not.associated(form_c%contr)) stop '1'
c      if (.not.associated(form_r2%contr)) stop '2'
c
c      write(luout,*) 'contr_in_contr(1): ',
c     &     contr_in_contr(form_r2%contr,form_c%contr)
c      write(luout,*) 'contr_in_contr(2): ',
c     &     contr_in_contr(form_r2%next%contr,form_c%contr)
c      write(luout,*) 'contr_in_contr(3): ',
c     &     contr_in_contr(form_r2%next%next%contr,form_c%contr)
c
c      call init_contr(ctest)
c      call split_contr2(ctest,form_r2%contr,
c     &     form_c%contr,
c     &     op_info)
c
c      call init_contr(ctest2)
c      call join_contr3(ctest2,ctest,form_r2%contr,
c     &     idxopr2,1,op_info)

      return
      end
