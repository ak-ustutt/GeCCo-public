      subroutine truncate_form(flist,op_info)
*----------------------------------------------------------------------*
*     Routine used to set up the deletion parameters for truncating a 
*     formula. Initially used to truncate the R12-amplitude equations in 
*     the CCSD-R12 model to produce the CCSD(R12) approach.
*     Should be extendable later to produce CCSD(T) etc. type ansaetze.
*     GWR February 2008
*----------------------------------------------------------------------*

      implicit none
      
      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'def_del_list.h'
      include 'par_opnames_gen.h'

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info

      type(del_cond_list) ::
     &     del_list
      type(del_cond), pointer ::
     &     temp_del

      if(ntest.ge.100)then
        write(luout,*) '--------------------'
        write(luout,*) ' Formula truncation '
        write(luout,*) '--------------------'
      endif

      ! Initially hard-wired for CCSD(R12).
      ! Will need to generalise this later.
      allocate(del_list%del_cond_item(3))
      del_list%or_dim = 3
      allocate(del_list%del_cond_item(1)%del_cond_arr(2),
     &     del_list%del_cond_item(2)%del_cond_arr(3),
     &     del_list%del_cond_item(3)%del_cond_arr(3))
      del_list%del_cond_item(1)%and_dim = 2
      del_list%del_cond_item(2)%and_dim = 3
      del_list%del_cond_item(3)%and_dim = 3

      ! Array 1.
      temp_del => del_list%del_cond_item(1)%del_cond_arr(1)
      temp_del%op_name = op_rba
      temp_del%num_op_restr(1:2) = 1
      temp_del%part_num_restr(1:2) = -1

      temp_del => del_list%del_cond_item(1)%del_cond_arr(2)
      temp_del%op_name = op_r12
      temp_del%num_op_restr(1) = 2
      temp_del%num_op_restr(2) = -1
      temp_del%part_num_restr(1:2) = -1
      
      ! Array 2.
      temp_del => del_list%del_cond_item(2)%del_cond_arr(1)
      temp_del%op_name = op_rba
      temp_del%num_op_restr(1:2) = 1
      temp_del%part_num_restr(1:2) = -1

      temp_del => del_list%del_cond_item(2)%del_cond_arr(2)
      temp_del%op_name = op_r12
      temp_del%num_op_restr(1:2) = 1
      temp_del%part_num_restr(1:2) = -1

      temp_del => del_list%del_cond_item(2)%del_cond_arr(3)
      temp_del%op_name = op_ham
      temp_del%num_op_restr(1:2) = 1
      temp_del%part_num_restr(1:2) = 2

      ! Array 3.
      temp_del => del_list%del_cond_item(3)%del_cond_arr(1)
      temp_del%op_name = op_rba
      temp_del%num_op_restr(1:2) = 1
      temp_del%part_num_restr(1:2) = -1

      temp_del => del_list%del_cond_item(3)%del_cond_arr(2)
      temp_del%op_name = op_r12
      temp_del%num_op_restr(1:2) = 1
      temp_del%part_num_restr(1:2) = -1

      temp_del => del_list%del_cond_item(3)%del_cond_arr(3)
      temp_del%op_name = op_top
      temp_del%num_op_restr(1) = 1
      temp_del%num_op_restr(2) = -1
      temp_del%part_num_restr(1:2) = -1


      call form_elem_del(flist,del_list,op_info)

      deallocate(del_list%del_cond_item(1)%del_cond_arr,
     &     del_list%del_cond_item(2)%del_cond_arr,
     &     del_list%del_cond_item(3)%del_cond_arr)
      deallocate(del_list%del_cond_item)

      return
      end
