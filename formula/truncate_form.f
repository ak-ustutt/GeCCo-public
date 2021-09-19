      subroutine truncate_form(flist,type,op_info,form_info)
*----------------------------------------------------------------------*
*     Routine used to set up the deletion parameters for truncating a 
*     formula. Initially used to truncate the R12-amplitude equations in 
*     the CCSD-R12 model to produce the CCSD(R12) approach.
*     Should be extendable later to produce CCSD(T) etc. type ansaetze.
*     GWR February 2008
*----------------------------------------------------------------------*

      implicit none
      
      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'mdef_formula_info.h'
      include 'def_formula_item.h'
      include 'def_del_list.h'
      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'

      type(formula_item), intent(inout), target::
     &     flist
      type(formula_info), intent(inout) ::
     &     form_info
      integer, intent(in)::
     &     type
      type(operator_info), intent(in) ::
     &     op_info

      type(del_cond_list) ::
     &     del_list
      type(del_cond), pointer ::
     &     temp_del
      type(formula_item) ::
     &     fl_intm
      type(filinf), pointer ::
     &     ffintm
      integer ::
     &     idx,narrays,nrpl,nspl
      integer, external ::
     &     idx_formlist

      if(ntest.ge.100)then
        write(lulog,*) '--------------------'
        write(lulog,*) ' Formula truncation '
        write(lulog,*) '--------------------'
      endif

      if(type.eq.0)then
        ! CCSD(R12).

        allocate(del_list%del_cond_item(3))
        del_list%or_dim = 3
        allocate(del_list%del_cond_item(1)%del_cond_arr(2),
     &       del_list%del_cond_item(2)%del_cond_arr(3),
     &       del_list%del_cond_item(3)%del_cond_arr(3))
        del_list%del_cond_item(1)%and_dim = 1
c        del_list%del_cond_item(1)%and_dim = 2
        del_list%del_cond_item(2)%and_dim = 3
        del_list%del_cond_item(3)%and_dim = 3

        ! Array 1.
c        temp_del => del_list%del_cond_item(1)%del_cond_arr(1)
c        temp_del%op_name = op_cba
c        temp_del%num_op_restr(1:2) = 1
c        temp_del%part_num_restr(1:2) = -1

c        temp_del => del_list%del_cond_item(1)%del_cond_arr(2)
        temp_del => del_list%del_cond_item(1)%del_cond_arr(1)
c        temp_del%op_name = op_c12
        temp_del%op_name = op_r12
        temp_del%transposed = .false.
        temp_del%num_op_restr(1) = 2
        temp_del%num_op_restr(2) = -1
        temp_del%part_num_restr(1:2) = -1
      
        ! Array 2.
        temp_del => del_list%del_cond_item(2)%del_cond_arr(1)
c       temp_del%op_name = op_cba
        temp_del%op_name = op_r12
        temp_del%transposed = .true.
        temp_del%num_op_restr(1:2) = 1
        temp_del%part_num_restr(1:2) = -1

        temp_del => del_list%del_cond_item(2)%del_cond_arr(2)
c     temp_del%op_name = op_c12
        temp_del%op_name = op_r12
        temp_del%transposed = .false.
        temp_del%num_op_restr(1:2) = 1
        temp_del%part_num_restr(1:2) = -1

        temp_del => del_list%del_cond_item(2)%del_cond_arr(3)
        temp_del%op_name = op_ham
        temp_del%transposed = .false.
        temp_del%num_op_restr(1:2) = 1
        temp_del%part_num_restr(1:2) = 2

        ! Array 3.
        temp_del => del_list%del_cond_item(3)%del_cond_arr(1)
c     temp_del%op_name = op_cba
        temp_del%op_name = op_r12
        temp_del%transposed = .true.
        temp_del%num_op_restr(1:2) = 1
        temp_del%part_num_restr(1:2) = -1

        temp_del => del_list%del_cond_item(3)%del_cond_arr(2)
c     temp_del%op_name = op_c12
        temp_del%op_name = op_r12
        temp_del%transposed = .false.
        temp_del%num_op_restr(1:2) = 1
        temp_del%part_num_restr(1:2) = -1

        temp_del => del_list%del_cond_item(3)%del_cond_arr(3)
        temp_del%op_name = op_top
        temp_del%transposed = .false.
        temp_del%num_op_restr(1) = 1
        temp_del%num_op_restr(2) = -1
        temp_del%part_num_restr(1:2) = -1

        narrays = 3

      elseif (type.eq.1)then
        ! Deletion of terms with 2 F12 vertices only.

        allocate(del_list%del_cond_item(1))
        del_list%or_dim = 1

        allocate(del_list%del_cond_item(1)%del_cond_arr(1))
        del_list%del_cond_item(1)%and_dim = 1
        temp_del => del_list%del_cond_item(1)%del_cond_arr(1)

        ! Properties of relevant operator.
        temp_del%op_name = op_r12
        temp_del%transposed = .false.
        temp_del%num_op_restr(1) = 2
        temp_del%num_op_restr(2) = -1
        temp_del%part_num_restr(1:2) = -1

        narrays = 1

      elseif(type.eq.2)then
        call quit(1,'truncate_form','do not use this method')
        ! Deletion of terms factorised by Z.

        ! Factor out Z first.
        idx = idx_formlist(form_r12_zint,form_info)
        if(idx.lt.0)
     &       call quit(1,'truncate_form','Z not on list')
        ffintm => form_info%form_arr(idx)%form%fhand
        call init_formula(fl_intm)
        call read_form_list(ffintm,fl_intm,.true.)
        call factor_out_subexpr2(flist,fl_intm,.false.,
     &                           nrpl,nspl,op_info)

        ! Now delete all terms including Z.
        allocate(del_list%del_cond_item(1))
        del_list%or_dim = 1

        allocate(del_list%del_cond_item(1)%del_cond_arr(1))
        del_list%del_cond_item(1)%and_dim = 1
        temp_del => del_list%del_cond_item(1)%del_cond_arr(1)

        ! Properties of relevant operator.
        temp_del%op_name = op_z_inter
        temp_del%transposed = .false.
        temp_del%num_op_restr(1) = 2
        temp_del%num_op_restr(2) = -1
        temp_del%part_num_restr(1:2) = -1

        narrays = 1

      elseif(type.eq.3)then
        call quit(1,'truncate_form','do not use this method')
        ! Deletion of terms factorised by Z and of terms containing
        ! 2 or more F12 terms.

        ! Factor out Z first.
        idx = idx_formlist(form_r12_zint,form_info)
        if(idx.lt.0)
     &       call quit(1,'truncate_form','Z not on list')
        ffintm => form_info%form_arr(idx)%form%fhand
        call init_formula(fl_intm)
        call read_form_list(ffintm,fl_intm,.true.)
        call factor_out_subexpr2(flist,fl_intm,.false.,
     &                           nrpl,nspl,op_info)

        ! Now delete all terms including Z and with >=2 F12 terms.
        allocate(del_list%del_cond_item(2))
        del_list%or_dim = 2

        allocate(del_list%del_cond_item(1)%del_cond_arr(1),
     &       del_list%del_cond_item(2)%del_cond_arr(1))
        del_list%del_cond_item(1)%and_dim = 1
        del_list%del_cond_item(2)%and_dim = 1

        temp_del => del_list%del_cond_item(1)%del_cond_arr(1)

        ! Properties of relevant operator.
        temp_del%op_name = op_z_inter
        temp_del%transposed = .false.
        temp_del%num_op_restr(1) = 1
        temp_del%num_op_restr(2) = -1
        temp_del%part_num_restr(1:2) = -1

        temp_del => del_list%del_cond_item(2)%del_cond_arr(1)

        ! Properties of relevant operator.
        temp_del%op_name = op_r12
        temp_del%transposed = .false.
        temp_del%num_op_restr(1) = 2
        temp_del%num_op_restr(2) = -1
        temp_del%part_num_restr(1:2) = -1

        narrays = 2

      else
        call quit(1,'truncate_form','unrecognised method')
      endif

      ! Call the truncation routine.
      call form_elem_del(flist,del_list,op_info)

      do idx = 1, narrays
        deallocate(del_list%del_cond_item(idx)%del_cond_arr)
      enddo
      deallocate(del_list%del_cond_item)

      return
      end
