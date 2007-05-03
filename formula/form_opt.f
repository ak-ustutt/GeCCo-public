*----------------------------------------------------------------------*
      subroutine form_opt(ffform_opt,
     &     nfcat,idxform,
     &     ffform,nform,
     &     op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     given a list of formulae, concatenate them into one formula
*     file, find optimal factorization and intermediates
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula.h'

      integer, intent(in) ::
     &     nfcat, nform, idxform(nfcat)
      type(filinf), intent(inout) ::
     &     ffform_opt
      type(file_array), intent(inout) ::
     &     ffform(*)
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      
      type(formula), pointer ::
     &     form_head, form_tail, form_ptr
      type(filinf), pointer ::
     &     cur_ffile

      integer ::
     &     icat, iprint

      logical, external ::
     &     rd_formula

      iprint = max(ntest,iprlvl)

      ! initialize list
      allocate(form_head)
      form_ptr => form_head
      form_ptr%command = command_end_of_formula
      nullify(form_ptr%next)
      nullify(form_ptr%contr)
      nullify(form_ptr%interm)

      ! ----------------------
      ! read in formula files:
      ! ----------------------
      do icat = 1, nfcat

        cur_ffile => ffform(idxform(icat))%fhand

        call read_form_list(cur_ffile,form_ptr)

        ! advance form_ptr to end of list
        ! (not possible via call list due to ifort problems)
        do while(associated(form_ptr%next))
          form_ptr => form_ptr%next
        end do
      end do
      form_tail => form_ptr

      ! ----------------------------------------
      ! find optimal factorization for each term
      ! ----------------------------------------
c      form_ptr => form_head
      call factorize(form_head,op_info,str_info,orb_info)

      if (iprint.ge.10) then
        write(luout,*) 'Optimized formula:'
        write(luout,*) '=================='
        call print_form_list(luout,form_head,op_info)
      end if

      call write_form_list(ffform_opt,form_head)
      
      return
      end
