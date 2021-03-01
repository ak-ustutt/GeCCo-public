*----------------------------------------------------------------------*
      subroutine print_itf(itflog,fl_head,itin,op_info,print_form,
     &                     formlog,tasks,itf_names)
*----------------------------------------------------------------------*
*     Print ITF lines to itflog
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &     itflog,                              ! Output file for ITF algo code
     &     formlog                              ! Output of GeCco formulae
      type(formula_item), intent(in), target ::
     &     fl_head                              ! Linked list of formulae
      type(operator_info), intent(in) ::
     &     op_info                              ! Operator info for printing formulae
      type(tensor_names), intent(in) ::
     &     itf_names                            ! contains renaming information
      logical, intent(in) ::
     &     itin,                                ! Create ITIN lines or symmetrise residual at the end
     &     print_form,                          ! Print to optional formulae file
     &     tasks

      type(formula_item), pointer ::
     &     fl_item                              ! Current formula_item
      integer ::
     &     counter(4), cnt                      ! Counter array, 1 helper variable
      integer ::
     &   inter_itype(MAXINT,INDEX_LEN),         ! Store intermediate index-type (itype) info from previous line
     &   i
      type(x_inter) ::
     &   x_dict(MAXX)
      type(inter_spin_cases) ::
     &   inter_spin_dict        ! Store intermediate spin cases from previous line


      ! Point to start of linked list of formulae
      fl_item => fl_head

      counter(1) = 1    ! contraction number
      counter(2) = 1    ! intermediate number
      counter(3) = 1    ! k4e
      counter(4) = 1    ! x intermediate
      inter_itype = 0

      do i = 1, MAXX
         x_dict(i)%label = ''
         x_dict(i)%ops = 0
      end do

      do i = 1, MAX_SPIN_CASES
         inter_spin_dict%names(i) = ''
      end do
      inter_spin_dict%ncase = 0

      ! Loop over formula_items, end of the list points to NULL
      do while (associated(fl_item%next))

         if (fl_item%command==command_add_intm .or.
     &       fl_item%command==command_cp_intm .or.
     &       fl_item%command==command_add_bc .or.
     &       fl_item%command==command_bc .or.
     &       fl_item%command==command_bc_reo) then

            if (fl_item%command==command_bc_reo) then
              ! patch the result of the contraction with the reordered intermediate info
              ! itf_index_info is already set up correctly
              if (fl_item%bcontr%nj_res < fl_item%reo%nj_out)
     &             call quit(1,'print_itf','dimension trouble')
              fl_item%bcontr%nj_res = fl_item%reo%nj_out
              fl_item%bcontr%occ_res(:,:,1:fl_item%reo%nj_out) =
     &             fl_item%reo%occ_opout(:,:,1:fl_item%reo%nj_out)
            end if

            call command_to_itf(fl_item%bcontr,itin,itflog,
     &                          fl_item%command, inter_itype,
     &                          itf_names,
     &                          counter,tasks,x_dict,
     &                          inter_spin_dict)

            ! Count the number of terms
            counter(1) = counter(1) + 1

         else if (fl_item%command==command_add_contribution) then
            write(itflog,*) '[CONTR]',fl_item%target  ! this case should not appear in an opt. formula
         else if (fl_item%command==command_add_bc_reo) then
            write(itflog,*) '[CONTRACT][REORDER][ADD]',
     &           fl_item%target
            call prt_bcontr(itflog,fl_item%bcontr)
            call prt_reorder(itflog,fl_item%reo)
            call warn('print_itf',
     &           'uncovered case appeared: [CONTRACT][REORDER][ADD] ')
         else if (fl_item%command==command_add_reo) then
            write(itflog,*) '[REORDER][ADD]',
     &           fl_item%target
            call prt_bcontr(itflog,fl_item%bcontr)
            call prt_reorder(itflog,fl_item%reo)
            call warn('print_itf',
     &           'uncovered case appeared: [REORDER][ADD] ')
         else if (fl_item%command==command_reorder) then
            write(itflog,*) '[REORDER]',
     &           fl_item%target
            call prt_reorder(itflog,fl_item%reo)
            call warn('print_itf',
     &           'uncovered case appeared: [REORDER] ')
         else if (fl_item%command==command_symmetrise) then
            write(itflog,*) '[SYMMETRISE]',fl_item%target
            call warn('print_itf',
     &           'uncovered case appeared: [SYMMETRISE] ')
         else if (fl_item%command==command_end_of_formula .or.
     &            fl_item%command==command_set_target_init .or.
     &            fl_item%command==command_set_target_update .or.
     &            fl_item%command==command_new_intermediate .or.
     &            fl_item%command==command_del_intermediate ) then
            ! Do nothing
            ! write(itflog,*) '[END]'
            ! write(itflog,*) '[INIT TARGET]',fl_item%target
            ! write(itflog,*) '[SET TARGET]',fl_item%target
            ! write(itflog,*) '[NEW INTERMEDIATE]',fl_item%target
            ! write(itflog,*) '[DELETE INTERMEDIATE]',fl_item%target
            ! write(itflog,*) '[REORDER]',fl_item%target
         else
            write(itflog,*) 'unknown command ',fl_item%command,
     &           fl_item%target
         end if

         ! Optionally print the formula items to another output file
         if (print_form) then
          write(formlog,*) "FORMULA NUMBER: ", counter(1)
          cnt = counter(1) ! beware of print_form_item2: it increments the counter
          call print_form_item2(formlog,'LONG',cnt,fl_item,
     &                          op_info)
         end if

         ! Check if at the end of the list, if not, point to the next item
         fl_item => fl_item%next

      end do

      return
      end
