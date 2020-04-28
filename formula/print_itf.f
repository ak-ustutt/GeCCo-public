*----------------------------------------------------------------------*
      subroutine print_itf(itflog,fl_head,itin,op_info,print_form,
     &                     formlog,tasks,taskslog)
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
     &     formlog,                             ! Output of GeCco formulae
     &     taskslog
      type(formula_item), intent(in), target ::
     &     fl_head                              ! Linked list of formulae
      type(operator_info), intent(in) ::
     &     op_info                              ! Operator info for printing formulae
      logical, intent(in) ::
     &     itin,                                ! Create ITIN lines or symmetrise residual at the end
     &     print_form,                          ! Print to optional formulae file
     &     tasks

      type(formula_item), pointer ::
     &     fl_item                              ! Current formula_item
      integer ::
     &     counter(4)                           ! Counter array
!     &     contr_no,                            ! Counter of contrations
!     &     k4e_no,                              ! Counter of K ext contractions (involving 4 ext ints)
!     &     x_no
      integer ::
     &   inter_itype(MAXINT,INDEX_LEN)                 ! Store intermediate index-type (itype) info from previous line


      ! Point to start of linked list of formulae
      fl_item => fl_head

      counter(1) = 1    ! contraction number
      counter(2) = 1    ! intermediate number
      counter(3) = 1    ! k4e
      counter(4) = 1    ! x intermediate
      !contr_no = 1
      !k4e_no = 1
      !x_no = 1

      inter_itype = 0

      ! Loop over formula_items, end of the list points to NULL
      do while (associated(fl_item%next))

         if (fl_item%command==command_add_intm .or.
     &       fl_item%command==command_cp_intm .or.
     &       fl_item%command==command_add_bc .or.
     &       fl_item%command==command_bc .or.
     &       fl_item%command==command_bc_reo) then

!            call command_to_itf(fl_item%bcontr,itin,itflog,
!     &                          fl_item%command, inter_itype,
!     &                          contr_no,k4e_no,tasks,taskslog,
!     &                          x_no)

            call command_to_itf(fl_item%bcontr,itin,itflog,
     &                          fl_item%command, inter_itype,
     &                          counter,tasks,taskslog)

            ! Count the number of terms
            counter(1) = counter(1) + 1

         else if (fl_item%command==command_add_contribution) then
            write(itflog,*) '[CONTR]',fl_item%target
         else if (fl_item%command==command_add_bc_reo) then
            write(itflog,*) '[CONTRACT][REORDER][ADD]',
     &           fl_item%target
            call prt_bcontr(itflog,fl_item%bcontr)
            call prt_reorder(itflog,fl_item%reo)
         else if (fl_item%command==command_add_reo) then
            write(itflog,*) '[REORDER][ADD]',
     &           fl_item%target
            call prt_bcontr(itflog,fl_item%bcontr)
            call prt_reorder(itflog,fl_item%reo)
         else if (fl_item%command==command_symmetrise) then
            write(itflog,*) '[SYMMETRISE]',fl_item%target
         else if (fl_item%command==command_end_of_formula .or.
     &            fl_item%command==command_set_target_init .or.
     &            fl_item%command==command_set_target_update .or.
     &            fl_item%command==command_new_intermediate .or.
     &            fl_item%command==command_del_intermediate .or.
     &            fl_item%command==command_reorder) then
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
          call print_form_item2(formlog,'LONG',counter(1),fl_item,
     &                          op_info)
         end if

         ! Check if at the end of the list, if not, point to the next item
         fl_item => fl_item%next

      end do

      return
      end
