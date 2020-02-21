*----------------------------------------------------------------------*
      pure function check_inter(label)
*----------------------------------------------------------------------*
!    Check if tensor is an intermediate
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &     label
      logical ::
     &     check_inter

      ! Assume these are the names of intermediates
      if (index(label, "STIN")>0 .or.
     &    index(label, "LTIN")>0) then
         check_inter=.true.
      else
         check_inter=.false.
      end if

      end function


*----------------------------------------------------------------------*
      subroutine print_itf(itflog,fl_head,itin,op_info,print_form,
     &                     formlog)
*----------------------------------------------------------------------*
*     Print ITF info to itflog
*----------------------------------------------------------------------*
      ! TODO: fix this...
      !use itf_utils !copied above from module
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &     itflog,            ! Output file for ITF algo code
     &     formlog            ! Output of GeCco formulae
      type(formula_item), intent(in), target ::
     &     fl_head            ! Linked list of formulae
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in) ::
     &     itin,           ! Create ITIN lines or symmetrise residual at the end
     &     print_form        ! Print to optional formulae file

      type(formula_item), pointer ::
     &     fl_item,   ! Current formula_item
     &     inter_start,    ! Mark start of intermediate search
     &     res_start,      ! Mark end of intermediate search
     &     summed_inter    ! Mark intermediate points in intermediate search
      type(spin_cases), dimension(MAXINT) ::
     &     spin_inters , ! Array of intermeidates with associated spin cases
     &     ospin_inters  ! Numerically ordered array of intermeidates with associated spin cases
      integer ::
     &     i,j,k,     ! Loop indcies
     &     contr_no, ! Counter of contrations
     &     kk,
     &     p_count
      logical ::
     &     check_inter,    ! Need to use itf module instead
     &     more_inter,     ! Check if more intermediates are needed
     &     finished_inter, ! Check if finished recursive intermediate search
     &     symm_res        ! True is intermediate contributes to a symmmetric residual
      integer ::
     &     tmp_case(INDEX_LEN),
     &     ninter,       ! Number of intermediates found in recursive search
     &     shift         ! Used to store sequentially
      character ::
     &     ch           ! Scratch
      integer ::
     &   inter_itype(INDEX_LEN)  ! Store inter itype info from previous line


      ! Point to start of linked list
      fl_item => fl_head
      contr_no = 0
      inter_itype = 0

      ! Loop over formula_items, end of the list points to NULL
      do while (associated(fl_item%next))

         if (fl_item%command==command_add_intm .or.
     &       fl_item%command==command_cp_intm .or.
     &       fl_item%command==command_add_bc .or.
     &       fl_item%command==command_bc .or.
     &       fl_item%command==command_bc_reo) then
            call command_to_itf(fl_item%bcontr,itin,itflog,
     &                          fl_item%command, inter_itype)
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
        call print_form_item2(formlog,'LONG',contr_no,fl_item,op_info)
      end if

      ! Check if at the end of the list, if not, point to the next item
      fl_item => fl_item%next

      ! Count the number of terms
      contr_no = contr_no+1

      end do

      return
      end
