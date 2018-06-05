*----------------------------------------------------------------------*
      subroutine print_itf(itflog,fl_head,op_info,print_form,formlog)
*----------------------------------------------------------------------*
*     Print ITF info to itflog
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_itf_tensor.h'

      integer, intent(in) ::
     &     itflog,
     &     formlog
      type(formula_item), intent(in), target ::
     &     fl_head
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in) ::
     &     print_form        ! Print to optional formulae file
      logical ::
     &     long
      integer ::
     &     inter=0    ! Counter of intermediates
      type(formula_item), pointer ::
     &     fl_item,  ! Current formula_item
     &     next_item,  ! Next formula_item
     &     prev_item, ! Previous formula_item
     &     next_item2,  ! Next formula_item
     &     prev_item2  ! Previous formula_item
      integer ::
     &     nops(4,2)  ! Matrix of index info
      character(len=maxlen_bc_label) ::
     &     old_res='>',     ! Name of tensors involved in the contraction
     &     contract_next='>', ! Name of tensors involved in the contraction
     &     rename_tensor
      integer ::
     &     i,j,      ! loop indcies
     &     contr_no
      character(len=8) ::
     &     istr1='        ',         ! Operator 1 index
     &     istr2='        ',         ! Operator 2 index
     &     istr3='        ',         ! Result index
     &     nstr1='        ',         ! Next operator 1 index
     &     nstr2='        ',         ! Next operator 2 index
     &     nstr3='        ',         ! Next result index
     &     contract_next_index='        ',         ! Next CONTRACT result index
     &     prev_str3='>       '      ! Previous result index
      type(itf_tensor) ::
     &     itf1, itf2, itf3,
     &     p1, p2, p3,
     &     n1, n2, n3

      ! Point to start of linked list
      fl_item=>fl_head
      contr_no=0

      ! Loop over formula_items
      do
      select case(fl_item%command)
      case(command_end_of_formula)
!        write(itflog,*) '[END]'
      case(command_set_target_init)
!        write(itflog,*) '[INIT TARGET]',fl_item%target
      case(command_set_target_update)
!        write(itflog,*) '[SET TARGET]',fl_item%target
      case(command_new_intermediate)
!        write(itflog,*) '[NEW INTERMEDIATE]',fl_item%target
!        write(itflog,'(2x,a)') trim(fl_item%interm%name)
!        write(itflog,'(2x,"attribute parentage: ",a," ",a)')
!     &                        trim(fl_item%parent1),
!     &                        trim(fl_item%parent2)
!        write(itflog,'(2x,"incore: ",i2)') fl_item%incore
!        call print_op_occ(itflog,fl_item%interm)

        inter=inter+1

      case(command_del_intermediate)
!        write(itflog,*) '[DELETE INTERMEDIATE]',fl_item%target
!        write(itflog,'(2x,a)') trim(fl_item%label)
      case(command_add_contribution)
        write(itflog,*) '[CONTR]',fl_item%target
      case(command_add_intm)
        write(itflog,*) '[ADD]',
     &       fl_item%target
        call prt_bcontr(itflog,fl_item%bcontr)
      case(command_cp_intm)
        write(itflog,*) '[COPY]',
     &       fl_item%target
        call prt_bcontr(itflog,fl_item%bcontr)
      case(command_add_bc)
!        write(itflog,*) '[CONTRACT][ADD]',
!     &       fl_item%target
!        call prt_bcontr(itflog,fl_item%bcontr)

        ! Get index for current contraction
        call assign_index(fl_item%bcontr,istr1,istr2,istr3,itflog)
        call assign_spin(istr1,istr2,istr3,
     &                   fl_item%bcontr%label_res,
     &                   fl_item%bcontr%label_op1,
     &                   fl_item%bcontr%label_op2,
     &                   fl_item%bcontr%fact,inter,itflog)

!        call itf_tensor_init(fl_item%bcontr,itf1,itf2,itf3)
        call clear_index(istr1,istr2, istr3)


      case(command_add_bc_reo)
        write(itflog,*) '[CONTRACT][REORDER][ADD]',
     &       fl_item%target
        call prt_bcontr(itflog,fl_item%bcontr)
        call prt_reorder(itflog,fl_item%reo)
      case(command_bc)
!        write(itflog,*) '[CONTRACT]',
!     &       fl_item%target
!        call prt_bcontr(itflog,fl_item%bcontr)

        ! Assuming that this is called only after NEW INTERMEDIATE
        call assign_index(fl_item%bcontr,istr1,istr2,istr3,itflog)
        call assign_spin(istr1,istr2,istr3,
     &                   fl_item%bcontr%label_res,
     &                   fl_item%bcontr%label_op1,
     &                   fl_item%bcontr%label_op2,
     &                   fl_item%bcontr%fact,inter,itflog)

        call clear_index(istr1,istr2, istr3)


      case(command_bc_reo)
!        write(itflog,*) '[CONTRACT][REORDER]',
!     &       fl_item%target
!        call prt_bcontr(itflog,fl_item%bcontr)
!        call prt_reorder(itflog,fl_item%reo)

        ! Get index for current contraction
        call assign_index(fl_item%bcontr,istr1,istr2,istr3,itflog)
        call assign_spin(istr1,istr2,istr3,
     &                   fl_item%bcontr%label_res,
     &                   fl_item%bcontr%label_op1,
     &                   fl_item%bcontr%label_op2,
     &                   fl_item%bcontr%fact,inter,itflog)

        call clear_index(istr1,istr2, istr3)


      case(command_reorder)
!        write(itflog,*) '[REORDER]',fl_item%target
!        call prt_reorder(itflog,fl_item%reo)
      case(command_add_reo)
        write(itflog,*) '[REORDER][ADD]',
     &       fl_item%target
        call prt_bcontr(itflog,fl_item%bcontr)
        call prt_reorder(itflog,fl_item%reo)
      case(command_symmetrise)
        write(itflog,*) '[SYMMETRISE]',fl_item%target
      case default
        write(itflog,*) 'unknown command ',fl_item%command,
     &       fl_item%target
      end select

      if (print_form) then
        call print_form_item2(formlog,'LONG',contr_no,fl_item,op_info)
      end if

      if (.not.associated(fl_item%next)) exit
      fl_item=>fl_item%next

      contr_no=contr_no+1

      end do
      
      end
