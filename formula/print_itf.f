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
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &     itflog,
     &     formlog
      type(formula_item), intent(in), target ::
     &     fl_head
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in) ::
     &     print_form        ! Print to optional formulae file

      integer ::
     &     inter=0    ! Counter of intermediates
      type(formula_item), pointer ::
     &     fl_item   ! Current formula_item
      integer ::
     &     i,j,      ! loop indcies
     &     contr_no,
     &     perm_array(4)
      type(itf_contr) ::
     &     itf_item

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

        ! Get permutaion factors, in an array
        ! Loop over array, send factor to init
        ! Spin sum
        ! Repeat for all permutaion factors

!        perm_array=0
!        call permute_tensors2(fl_item%bcontr,perm_array,itflog)
!        if (perm_array(1)==0) then
!           ! No permutations
!           call itf_contr_init(fl_item%bcontr,itf_item,perm_array(1),
!     &                         itflog)
!           call assign_spin(itf_item)
!        else
!           do i=1, size(perm_array)
!              ! Loop over permuation cases and send seperatley to
!              ! assign_spin
!              call itf_contr_init(fl_item%bcontr,itf_item,perm_array(i),
!     &                            itflog)
!              call assign_spin(itf_item)
!              if (perm_array(i+1)==0) exit
!           end do
!        end if

        call itf_contr_init(fl_item%bcontr,itf_item,1,itflog)
        call assign_spin(itf_item)

      case(command_add_bc_reo)
        write(itflog,*) '[CONTRACT][REORDER][ADD]',
     &       fl_item%target
        call prt_bcontr(itflog,fl_item%bcontr)
        call prt_reorder(itflog,fl_item%reo)
      case(command_bc)
!        write(itflog,*) '[CONTRACT]',
!     &       fl_item%target
!        call prt_bcontr(itflog,fl_item%bcontr)

        call itf_contr_init(fl_item%bcontr,itf_item,1,itflog)
        call assign_spin(itf_item)


      case(command_bc_reo)
!        write(itflog,*) '[CONTRACT][REORDER]',
!     &       fl_item%target
!        call prt_bcontr(itflog,fl_item%bcontr)
!        call prt_reorder(itflog,fl_item%reo)

        call itf_contr_init(fl_item%bcontr,itf_item,1,itflog)
        call assign_spin(itf_item)


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
