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
     &     itflog,            ! Output file for ITF algo code
     &     formlog            ! Output of GeCco formulae
      type(formula_item), intent(in), target ::
     &     fl_head            ! Linked list of formulae
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in) ::
     &     print_form        ! Print to optional formulae file

      type(formula_item), pointer ::
     &     fl_item,   ! Current formula_item
     &     inter_start,
     &     res_start
      type(spin_cases), dimension(4) ::
     &     spin_inters  ! Array of intermeidates with associated spin cases
      type(itf_intermediate), pointer ::
     &     inter
      type(itf_intermediate_spin), pointer ::
     &     spin
      integer ::
     &     i,j,      ! Loop indcies
     &     contr_no  ! Counter of contrations

      ! Point to start of linked list
      fl_item=>fl_head
      contr_no=0

      !nullify(inter_head)
      !allocate(inter_head)
      !inter_head%val=0

      !tmp_inter => inter_head
      !i = 1
      !do j=0, 4
      !   write(itflog,*) "Value: ", tmp_inter%val
      !   allocate(inter_next)
      !   inter_next%val = i
      !   tmp_inter%next_inter=>inter_next
      !   i = i + 1
      !   tmp_inter = inter_next
      !end do
      !deallocate(inter_head)

      allocate(inter)
      allocate(inter%interm(20))
      inter%size=0
      j=1
      



      ! Loop over formula_items, end of the list points to NULL
      do while (associated(fl_item%next))

      
      if (associated(fl_item%interm)) then
!         ! if fl_item%interm == true, then the next item in the list
!         ! will be the intermediate contraction; so move to the next
!         ! item in the list
!         if (.not.associated(fl_item%next)) exit
!         fl_item=>fl_item%next
!
!         allocate(inter%interm(j)%spin_case(6))
!         inter%size = j
!         j = j + 1
!         ! Set up index of intermediate
!         ! Spin summ for all cases and store results
!         ! continue until we hit the next result
!         ! Then assign index and spin sume result
!         ! See which interemdiates are needed
!         ! Then print the interediates, then the result
!         ! dellocate inter%interm array
!         call intermediate_to_itf(fl_item%bcontr,itflog,fl_item%command,
!     &                            inter%interm(j))


         ! New idea, recursive search back along the list
         ! Mark point where intermediates start
         if (.not.associated(fl_item%next)) exit
         inter_start => fl_item%next

         ! Cycle through list until we hit the next result which isn't
         ! an intermediate. All previous intermediates are assumed to
         ! contribute to the following result lines
         do
            fl_item => fl_item%next

            if(associated(fl_item%interm)) cycle
            if(.not.associated(fl_item%next)) exit
            if(scan(fl_item%bcontr%label_res, "STIN")==0) then
               res_start => fl_item%next
               exit
            end if

         end do

         ! We want to build an array of intermediate names and their
         ! various spin cases here, for now each line only depends on 1
         ! intermediate, so spin_inters(1) is all we need.
         ! TODO: fix this
         call intermediate_to_itf(fl_item%bcontr,itflog,fl_item%command,
     &                            spin_inters)
         write(itflog,*) "Hi: ", spin_inters(1)%cases
         write(itflog,*) "Hi2: ", spin_inters(1)%name

         ! Go back to inter_start and look for the intermediates
         fl_item => inter_start
         
         do
            if (fl_item%bcontr%label_res == spin_inters(1)%name) then
               !call check_inter(fl_item%bcontr%label_res,logic) 
               exit
            end if
            fl_item => fl_item%next
         end do

         ! Need to check if they depend on any other intermediates
         ! subroutine to check
         ! call intermediate_to_itf() to update spin_inters
         ! back to do loop and repeat with new name

         ! Once we have the complete spin_inters array - this contains
         ! all the info for all the intermediates used in the future
         ! residual. So now we can spin summ these and print them out, +
         ! change there names STIN001_aaaa
         ! Then print out the residual

         ! Check if next residual needs intermdiates and which spin
         ! cases are needed.
         ! Exit if next resdiual is different or a new intermediate is
         ! declared


      else
         ! Not an intermediate, so select the correct command case

         select case(fl_item%command)
         case(command_end_of_formula)
!           write(itflog,*) '[END]'
         case(command_set_target_init)
!           write(itflog,*) '[INIT TARGET]',fl_item%target
         case(command_set_target_update)
!           write(itflog,*) '[SET TARGET]',fl_item%target
         case(command_new_intermediate)
!           write(itflog,*) '[NEW INTERMEDIATE]',fl_item%target
!           write(itflog,'(2x,a)') trim(fl_item%interm%name)
!           write(itflog,'(2x,"attribute parentage: ",a," ",a)')
!        &                        trim(fl_item%parent1),
!        &                        trim(fl_item%parent2)
!           write(itflog,'(2x,"incore: ",i2)') fl_item%incore
!           call print_op_occ(itflog,fl_item%interm)
         case(command_del_intermediate)
!           write(itflog,*) '[DELETE INTERMEDIATE]',fl_item%target
!           write(itflog,'(2x,a)') trim(fl_item%label)
         case(command_add_contribution)
           write(itflog,*) '[CONTR]',fl_item%target
         case(command_add_intm)
!           write(itflog,*) '[ADD]',
!        &       fl_item%target
!           call prt_bcontr(itflog,fl_item%bcontr)

           call command_to_itf(fl_item%bcontr,itflog,fl_item%command)

         case(command_cp_intm)
!           write(itflog,*) '[COPY]',
!        &       fl_item%target
!           call prt_bcontr(itflog,fl_item%bcontr)

           ! Assume that copy means alloc, then :=
           call command_to_itf(fl_item%bcontr,itflog,fl_item%command)

         case(command_add_bc)
!           write(itflog,*) '[CONTRACT][ADD]',
!        &       fl_item%target
!           call prt_bcontr(itflog,fl_item%bcontr)

           call command_to_itf(fl_item%bcontr,itflog,fl_item%command)

         case(command_add_bc_reo)
           write(itflog,*) '[CONTRACT][REORDER][ADD]',
     &          fl_item%target
           call prt_bcontr(itflog,fl_item%bcontr)
           call prt_reorder(itflog,fl_item%reo)
         case(command_bc)
!           write(itflog,*) '[CONTRACT]',
!        &       fl_item%target
!           call prt_bcontr(itflog,fl_item%bcontr)

           call command_to_itf(fl_item%bcontr,itflog,fl_item%command)

         case(command_bc_reo)
!           write(itflog,*) '[CONTRACT][REORDER]',
!        &       fl_item%target
!           call prt_bcontr(itflog,fl_item%bcontr)
!           call prt_reorder(itflog,fl_item%reo)

           call command_to_itf(fl_item%bcontr,itflog,fl_item%command)

         case(command_reorder)
!           write(itflog,*) '[REORDER]',fl_item%target
!           call prt_reorder(itflog,fl_item%reo)
         case(command_add_reo)
           write(itflog,*) '[REORDER][ADD]',
     &          fl_item%target
           call prt_bcontr(itflog,fl_item%bcontr)
           call prt_reorder(itflog,fl_item%reo)
         case(command_symmetrise)
           write(itflog,*) '[SYMMETRISE]',fl_item%target
         case default
           write(itflog,*) 'unknown command ',fl_item%command,
     &          fl_item%target
         end select

      end if

      ! Optionally print the formula items to another output file
      if (print_form) then
        call print_form_item2(formlog,'LONG',contr_no,fl_item,op_info)
      end if

      ! Check if at the end of the list, if not, point to the next item
      !if (.not.associated(fl_item%next)) exit
      fl_item=>fl_item%next

      ! Count the number of terms
      contr_no=contr_no+1

      ! Deallocated spin cases of an intermediate
      if (inter%size /= 0 .and. .not.associated(fl_item%interm)) then
         do i = 1, inter%size
            deallocate(inter%interm(i)%spin_case)
         end do
         j = 1
         inter%size = 0
      end if

      end do


      deallocate(inter%interm)
      deallocate(inter)

      end
