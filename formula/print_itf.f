*----------------------------------------------------------------------*
      pure function check_inter(label)
*----------------------------------------------------------------------*
!    Check if tensor is an intermediate
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'
      
      character(len=maxlen_bc_label), intent(in) ::
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
      subroutine print_itf(itflog,fl_head,op_info,print_form,formlog)
*----------------------------------------------------------------------*
*     Print ITF info to itflog
*----------------------------------------------------------------------*
      !use itf_utils copied above from module
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
     &     inter_start,    ! Mark start of intermediate search
     &     res_start,      ! Mark end of intermediate search
     &     summed_inter    ! Mark intermediate points in intermediate search
      type(spin_cases), dimension(4) ::
     &     spin_inters  ! Array of intermeidates with associated spin cases
      integer ::
     &     i,j,k,     ! Loop indcies
     &     contr_no  ! Counter of contrations
      logical ::
     &     check_inter,    ! Need to use itf module instead
     &     more_inter,     ! Check if more intermediates are needed
     &     finished_inter  ! Check if finished recursive intermediate search
      integer ::
     &     tmp_case(4),
     &     ninter       ! Number of intermediates found in recursive search

      ! Point to start of linked list
      fl_item => fl_head
      contr_no = 0


      ! Loop over formula_items, end of the list points to NULL
      do while (associated(fl_item%next))
      
      ! Check if formula item is an intermediate
      if (associated(fl_item%interm)) then

         ! Recursive search back along the list.
         ! Mark point where intermediates start
!         write(itflog,*) "Starting intermediate search"
         if (.not.associated(fl_item%next)) exit
         inter_start => fl_item%next

         ! Cycle through list until we hit the next result which isn't
         ! an intermediate. All previous intermediates are assumed to
         ! contribute to the following result lines
         do
            fl_item => fl_item%next

            if(associated(fl_item%interm)) cycle
            if(.not.associated(fl_item%next)) then
!               write(itflog,*) "ERROR: intermediate was declared, but
!     &                          not used!"
               exit
            end if
            if(scan(fl_item%bcontr%label_res, "STIN")==0) then
!               write(itflog,*) "Found next result"
               res_start => fl_item
               exit
            end if

         end do

!         write(itflog,*) "res_start ", res_start%bcontr%label_res

         ! We want to build an array of intermediate names and their
         ! various spin cases here
         ! Global variable ninter is the number of intermediates which
         ! we need to deal with
         ninter = 0
         call intermediate_to_itf(fl_item%bcontr,itflog,fl_item%command,
     &                            spin_inters, ninter)

!         do i = 1, ninter
!            write(itflog,*) "spin_cases: ", spin_inters(ninter)%cases
!            write(itflog,*) "name: ", spin_inters(ninter)%name
!         end do

         ! Go back to inter_start and look for the intermediates
         fl_item => inter_start
         ! Marker needs to be set to where we have already check for
         ! intermediates
         summed_inter => res_start
         finished_inter = .false.
         

         ! Recursive search through list to get all information about
         ! every intermediate used to produce a result
         do while (.not.finished_inter)
            ! Check if infomation about the intermediate is needed, all required
            ! intermediates are stored in spin_inters
            do i = 1, ninter
               if (fl_item%bcontr%label_res == spin_inters(i)%name)
     &                      more_inter = .true.
            end do

            !if (fl_item%bcontr%label_res == spin_inters(1)%name) then
            ! If the intermediate is needed, check if it depends on any
            ! intermediates + find out their names/ spin cases
            if (more_inter) then
               if (check_inter(fl_item%bcontr%label_op1) .or.
     &             check_inter(fl_item%bcontr%label_op2)) then
                  call intermediate_to_itf(fl_item%bcontr,itflog,
     &                               fl_item%command,spin_inters,ninter)

                  ! Mark the position of the previously checked
                  ! intermediates
                  summed_inter => fl_item

                  ! Move back to the start and start search for next
                  ! intermediate
                  fl_item => inter_start
               else
                  ! If it doesn't have any intermediates, move onto the
                  ! next item
                  fl_item => fl_item%next

                  ! If the next item is has reached the begining of the
                  ! previously summed intermediate, then we need to
                  ! break out this loop.
                  if (associated(fl_item,summed_inter)) finished_inter =
     &                                                     .true.
               end if
            else
               ! If the interemediate isn't needed, move onto the next
               ! item
               fl_item => fl_item%next
               if (associated(fl_item,summed_inter)) finished_inter =
     &                                                  .true.
            end if

            more_inter = .false.

            ! Check we haven't reached the residual result
            if (associated(fl_item,res_start)) then
!               write(itflog,*) "Found the end"
               fl_item => res_start
               finished_inter = .true.
            end if
         end do

         ! Once we have the complete spin_inters array - this contains
         ! all the info for all the intermediates used in the future
         ! residual. So now we can spin summ these and print them out, +
         ! change there names STIN001aaaa, then print out the residual
         fl_item => inter_start

         ! Loop over intermediates.
         ! We want to loop over all lines of intermediate, before doing
         ! another/ printing the next spin case. I_aaaa, then I_abab
         do k = 1, ninter
            ! Loop over the number of spin cases for an intermediate
            ! For each spin case of an intermediate, we want to loop though
            ! the list and print out all the lines which contribute
            do i = 1, spin_inters(k)%ncase
               do ! Loop through the list
                  if (fl_item%bcontr%label_res == spin_inters(k)%name)
     &            then
                     ! Send off specific spin case to be summed and printed
                     do j = 1, 4
                        tmp_case(j) = spin_inters(k)%cases(j,i)
                     end do
                     call intermediate2_to_itf(fl_item%bcontr,itflog,
     &                                      fl_item%command,tmp_case)
                  end if

                  ! Move onto next item and repeat
                  fl_item => fl_item%next 

                  if (associated(fl_item,res_start)) then
!                     write(itflog,*) "Finished one spin inter block"
                     ! Go back to start and print out remaing spin cases +
                     ! reapeat
                     fl_item => inter_start
                     exit
                  end if
               end do
            end do
         end do

         ! Spin summ and print residual which uses the above
         ! intermediates
         fl_item => res_start
         call command_to_itf(fl_item%bcontr,itflog,fl_item%command)


         ! Not needed for now, but maybe in the future:
         ! Check if next residual needs intermdiates and which spin
         ! cases are needed.
         ! Exit if next resdiual is different or a new intermediate is
         ! declared. 

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
      fl_item => fl_item%next

      ! Count the number of terms
      contr_no = contr_no+1

      end do

      return
      end
