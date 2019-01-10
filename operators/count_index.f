*----------------------------------------------------------------------*
      module itf_utils
*----------------------------------------------------------------------*
!     Contains functions used throughout the code
*----------------------------------------------------------------------*
         implicit none
      contains
*----------------------------------------------------------------------*
      function trimal(str1) result(str2)
*----------------------------------------------------------------------*
!     trim() and adjustl() a string
*----------------------------------------------------------------------*
      implicit none

      character(len=*), intent(in) ::
     &     str1
      character(len=:), allocatable ::
     &     str2

      str2 = trim(adjustl(str1))

      end function trimal

*----------------------------------------------------------------------*
      pure function rename_tensor(string, rank)
*----------------------------------------------------------------------*
!     Rename tensor acording to taste
!     This should be expanded to rename all tensors so they correspond
!     with the ITF algo file
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'

      character(len=maxlen_bc_label), intent(in) ::
     &    string
      integer, intent(in) ::
     &    rank          ! Rank of tensor
      character(len=maxlen_bc_label) ::
     &    rename_tensor

      if (trim(string).eq.'O2g' .or. trim(string).eq.'O2') then
          rename_tensor='R'
      else if (trim(string).eq.'O1') then
          rename_tensor='R'
      else if (trim(string).eq.'T2g' .or. trim(string).eq.'T2') then
          rename_tensor='T'
      else if (trim(string).eq.'T1') then
          rename_tensor='T'
      else if (trim(string).eq.'H') then
          if (rank==2) then
             rename_tensor='f'
          else
             rename_tensor='K'
          end if
      else if (trim(string).eq.'GAM0') then
          if (rank==2) then
             rename_tensor='Dm1'
          else if (rank==4) then
             rename_tensor='Dm2'
          else if (rank==6) then
             rename_tensor='Dm3'
          !else if (rank==0)
          !   rename_tensor=''
          else
             rename_tensor='Dm'
          end if
      else if (trim(string).eq.'ECCD' .or. trim(string).eq.'MRCC_LAG')
     & then
          rename_tensor='ECC'
      else
          rename_tensor=trim(string)
      end if

      if (rename_tensor(1:1)=='_') rename_tensor(1:1)=''

      end function

*----------------------------------------------------------------------*
      pure function check_inter(label)
*----------------------------------------------------------------------*
!     Check if tensor is an intermediate
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
     &    index(label, "LTIN")>0 .or.
     &    index(label, "PTIN")>0) then
         check_inter=.true.
      else
         check_inter=.false.
      end if

      end function


*----------------------------------------------------------------------*
      pure function check_int(label)
*----------------------------------------------------------------------*
!     Check if tensor is an integral
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      !include 'def_itf_contr.h'

      character(len=maxlen_bc_label), intent(in) ::
     &     label

      logical ::
     &     check_int

      ! Assume these are the names of intermediates
      if (index(label, "H")>0) then
         check_int=.true.
      else
         check_int=.false.
      end if

      end function


*----------------------------------------------------------------------*
      pure function t_index(index, upper)
*----------------------------------------------------------------------*
!     Transpose ITF index string, abij => abji
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=index_len), intent(in) ::
     &     index       ! ITF index string
      logical, optional, intent(in) ::
     &     upper
      character(len=index_len) ::
     &     t_index      ! Transpose of ITF index string

      logical ::
     &     contra

      if (present(upper)) then
         contra = upper
      else
         contra = .false.
      end if

      t_index=index

      if (contra) then
         t_index(3:3)=index(4:4)
         t_index(4:4)=index(3:3)
      else
         t_index(1:1)=index(2:2)
         t_index(2:2)=index(1:1)
      end if

      end function


      end module itf_utils

*----------------------------------------------------------------------*
      subroutine line_error(error,item)
*----------------------------------------------------------------------*
!     Print error message in the itf logfile
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=*), intent(in) ::
     &     error
      type(itf_contr), intent(in) ::
     &     item

      write(item%logfile,*) "ERROR: ", error
      write(item%logfile,*) "================================="
      write(item%logfile,*) "Result: ", item%label_res, item%idx3
      write(item%logfile,*) "Tensor1: ", item%label_t1, item%idx1
      write(item%logfile,*) "Tensor2: ", item%label_t2, item%idx2
      write(item%logfile,*)

      return
      end

*----------------------------------------------------------------------*
      subroutine count_index(idx,iocc,irstr,ngas,nspin,nops)
*----------------------------------------------------------------------*
!     Count index of operator
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     iocc(ngastp,2), ngas, nspin,
     &     irstr(2,ngas,2,2,nspin), idx
      integer, intent(inout) ::
     &     nops(4,2)                     ! Matrix of index info

      ! Creation operators
      ! Particle, ijkl
      nops(1,1)=nops(1,1) + iocc(1,1)
      ! Hole, abcd
      nops(2,1)=nops(2,1) + iocc(2,1)
      ! Valence, pqrs
      nops(3,1)=nops(3,1) + iocc(3,1)
      ! Explicit, x
      nops(4,1)=nops(4,1) + iocc(4,1)

      ! Annhilation operators as above
      nops(1,2)=nops(1,2) + iocc(1,2)
      nops(2,2)=nops(2,2) + iocc(2,2)
      nops(3,2)=nops(3,2) + iocc(3,2)
      nops(4,2)=nops(4,2) + iocc(4,2)

      return
      end


*----------------------------------------------------------------------*
      subroutine command_to_itf(contr_info,itflog,command)
*----------------------------------------------------------------------*
!     Take GeCco binary contraction and produce ITF algo code.
!     Includes antisymmetry of residual equations and spin summation.
*----------------------------------------------------------------------*

      use itf_utils
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h' ! For def_formular_item.h
      include 'def_contraction.h'
      include 'def_formula_item.h' ! For command parameters
      include 'def_itf_contr.h'

      type(binary_contr), intent(inout) ::
     &     contr_info      ! Inofrmation about binary contraction
      integer, intent(in) ::
     &     itflog,         ! Output file
     &     command         ! Type of formula item command, ie. contraction, copy etc.

      type(itf_contr) ::
     &     itf_item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
!     &    perm_array(4),   ! Info of permutation factors
     &    perm_case,   ! Info of permutation factors
     &    i                ! Loop index
      logical ::
     &    inter            ! True if result is an intermediate
      character(len=maxlen_bc_label) ::
     &    old_name

      ! Initalise permutation factors to 0 == no permutation
      perm_case = 0

      ! Determine if result needs permuting
      inter = check_inter(contr_info%label_res)

      if (.not.inter) then
         call permute_tensors(contr_info,perm_case,itflog)
      end if

      ! If the perm_array doesn't contain any zeros, then we should
      ! introduce an interemediate which collects half of the different
      ! permutation cases, then do:
      ! .R[abij] += I[abij]
      ! .R[abij] += I[baji]
      ! Save old name and replace it with a new one
      if (perm_case > 0) then
         old_name = contr_info%label_res
         contr_info%label_res = "ITIN"
         itf_item%symm = .true.
      end if


      ! Pick out specific commands, form the itf_contr object, spin sum
      ! and print out contraction line
      if (command==command_add_intm .or. command==command_cp_intm) then
         ! For [ADD] and [COPY] cases
         call itf_contr_init(contr_info,itf_item,0,command,itflog)
         call print_itf_line(itf_item,.false.,.false.)
      else
         ! For other binary contractions
         if (perm_case == 0) then
            ! No permutations
            call itf_contr_init(contr_info,itf_item,0,command,itflog)
            call assign_spin(itf_item)
         else
            do i=1, perm_case
               ! Loop over permuation cases and send seperatley to
               ! assign_spin. For most cases this is just one, however
               ! for (1-Pij)(1-Pab), we need to generate one of these
               ! permutaions before symmetrising
               call itf_contr_init(contr_info,itf_item,i,command,itflog)

               ! Mutliply factor by -1.0 due to permutation
               if (i == 2) then
                  itf_item%fact = itf_item%fact * -1.0

                  ! Need to transpose by tensors after permutation, to
                  ! avoid symmetry problem when using (1 + Pabij)
                  itf_item%idx1 = t_index(itf_item%idx1)
                  itf_item%idx2 = t_index(itf_item%idx2)
               end if

               call assign_spin(itf_item)
            end do

            ! If created a perm intermedite, print the symmetrised lines
            call print_symmetrise(old_name,itf_item)
         end if
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine intermediate_to_itf(contr_info,itflog,command,
     &                               spin_inters,n_inter,permute)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*
      use itf_utils

      implicit none
      include 'opdim.h'
      include 'mdef_operator_info.h' ! For def_formular_item.h
      include 'def_contraction.h'
      include 'def_formula_item.h' ! For command parameters
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &     contr_info      ! Inofrmation about binary contraction
      integer, intent(in) ::
     &     itflog,         ! Output file
     &     command         ! Type of formula item command, ie. contraction, copy etc.
      integer, intent(inout) ::
     &     n_inter         ! Overall number of intermediates needed by result
      type(spin_cases), dimension(MAXINT), intent(inout) ::
     &     spin_inters
      integer, intent(inout) ::
     &     permute         ! 2 = Need to permute index

      type(itf_contr) ::
     &     itf_item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &    perm_array(4),   ! Info of permutation factors
     &    tmp_case(4),
     &    i,j,k
      logical ::
     &    summed

      perm_array=0

!      call itf_contr_init(contr_info,itf_item,perm_array(1),
!     &                    command,itflog)
      call itf_contr_init(contr_info,itf_item,permute,
     &                    command,itflog)

      !! Multiply factor by -1.0 due to permutation
      if (permute == 2) then
         itf_item%fact = itf_item%fact * -1.0

         ! Need to transpose by tensors after permutation, to
         ! avoid symmetry problem when using (1 + Pabij)
         itf_item%idx1 = t_index(itf_item%idx1)
         itf_item%idx2 = t_index(itf_item%idx2)
      end if

      ! Allocate space to store information about intermediates and
      ! their spin cases. Only allocate 2 objects as there can only be
      ! at most two intermediates on a line
      allocate(itf_item%inter_spins(2))

      ! Do not want to print out the lines while gathering info about
      ! intermediates
      itf_item%print_line = .false.


      ! Need to catch lines which don't need to be spin summed
      call assign_simple_spin(itf_item, summed)

      ! Need to set spin result for intermediate that depends on
      ! another intermediate
      if (.not. summed) then

         if (itf_item%inter(3)) then
            do i = 1, MAXINT
               if (itf_item%label_res == spin_inters(i)%name) then
                  do k = 1, spin_inters(i)%ncase
                     do j = 1, 4
                        tmp_case(j) = spin_inters(i)%cases(j,k)
                     end do
                     itf_item%spin_case = tmp_case

                     call assign_spin(itf_item)
                  end do
               end if
            end do
         else
            ! Result is not an intermediate, ie. its a residual, so lets find
            ! out what spin intermediates it needs.
            call assign_spin(itf_item)
         end if
      end if

      itf_item%print_line = .true.

      ! TODO: Ignore rank 6 cases and tensor product cases for now...
      if (itf_item%rank3 /= 4 .and. itf_item%rank1 /= 2 .and.
     &    itf_item%rank2 /= 2) then
      if (itf_item%rank3 == 6 .or. itf_item%rank1 == 6 .or.
     &    itf_item%rank2 == 6) then
      if (itf_item%ninter == 0) call line_error("Couldn't find
     & intermediate", itf_item)
      end if
      end if


      ! Copy information back to array in print_itf()
      do i = 1, itf_item%ninter
         spin_inters(i+n_inter)=itf_item%inter_spins(i)
         !write(itf_item%logfile,*) "INER SPIN: ",itf_item%inter_spins(i)
      end do

      ! Overall number of intermediates used to index spin_inters
      n_inter = n_inter + itf_item%ninter

      deallocate(itf_item%inter_spins)

      return
      end


*----------------------------------------------------------------------*
      subroutine find_spin_intermediate(contr_info,itflog,command,
     &                                  spin_inters,n_inter)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*
      use itf_utils

      implicit none
      include 'opdim.h'
      include 'mdef_operator_info.h' ! For def_formular_item.h
      include 'def_contraction.h'
      include 'def_formula_item.h' ! For command parameters
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &     contr_info      ! Inofrmation about binary contraction
      integer, intent(in) ::
     &     itflog,         ! Output file
     &     command         ! Type of formula item command, ie. contraction, copy etc.
      integer, intent(inout) ::
     &     n_inter         ! Overall number of intermediates needed by result
      type(spin_cases), dimension(MAXINT), intent(inout) ::
     &     spin_inters

      integer ::
     &    perm_case,   ! Info of permutation factors
     &    i                ! Loop index
      logical ::
     &    inter            ! True if result is an intermediate
      character(len=maxlen_bc_label) ::
     &    old_name


      ! Initalise permutation factors to 0 == no permutation
      perm_case = 0

      ! Determine if result needs permuting
      inter = check_inter(contr_info%label_res)

      if (.not.inter) then
         call permute_tensors(contr_info,perm_case,itflog)
      end if

      if (perm_case == 0) then
         ! No permutations
         call intermediate_to_itf(contr_info,itflog,command,
     &                            spin_inters,n_inter,1)
      else
         do i=1, perm_case
            call intermediate_to_itf(contr_info,itflog,command,
     &                               spin_inters,n_inter,i)
         end do
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine assign_simple_spin(item, summed)
*----------------------------------------------------------------------*
!     Assign spin to simple cases using logic conditions
!     This avoid the spin sum routine
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &    item
      logical, intent(inout) ::
     &    summed

      integer ::
     &    i,j


      summed = .false.

!      if (item%rank3 == 2 .and. item%rank1 + item%rank2 == 2 .or.
      if (item%rank3 + item%rank1 + item%rank2 == 4 .or.
     &    item%rank3 + item%rank1 + item%rank2 == 6) then
         j = 1
         do i = 1, 2
            if (item%inter(i)) then
               ! Only need the aa case
               item%inter_spins(j)%cases(:,j) = (/ 1, 0, 1, 0 /)
               item%inter_spins(j)%ncase = 1
               if (i == 1) item%inter_spins(j)%name = item%label_t1
               if (i == 2) item%inter_spins(j)%name = item%label_t2
               item%ninter = item%ninter + 1
               j = j + 1
            end if
         end do

         !write(item%logfile,*) "SIMPLE SPIN: ", item%inter_spins
         summed = .true.
      else if (item%rank3 + item%rank1 + item%rank2 == 0) then
         ! Do nothing for scalar tensors
         j = 1
         do i = 1, 2
            if (item%inter(i)) then
               ! Only need the aa case
               item%inter_spins(j)%cases(:,j) = (/ 0, 0, 0, 0 /)
               item%inter_spins(j)%ncase = 1
               if (i == 1) item%inter_spins(j)%name = item%label_t1
               if (i == 2) item%inter_spins(j)%name = item%label_t2
               item%ninter = item%ninter + 1
               j = j + 1
            end if
         end do

         !write(item%logfile,*) "SIMPLE SPIN: ", item%inter_spins
         summed = .true.
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine intermediate2_to_itf(contr_info,itflog,command,
     &                               label,spin_case)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'mdef_operator_info.h' ! For def_formular_item.h
      include 'def_contraction.h'
      include 'def_formula_item.h' ! For command parameters
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &     contr_info      ! Inofrmation about binary contraction
      integer, intent(in) ::
     &     itflog,         ! Output file
     &     command         ! Type of formula item command, ie. contraction, copy etc.
      integer, intent(in) ::
     &     spin_case(4)
      character(len=maxlen_bc_label), intent(in) ::
     &     label

      type(itf_contr) ::
     &     itf_item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &    perm_case,
     &    i,j
      character(len=4) ::
     &    spin_name


      perm_case = 0

      call itf_contr_init(contr_info,itf_item,perm_case,
     &                    command,itflog)

      ! Set overall spin case of result
      itf_item%spin_case = spin_case

      ! Change intermediate name to reflect spin case
      j = 1
      spin_name = ''

      do i = 1, size(spin_case)
         if (spin_case(i)==1) then
            spin_name(j:j) = 'a'
            j = j + 1
         else if (spin_case(i)==2) then
            spin_name(j:j) = 'b'
            j = j + 1
         end if
      end do

!      write(itf_item%logfile,*) "NAME: ", spin_name
!      write(itf_item%logfile,*) "NAME: ", spin_case
!      write(itf_item%logfile,*) "NAME: ", itf_item%idx3
!      write(itf_item%logfile,*) "NAME: ", itf_item%idx1
!      write(itf_item%logfile,*) "NAME: ", itf_item%idx2

      ! If an intermediate arises as the result of a permutation, we
      ! need to create this new intermeidate. This requires the
      ! transpose
      if (scan('P', label)) then
         itf_item%idx3 = t_index(itf_item%idx3)
         itf_item%idx2 = t_index(itf_item%idx2)
         itf_item%idx1 = t_index(itf_item%idx1,.true.)
      end if

      itf_item%label_res = trim(itf_item%label_res)//trim(spin_name)
      call assign_spin(itf_item)

      return
      end

*----------------------------------------------------------------------*
      subroutine print_itf_line(item,s1,s2)
*----------------------------------------------------------------------*
!     Print line of ITF code
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'mdef_operator_info.h' ! For def_formular_item.h
      include 'def_contraction.h'
      include 'def_formula_item.h' ! For command parameters
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &     item
      logical, intent(in) ::
     &     s1,s2

      character(len=maxlen_bc_label) ::
     &     nres, nt1, nt2          ! Name of tensors involved in the contraction
      character(len=5) ::
     &     s_int                 ! Intermdiate tensor number
      character(len=70) ::
     &     itf_line,          ! Line of ITF code
     &     st1, st2           ! Name of spin summed tensors + index
      character(len=2) ::
     &     equal_op           ! ITF contraction operator; ie. +=, -=, :=
      character(len=25) ::
     &     sfact,             ! String representation of factor
     &     sfact_star         ! String representation of factor formatted for output
      integer ::
     &     i

      ! Change names of specific tensors
      nres=rename_tensor(item%label_res, item%rank3)
      nt1=rename_tensor(item%label_t1, item%rank1)
      nt2=rename_tensor(item%label_t2, item%rank2)

      ! When permuting, the intermediate name will flip spin. We need it
      ! to be the same as the previously declared intermediates, so this
      ! flips it back.
      ! TODO: Is this needed anymore, new permutation intermediate uses
      ! its own spin label
      !if (item%permute>1) then
      !   if (item%inter(1)) then
      !      do i = 1, item%rank1
      !         if (item%inter1(i:i)=='a') then
      !            item%inter1(i:i)='b'
      !         else if (item%inter1(i:i)=='b') then
      !            item%inter1(i:i)='a'
      !         end if
      !      end do
      !   end if

      !   if (item%inter(2)) then
      !      do i = 1, item%rank2
      !         if (item%inter2(i:i)=='a') then
      !            item%inter2(i:i)='b'
      !         else if (item%inter2(i:i)=='b') then
      !            item%inter2(i:i)='a'
      !         end if
      !      end do
      !   end if
      !end if

      if (item%permute>1) then
         if (item%inter(1) .or. item%inter(2)) then
            ! Because we take the transpose of the spin orbial eqns twice,
            ! we retain the negative sign from the permutation. This is
            ! only a problem when we define a new spin intermediate
            ! which results from the permutation of a result line
            item%fact = item%fact * -1.0
         end if
      end if

      ! Add intermediate spin strings to names
      if (item%inter(1)) nt1 = trim(nt1)//trim(item%inter1)
      if (item%inter(2)) nt2 = trim(nt2)//trim(item%inter2)

      ! Change tensor to spatial orbital quantity, unless it is an
      ! intermediate
      if (s1 .and. .not.item%inter(1)) then
         ! Pure spin
         st1='('//trimal(nt1)//'['//trim(item%idx1)//']'//
     &       ' - '//trimal(nt1)//'['//trim(t_index(item%idx1))//']'//')'
      else
         st1=trimal(nt1)//'['//trim(item%idx1)//']'
      end if

      if (s2 .and. .not.item%inter(2)) then
         ! Pure spin
         st2='('//trimal(nt2)//'['//trim(item%idx2)//']'//
     &       ' - '//trimal(nt2)//'['//trim(t_index(item%idx2))//']'//')'
      else
         if (item%command==command_add_intm .or.
     &       item%command==command_cp_intm) then
            ! Don't need second operator for [ADD] or [COPY]
            st2=''
         else
            st2=trimal(nt2)//'['//trim(item%idx2)//']'
         end if
      end if

      ! Convert factor to string, ignore if 1.0 or -1.0
      sfact=''
      sfact_star=''
      if (item%fact.ne.1.0) then
          if (item%fact.ne.-1.0) then
              write(sfact,*) item%fact
              do i=1, len(sfact)
                 ! Start after decimal point
                 if (sfact(i:i)=='0' .or. sfact(i:i)=='-') then
                    sfact(i:i)=' '
                 end if
              end do
              sfact_star=trimal(sfact)//'*'
          end if
      end if

      ! Deterime what the contraction operator looks like
      equal_op='  '
      if (item%fact.lt.0.0) then
         equal_op='-='
      else if (item%command==command_cp_intm) then
         equal_op=':='
      else
         equal_op='+='
      end if

      ! Construct complete itf algo line from the above parts
      itf_line='.'//trimal(nres)//
     &    '['//trim(item%idx3)//'] '//equal_op//' '//
     &    trim(sfact_star)//trimal(st1)//' '//trimal(st2)

      ! Print it to bcontr.tmp
      write(item%logfile,*) trim(itf_line)

      return
      end

*----------------------------------------------------------------------*
      subroutine assign_add_index(contr_info,item)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &     contr_info   ! Information about binary contraction
      type(itf_contr), intent(inout) ::
     &     item
      integer ::
     &     o1(4,2)      ! Occupations of first tensor
      integer ::
     &     i,j
      character, dimension(4) ::
     &     hol=(/ 'i','j','k','l' /),
     &     par=(/ 'a','b','c','d' /)
      character, dimension(8) ::
     &     val=(/ 'p','q','r','s','t','u','v','w' /)
      character(len=8) ::
     &     c1='        ',    ! Workspace to assign index before array
     &     c2='        ',
     &     c3='        ',
     &     c4='        ',
     &     a1='        ',
     &     a2='        ',
     &     a3='        ',
     &     a4='        '
      character(len=8), dimension(8) ::
     &     o1_array     ! Index of operator 1
      integer, dimension(3) ::
     &     idx_type     ! Info about index convention
      integer, parameter ::
     &     t_amp = 0,       ! [apij] (aacc)
     &     ham   = 1        ! [abip]
      character(len=20), dimension(4) ::
     &     tensor_ham=(/ 'H', 'INT_D', 'INT_HH', 'INT_PP' /)  ! Tensor names to use ham index convention

      ! Set index convention
      idx_type=(/ 0, 0, 0 /)
      do i=1, len(tensor_ham)
          if (contr_info%label_op1.eq.trim(tensor_ham(i))) then
              ! Use default convention for now
              idx_type(1)=0
          end if
          if (contr_info%label_res.eq.trim(tensor_ham(i))) then
              idx_type(3)=0
          end if
      end do

      o1=0

      ! Get occuation info
      do i = 1, contr_info%nj_op1
        call count_index(i,
     &     contr_info%occ_op1(1:,1:,i),
     &     contr_info%rst_op1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,o1)
      end do

      ! Order in ITF usually follows: apij
      ! Defualt [ccaa] as in the case of T[abij]

      ! Assign e1 (external indicies of t1)
      do i=1, o1(2,1)
          c1(i:)=par(i)
      end do
      o1_array(1)=c1
      do i=1, o1(3,1)
          c2(i:)=val(i)
      end do
      o1_array(2)=c2
      do i=1, o1(1,1)
          c3(i:)=hol(i)
      end do
      o1_array(3)=c3

      ! Need to to be shifted to not match assignment of creations above
      do i=1, o1(2,2)
          a1(i:)=par(i+o1(2,1))
      end do
      o1_array(5)=a1
      do i=1, o1(3,2)
          a2(i:)=val(i+o1(3,1))
      end do
      o1_array(6)=a2
      do i=1, o1(1,2)
          a3(i:)=hol(i+o1(1,1))
      end do
      o1_array(7)=a3

      c1='        '
      c2='        '
      c3='        '
      a1='        '
      a2='        '
      a3='        '

      ! Construct final index strings
      ! Operator 1
      select case(idx_type(1))
      case(ham)
      ! Hamiltonian/integral convention
      item%idx1=trimal(o1_array(1))//trimal(o1_array(2))//
     &          trimal(o1_array(3))//trimal(o1_array(5))//
     &          trimal(o1_array(6))//trimal(o1_array(7))
      case default
      ! [apij] (aacc), ie. T[abij]
      item%idx1=trimal(o1_array(1))//trimal(o1_array(2))//
     &          trimal(o1_array(3))//trimal(o1_array(5))//
     &          trimal(o1_array(6))//trimal(o1_array(7))
      end select

      item%idx3=item%idx1

      return
      end


*----------------------------------------------------------------------*
      subroutine construct_slot(occ_array,slot)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &     occ_array(4,2)
      type(tensor_slot), intent(inout) ::
     &     slot

      integer :: i,j

      j = 0
      do i = 1, 4
         j = j + occ_array(i,1)
      end do

      allocate(slot%cre(j))
      slot%cslots = j

      do i = 1, j
         slot%cre(i) = ' '
      end do

      j = 0
      do i = 1, 4
         j = j +occ_array(i,2)
      end do

      allocate(slot%ann(j))
      slot%aslots = j

      do i = 1, j
         slot%ann(i) = ' '
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine desctruct_slot(slot)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(tensor_slot), intent(inout) ::
     &     slot

      deallocate(slot%ann)
      deallocate(slot%cre)

      return
      end


*----------------------------------------------------------------------*
      subroutine assign_index(contr_info,item)
*----------------------------------------------------------------------*
!     Assign an letter index to each tensor in a line
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &     contr_info   ! Information about binary contraction
      type(itf_contr), intent(inout) ::
     &     item
      integer ::
     &     c(4,2),       ! Occupations of contraction index
     &     e1(4,2),      ! Occupations of external index 1
     &     e2(4,2),     ! Occupations of external index 2
     &     e3(4,2),     ! Occupations of external index 2
     &     nc(4,2),       ! Copies of above but in different order
     &     ne1(4,2),      ! Copies
     &     ne2(4,2),    ! Copies
     &     nr(4,2)    ! Copies of result ops
      integer ::
     &     i,j,k,l,m,n,ii
      character, dimension(3) ::
     &    chol=(/ 'k','l','m' /)
      character, dimension(4) ::
     &     hol=(/ 'i','j','k','l' /),
     &     par=(/ 'a','b','c','d' /),
     &     cpar=(/ 'c','d','e','f' /),
     &     val2=(/ 'p','q','r','s' /),
     &     cval2=(/ 'r','s','t','u' /)
      character, dimension(8) ::
     &     val=(/ 'p','q','r','s','t','u','v','w' /)
      character, dimension(16) ::
     &     ind=(/ 'a','b','c','d','i','j','k','l','p','q','r','s','t',
     &            'u','v','w' /)
      character(len=8) ::
     &     c1, c2, c3, c4,  ! Workspace
     &     a1, a2, a3, a4
      character(len=8), dimension(8) ::
      ! x_array(1:4) = creation operators (par/val/hol/f12)
      ! x_array(5:8) = annihilation operators (par/val/hol/f12)
     &     c_array,     ! Contraction index array
     &     e1_array,    ! External index of operator 1 array
     &     e2_array,    ! External index of operator 2 array
     &     t1_array,
     &     t2_array,
     &     tmp1_array, tmp2_array
      integer, dimension(2) ::
     &     idx_type     ! Info about index convention
      integer, parameter ::
     &     t_amp = 0,       ! [apij] (aacc)
     &     ham   = 1        ! [abip]
      character(len=20), dimension(4) ::
     &     tensor_ham=(/ 'H', 'INT_D', 'INT_HH', 'INT_PP' /)  ! Tensor names to use ham index convention
      integer ::
     &     conv(6), conv2(6)  ! Index convention arrays
      character(len=1) :: tmp ! Used in swapping pairs around
      integer :: ntmp
      real(4) :: factor

      type(pair) ::  test1
      integer :: shift, sp, ncre1, nann1, ncre2, nann2, ncre3, nann3
      integer :: ncre4, nann4 ! Number of ops in result tensor
      integer :: i1, i2 ! Search creation or annihilation first
      integer :: nc1, nc2, n1, n2, n3 ! Used to design generic function
      integer :: nloop ! Number of contraction loops
      integer :: r1, r2, r3 ! Ranks of tensors
      integer :: ext1, ext2 ! External index shifts
      logical :: found, sort
      type(pair_list) :: p_list, t1_list, t2_list, r_list, tmp_list
      ! Try new index assigment
!      type(tensor_slot) ::
!     &     op1, op2, ep1, ep2, cp1, cp2




      e1_array = '        '
      e2_array = '        '
      t1_array = '        '
      t2_array = '        '
      c_array = '        '



      ! Initalize tensor slots (old idea)
!      e1 = 0
!      do i = 1, contr_info%nj_op1
!        call count_index(i,
!     &     contr_info%occ_op1(1:,1:,i),
!     &     contr_info%rst_op1(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,e1)
!      end do
!      call construct_slot(e1,op1)
!
!      e1 = 0
!      do i = 1, contr_info%nj_op2
!        call count_index(i,
!     &     contr_info%occ_op2(1:,1:,i),
!     &     contr_info%rst_op2(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,e1)
!      end do
!      call construct_slot(e1,op2)
!
!      e1 = 0
!      do i = 1, contr_info%nj_op1
!        call count_index(i,
!     &     contr_info%occ_ex1(1:,1:,i),
!     &     contr_info%rst_ex1(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,e1)
!      end do
!      call construct_slot(e1,ep1)
!
!      e1 = 0
!      do i = 1, contr_info%nj_op2
!        call count_index(i,
!     &     contr_info%occ_ex2(1:,1:,i),
!     &     contr_info%rst_ex2(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,e1)
!      end do
!      call construct_slot(e1,ep2)
!
!      e1 = 0
!      do i = 1, contr_info%n_cnt
!        call count_index(i,
!     &     contr_info%occ_cnt(1:,1:,i),
!     &     contr_info%rst_cnt(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,e1)
!      end do
!      call construct_slot(e1,cp1)
!
!      e1 = 0
!      do i = 1, contr_info%nj_res
!        call count_index(i,
!     &     contr_info%occ_res(1:,1:,i),
!     &     contr_info%rst_res(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,e1)
!      end do
!      call construct_slot(e1,cp2)





      ! Set index convention
      idx_type=(/ 0, 0 /)
      do i=1, len(tensor_ham)
          if (contr_info%label_op1.eq.trim(tensor_ham(i))) then
              idx_type(1)=1
          end if
          if (contr_info%label_op2.eq.trim(tensor_ham(i))) then
              idx_type(2)=1
          end if
      end do

      c=0
      e1=0
      e2=0
      e3=0
      c1='        '
      c2='        '
      c3='        '
      a1='        '
      a2='        '
      a3='        '

      ! Get occuation info
      do i = 1, contr_info%n_cnt
        call count_index(i,
     &     contr_info%occ_cnt(1:,1:,i),
     &     contr_info%rst_cnt(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,c)
      end do
      do i = 1, contr_info%nj_op1
        call count_index(i,
     &     contr_info%occ_ex1(1:,1:,i),
     &     contr_info%rst_ex1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e1)
      end do
      do i = 1, contr_info%nj_op2
        call count_index(i,
     &     contr_info%occ_ex2(1:,1:,i),
     &     contr_info%rst_ex2(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e2)
      end do
      do i = 1, contr_info%nj_res
        call count_index(i,
     &     contr_info%occ_res(1:,1:,i),
     &     contr_info%rst_res(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e3)
      end do

      ! Figure out factor from equivalent lines
      factor = 1.0
      do j = 1, 2
         do i = 1, 4
            if (c(i,j) == 0) cycle
            if (mod(c(i,j),2) == 0) then
               factor = factor * (1.0/real(c(i,j)))
               !write(item%logfile,*) "FACT ", 1.0/c(i,j)
            end if
         end do
      end do
      !write(item%logfile,*) "FACT ", item%fact
      !write(item%logfile,*) "FACTOR ", factor
      !write(item%logfile,*) "FACTOR ", item%fact * factor
      item%fact = item%fact * factor


!      if (ep1%cslots /= 0) then
!         j = 1
!         do i=1, e1(2,1)
!             ep1%cre(i) = par(i)
!             j = j + 1
!         end do
!         do i=1, e1(3,1)
!             ep1%cre(j)=val(i)
!             j = j + 1
!         end do
!         do i=1, e1(1,1)
!             ep1%cre(j)=hol(i)
!             j = j + 1
!         end do
!      end if
!
!      if (ep1%aslots /= 0) then
!         j = 1
!         do i=1, e1(2,2)
!             ep1%ann(i) = par(i+e1(2,1))
!             j = j + 1
!         end do
!         do i=1, e1(3,2)
!             ep1%ann(j)=val(i+e1(3,1))
!             j = j + 1
!         end do
!         do i=1, e1(1,2)
!             ep1%ann(j)=hol(i+e1(1,1))
!             j = j + 1
!         end do
!      end if
!
!      if (ep2%cslots /= 0) then
!         j = 1
!         do i=1, e2(2,1)
!             ep2%cre(i)=par(i+e1(2,1)+e1(2,2))
!             j = j + 1
!         end do
!         do i=1, e2(3,1)
!             ep2%cre(j)=val(i+e1(3,1)+e1(3,2))
!             j = j + 1
!         end do
!         do i=1, e2(1,1)
!             ep2%cre(j)=hol(i+e1(1,1)+e1(1,2))
!             j = j + 1
!         end do
!      end if
!
!      if (ep2%aslots /= 0) then
!         j = 1
!         do i=1, e2(2,2)
!             ep2%ann(i)=par(i+e2(2,1)+e1(2,1)+e1(2,2))
!             j = j + 1
!         end do
!         do i=1, e2(3,2)
!             ep2%ann(j)=val(i+e2(3,1)+e1(3,1)+e1(3,2))
!             j = j + 1
!         end do
!         do i=1, e2(1,2)
!             ep2%ann(j)=hol(i+e2(1,1)+e1(1,1)+e1(1,2))
!             j = j + 1
!         end do
!      end if
!
!
!      write(item%logfile,*) "EP1 ", ep1%ann, ep1%aslots
!      write(item%logfile,*) "EP1 ", ep1%cre, ep1%cslots
!      write(item%logfile,*) "EP2 ", ep2%ann, ep2%aslots
!      write(item%logfile,*) "EP2 ", ep2%cre, ep2%cslots
!      write(item%logfile,*)
!
!      if (cp1%cslots /= 0) then
!         j = 1
!         do i=1, c(2,1)
!            cp1%cre(i) = cpar(i)
!            j = j + 1
!         end do
!         do i=1, c(3,1)
!            cp1%cre(j)=cval2(i)
!            j = j + 1
!         end do
!         do i=1, c(1,1)
!            cp1%cre(j)=chol(i)
!            j = j + 1
!         end do
!      end if
!
!      if (cp1%aslots /= 0) then
!         j = 1
!         do i=1, c(2,2)
!            cp1%ann(i) = cpar(i)
!            j = j + 1
!         end do
!         do i=1, c(3,2)
!            cp1%ann(j)= cval2(i)
!            j = j + 1
!         end do
!         do i=1, c(1,2)
!            cp1%ann(j)= chol(i)
!            j = j + 1
!         end do
!      end if


!      if (cp1%cslots /= 0) then
!         j = 1
!         do i=1, c(2,1)
!            k=i+e1(2,1)+e1(2,2)+e2(2,1)+e2(2,2)
!            cp1%cre(i) = par(k)
!            j = j + 1
!         end do
!         do i=1, c(3,1)
!            k=i+e1(3,1)+e1(3,2)+e2(3,1)+e2(3,2)
!            cp1%cre(j)=val(k)
!            j = j + 1
!         end do
!         do i=1, c(1,1)
!            k=i+e1(1,1)+e1(1,2)+e2(1,1)+e2(1,2)
!            cp1%cre(j)=hol(k)
!            j = j + 1
!         end do
!      end if
!
!      if (cp1%aslots /= 0) then
!         j = 1
!         do i=1, c(2,2)
!            k=i+c(2,1)+e1(2,1)+e1(2,2)+e2(2,1)+e2(2,2)
!            cp1%ann(i) = par(k)
!            j = j + 1
!         end do
!         do i=1, c(3,2)
!            k=i+c(3,1)+e1(3,1)+e1(3,2)+e2(3,1)+e2(3,2)
!            cp1%ann(j)=val(k)
!            j = j + 1
!         end do
!         do i=1, c(1,2)
!            k=i+c(1,1)+e1(1,1)+e1(1,2)+e2(1,1)+e2(1,2)
!            cp1%ann(j)=hol(k)
!            j = j + 1
!         end do
!      end if
!
!
!      if (ep1%cslots /=0) then
!         do i = 1, ep1%cslots
!            op1%cre(i) = ep1%cre(i)
!         end do
!      end if
!      if (ep1%aslots /=0) then
!         do i = 1, ep1%aslots
!            op1%ann(i) = ep1%ann(i)
!         end do
!      end if
!
!      if (ep2%cslots /=0) then
!         do i = 1, ep2%cslots
!            op2%cre(i) = ep2%cre(i)
!         end do
!      end if
!      if (ep2%aslots /=0) then
!         do i = 1, ep2%aslots
!            op2%ann(i) = ep2%ann(i)
!         end do
!      end if
!
!
!      if (cp1%cslots /=0) then
!         do i = 1, cp1%cslots
!            op1%cre(i+ep1%cslots) = cp1%cre(i)
!         end do
!      end if
!      if (cp1%aslots /=0) then
!         do i = 1, cp1%aslots
!            op1%ann(i+ep1%aslots) = cp1%ann(i)
!         end do
!      end if
!
!      if (cp1%aslots /=0) then
!         do i = 1, cp1%aslots
!            op2%cre(i+ep2%cslots) = cp1%ann(i)
!         end do
!      end if
!      if (cp1%cslots /=0) then
!         do i = 1, cp1%cslots
!            op2%ann(i+ep2%aslots) = cp1%cre(i)
!         end do
!      end if
!
!
!      if (cp2%cslots /= 0) then
!         j = 1
!         do i=1, e3(2,1)
!             cp2%cre(i) = par(i)
!             j = j + 1
!         end do
!         do i=1, e3(3,1)
!             cp2%cre(j)=val(i)
!             j = j + 1
!         end do
!         do i=1, e3(1,1)
!             cp2%cre(j)=hol(i)
!             j = j + 1
!         end do
!      end if
!
!      if (cp2%aslots /= 0) then
!         j = 1
!         do i=1, e3(2,2)
!             cp2%ann(i) = par(i+e3(2,1))
!             j = j + 1
!         end do
!         do i=1, e3(3,2)
!             cp2%ann(j)=val(i+e3(3,1))
!             j = j + 1
!         end do
!         do i=1, e3(1,2)
!             cp2%ann(j)=hol(i+e3(1,1))
!             j = j + 1
!         end do
!      end if


!      if (ep1%cslots + ep1%aslots > 0 .and.
!     &    ep2%cslots + ep2%aslots > 0) then
!         write(item%logfile,*) "pair on different ops, need contrac"
!      else
!         write(item%logfile,*) "pair on same ops, no contrac"
!      end if
!
!      ! Just use op1 and op2???
!      write(item%logfile,*) "RES ", cp2%ann, cp2%aslots
!      write(item%logfile,*) "RES ", cp2%cre, cp2%cslots
!      write(item%logfile,*) "OP1 ", op1%ann, op1%aslots
!      write(item%logfile,*) "OP1 ", op1%cre, op1%cslots
!      write(item%logfile,*) "OP2 ", op2%ann, op2%aslots
!      write(item%logfile,*) "OP2 ", op2%cre, op2%cslots

!!!      ! Need to think about multiple contraction pairs
!!!      write(item%logfile,*) "RES ", cp1%ann, cp1%aslots
!!!      write(item%logfile,*) "RES ", cp1%cre, cp1%cslots
!!!      cp1%pair(1) = ' '
!!!      cp1%pair(2) = ' '
!!!      if (cp1%cslots/=0 .or. cp1%aslots/=0) then
!!!      if (cp1%cslots + cp1%aslots > 1) then
!!!      do i = 1, cp1%cslots
!!!         if (cp1%cre(i) /= ' ') then
!!!            cp1%pair(1) = cp1%cre(i)
!!!            cp1%cre(i) = ' '
!!!!            cp1%cslots = cp1%cslots - 1
!!!            exit
!!!         end if
!!!      end do
!!!      do i = 1, cp1%aslots
!!!         if (cp1%ann(i) /= ' ') then
!!!            cp1%pair(2) = cp1%ann(i)
!!!            cp1%ann(i) = ' '
!!!!            cp1%aslots = cp1%aslots - 1
!!!            exit
!!!         end if
!!!      end do
!!!
!!!      write(item%logfile,*) "Pair ", cp1%pair
!!!      end if
!!!      end if
!!!
!!!      write(item%logfile,*) "RES ", cp1%ann, cp1%aslots
!!!      write(item%logfile,*) "RES ", cp1%cre, cp1%cslots
!!!
!!!      ! Check if external index on different op1s
!!!      if (ep1%cslots + ep1%aslots > 0 .and.
!!!     &    ep2%cslots + ep2%cslots > 0) then
!!!
!!!         write(item%logfile,*) "EP1 ", ep1%ann, ep1%aslots
!!!         write(item%logfile,*) "EP1 ", ep1%cre, ep1%cslots
!!!         write(item%logfile,*) "EP2 ", ep2%ann, ep2%aslots
!!!         write(item%logfile,*) "EP2 ", ep2%cre, ep2%cslots
!!!
!!!         ! If they are, assign a contraction index
!!!         cp1%pair2(1) = ' '
!!!         cp1%pair2(2) = ' '
!!!         if (ep1%cslots > ep2%cslots) then
!!!         do i = 1, ep1%cslots
!!!            if (ep1%cre(i) /= ' ') then
!!!               cp1%pair2(1) = ep1%cre(i)
!!!               ep1%cre(i) = ' '
!!!               ep1%cslots = ep1%cslots - 1
!!!               exit
!!!            end if
!!!         end do
!!!         do i = 1, ep2%aslots
!!!            if (ep2%ann(i) /= ' ') then
!!!               cp1%pair2(2) = ep2%ann(i)
!!!               ep2%ann(i) = ' '
!!!               ep2%aslots = ep2%aslots - 1
!!!               exit
!!!            end if
!!!         end do
!!!
!!!         else
!!!         do i = 1, ep2%cslots
!!!            if (ep2%cre(i) /= ' ') then
!!!               cp1%pair2(1) = ep2%cre(i)
!!!               ep2%cre(i) = ' '
!!!               ep2%cslots = ep2%cslots - 1
!!!               exit
!!!            end if
!!!         end do
!!!         do i = 1, ep1%aslots
!!!            if (ep1%ann(i) /= ' ') then
!!!               cp1%pair2(2) = ep1%ann(i)
!!!               ep1%ann(i) = ' '
!!!               ep1%aslots = ep1%aslots - 1
!!!               exit
!!!            end if
!!!         end do
!!!         write(item%logfile,*) "External Pair ", cp1%pair2
!!!         end if
!!!
!!!         ! Assign a contraction index to the external pair
!!!         cp1%contract = ' '
!!!         if (ep1%aslots > 0) then
!!!            write(item%logfile,*) "Hello1"
!!!            do i = 1, cp1%cslots
!!!               if (cp1%cre(i) /= ' ') then
!!!                  cp1%contract = cp1%cre(i)
!!!               end if
!!!            end do
!!!         else
!!!            write(item%logfile,*) "Hello2"
!!!            do i = 1, cp1%aslots
!!!               if (cp1%ann(i) /= ' ') then
!!!                  cp1%contract = cp1%ann(i)
!!!               end if
!!!            end do
!!!         end if
!!!         write(item%logfile,*) "Contraction ", cp1%contract
!!!
!!!      end if

!      call desctruct_slot(op1)
!      call desctruct_slot(op2)
!      call desctruct_slot(ep1)
!      call desctruct_slot(ep2)
!      call desctruct_slot(cp1)
!      call desctruct_slot(cp2)



      ! Order in ITF usually follows: apij
      ! Defualt [ccaa] as in the case of T[abij]

      ! Assign e1 (external indicies of t1)
      shift = 1
      do i=1, e1(2,1)
          c1(i:)=par(i)
          test1%pindex(shift) = par(i)
          shift = shift + 1
      end do
      e1_array(1)=c1
      do i=1, e1(3,1)
          c2(i:)=val(i)
          test1%pindex(shift) = val(i)
          shift = shift + 1
      end do
      e1_array(2)=c2
      do i=1, e1(1,1)
          c3(i:)=hol(i)
          test1%pindex(shift) = hol(i)
          shift = shift + 1
      end do
      e1_array(3)=c3

      ! Need to to be shifted to not match assignment of creations above
      do i=1, e1(2,2)
          a1(i:)=par(i+e1(2,1))
          test1%pindex(shift) = par(i)
          shift = shift + 1
      end do
      e1_array(5)=a1
      do i=1, e1(3,2)
          a2(i:)=val(i+e1(3,1))
          test1%pindex(shift) = val(i)
          shift = shift + 1
      end do
      e1_array(6)=a2
      do i=1, e1(1,2)
          a3(i:)=hol(i+e1(1,1))
          test1%pindex(shift) = hol(i)
          shift = shift + 1
      end do
      e1_array(7)=a3

!      write(item%logfile,*) "e1_array ", e1_array

      c1='        '
      c2='        '
      c3='        '
      a1='        '
      a2='        '
      a3='        '

      ! Shifted so as not to match e1 index
      do i=1, e2(2,1)
          c1(i:)=par(i+e1(2,1)+e1(2,2))
          test1%pindex(shift) = c1(i:)
          shift = shift + 1
      end do
      e2_array(1)=c1
      do i=1, e2(3,1)
          c2(i:)=val(i+e1(3,1)+e1(3,2))
          test1%pindex(shift) = c2(i:)
          shift = shift + 1
      end do
      e2_array(2)=c2
      do i=1, e2(1,1)
          c3(i:)=hol(i+e1(1,1)+e1(1,2))
          test1%pindex(shift) = c3(i:)
          shift = shift + 1
      end do
      e2_array(3)=c3

      ! Shifted so as not to match e1 index or above creations
      do i=1, e2(2,2)
          a1(i:)=par(i+e2(2,1)+e1(2,1)+e1(2,2))
          test1%pindex(shift) = a1(i:)
          shift = shift + 1
      end do
      e2_array(5)=a1
      do i=1, e2(3,2)
          a2(i:)=val(i+e2(3,1)+e1(3,1)+e1(3,2))
          test1%pindex(shift) = a2(i:)
          shift = shift + 1
      end do
      e2_array(6)=a2
      do i=1, e2(1,2)
          a3(i:)=hol(i+e2(1,1)+e1(1,1)+e1(1,2))
          test1%pindex(shift) = a2(i:)
          shift = shift + 1
      end do
      e2_array(7)=a3

      c1='        '
      c2='        '
      c3='        '
      a1='        '
      a2='        '
      a3='        '

      ! Permute indicies to get antisymm tensors
      if (item%permute==0 .or. item%permute==1) then
         ! No permutations
         do i=1, size(t1_array)
            t1_array(i)=e1_array(i)
            t2_array(i)=e2_array(i)
         end do
      else if (item%permute==2) then
         ! Permute creations
         do i=1, size(t1_array)/2
            !t1_array(i)=e2_array(i)
            !t1_array(i+4)=e1_array(i+4)
            !t2_array(i)=e1_array(i)
            !t2_array(i+4)=e2_array(i+4)

            t1_array(i)=e1_array(i)
            t1_array(i+4)=e2_array(i+4)
            t2_array(i)=e2_array(i)
            t2_array(i+4)=e1_array(i+4)
         end do
      end if


      ! Assign c (contracted by)
      ! These need to be shifted, so as not to match e1 or c2
      do i=1, c(2,1)
          c1(i:)=par(i+e1(2,1)+e1(2,2)+e2(2,1)+e2(2,2))
      end do
      c_array(1)=c1
      do i=1, c(3,1)
          c2(i:)=val(i+e1(3,1)+e1(3,2)+e2(3,1)+e2(3,2))
      end do
      c_array(2)=c2
      do i=1, c(1,1)
          c3(i:)=hol(i+e1(1,1)+e1(1,2)+e2(1,1)+e2(1,2))
      end do
      c_array(3)=c3

      ! Final shift so as not to match above creations
      do i=1, c(2,2)
          a1(i:)=par(i+c(2,1)+e1(2,1)+e1(2,2)+e2(2,1)+e2(2,2))
      end do
      c_array(5)=a1
      do i=1, c(3,2)
          a2(i:)=val(i+c(3,1)+e1(3,1)+e1(3,2)+e2(3,1)+e2(3,2))
      end do
      c_array(6)=a2
      do i=1, c(1,2)
          a3(i:)=hol(i+c(1,1)+e1(1,1)+e1(1,2)+e2(1,1)+e2(1,2))
      end do
      c_array(7)=a3




!======================================================================
      ! Start new idea for index (Jan 2019)
!======================================================================
      write(item%logfile,*) "E1:", e1_array
      write(item%logfile,*) "E2:", e2_array
      write(item%logfile,*) "C: ", c_array
      test1%pindex = ' '

      ! Allocate pair list, at most 8 pairs
      ! TODO: allocate based on ranks of tensors
      allocate(p_list%plist(8))
      allocate(t1_list%plist(4))
      allocate(t2_list%plist(4))
      allocate(r_list%plist(4))
      allocate(tmp_list%plist(4))

      do i = 1, 8
         p_list%plist(i)%pindex = ''
         p_list%plist(i)%link = ''
      end do
      do i = 1, 4
         t1_list%plist(i)%pindex = ''
         t1_list%plist(i)%link = ''
         t2_list%plist(i)%pindex = ''
         t2_list%plist(i)%link = ''
         r_list%plist(i)%pindex = ''
         r_list%plist(i)%link = ''
         tmp_list%plist(i)%pindex = ''
         tmp_list%plist(i)%link = ''
      end do

      ! Move occupation arrays to canonical order, p/h/v
      nc(1,1) = c(2,1)
      nc(2,1) = c(1,1)
      nc(3,1) = c(3,1)
      nc(4,1) = c(4,1)

      nc(1,2) = c(2,2)
      nc(2,2) = c(1,2)
      nc(3,2) = c(3,2)
      nc(4,2) = c(4,2)

      ne1(1,1) = e1(2,1)
      ne1(2,1) = e1(1,1)
      ne1(3,1) = e1(3,1)
      ne1(4,1) = e1(4,1)

      ne1(1,2) = e1(2,2)
      ne1(2,2) = e1(1,2)
      ne1(3,2) = e1(3,2)
      ne1(4,2) = e1(4,2)

      ne2(1,1) = e2(2,1)
      ne2(2,1) = e2(1,1)
      ne2(3,1) = e2(3,1)
      ne2(4,1) = e2(4,1)

      ne2(1,2) = e2(2,2)
      ne2(2,2) = e2(1,2)
      ne2(3,2) = e2(3,2)
      ne2(4,2) = e2(4,2)

      nr(1,1) = e3(2,1)
      nr(2,1) = e3(1,1)
      nr(3,1) = e3(3,1)
      nr(4,1) = e3(4,1)

      nr(1,2) = e3(2,2)
      nr(2,2) = e3(1,2)
      nr(3,2) = e3(3,2)
      nr(4,2) = e3(4,2)


      ! Find out number of creation/annihilation operators per operator
      ncre1 = 0
      ncre2 = 0
      nann1 = 0
      nann2 = 0
      ncre3 = 0
      nann3 = 0
      ncre4 = 0
      nann4 = 0
      do i = 1, 4
         ncre1 = ncre1 + ne1(i,1)
         ncre2 = ncre2 + ne2(i,1)
         nann1 = nann1 + ne1(i,2)
         nann2 = nann2 + ne2(i,2)
         ncre3 = ncre3 + nc(i,1)
         nann3 = nann3 + nc(i,2)
         ncre4 = ncre4 + nr(i,1)
         nann4 = nann4 + nr(i,2)
      end do

!      write(item%logfile,*) "ncre1: ", ncre1
!      write(item%logfile,*) "nann1: ", nann1
!      write(item%logfile,*) "ncre2: ", ncre2
!      write(item%logfile,*) "nann2: ", nann2
!      write(item%logfile,*) "ncre3: ", ncre3
!      write(item%logfile,*) "nann3: ", nann3

      ! Set ranks of tensors
      r1 = ncre1 + nann1 + ncre3 + nann3
      r2 = ncre2 + nann2 + ncre3 + nann3
      r3 = ncre4 + nann4
      write(item%logfile,*) "r1: ", r1
      write(item%logfile,*) "r2: ", r2
      write(item%logfile,*) "r3: ", r3


      ! Assign contraction loops
      ! Need to look for creations/annihilation first...

      ! Search the creation ops or annihilation ops first?
      if (ncre3 >= nann3) then
         ! Loop through creations first
         i1 = 1
         i2 = 2
      else
         ! Loop through annihilation first
         i1 = 2
         i2 = 1
      end if

      ! Pair list shift (shift pair)
      sp = 1

      do while (ncre3 /= 0 .and. nann3 /= 0)
         ! Loop over P/H/V
         do i=1, 3
            shift = 1
            do j=1, nc(i,i1)
               ii = 1+(4*(i-1))+ne1(i,1)+ne1(i,2)+ne2(i,1)+ne2(i,2)+sp-1
               ! Need extra shift if annihilator
               if (i1 == 2) ii = ii + nc(i,1)

               p_list%plist(sp)%pindex(shift) = ind(ii)
               shift = shift + 1

               if (i1 == 2) then
                  nann3 = nann3 - 1
               else
                  ncre3 = ncre3 - 1
               end if

               ! Mark which operator this index belongs to
               ! Contraction index always wrt the first operator
               p_list%plist(sp)%ops(1) = 1
!               write(item%logfile,*) "WHAT: ", ii, ind(ii)

               ! Look for matching operator
               do k=1, 3
                  do l=1, nc(k,i2)
                     ii = 1+(4*(k-1))+ne1(k,1)+ne1(k,2)+
     &                    ne2(k,1)+ne2(k,2) + sp-1
                     ! Need extra shift if annihilator
                     if (i1 == 1) ii = ii + nc(k,1)

                     p_list%plist(sp)%pindex(shift) = ind(ii)

                     if (i1 == 2) then
                        ncre3 = ncre3 - 1
                     else
                        nann3 = nann3 - 1
                     end if

                     p_list%plist(sp)%linked = .false.
                     p_list%plist(sp)%ops(2) = 1
                     write(item%logfile,*) "WHAT: ", ii, ind(ii)

                     ! Found a pair, so increment pair list index
                     sp = sp + 1
                     exit
                  end do
               end do

               ! A pair loop has been found so exit
               exit
            end do
            ! Can't match creation or annihilation so exit
            write(item%logfile,*) "ncre3: ", ncre3, "nann3", nann3
            if (ncre3 == 0 .or. nann3 == 0) exit
         end do
      end do

      ! Set number of contraction loops
      nloop = sp-1

      do i = 1, nloop
      write(item%logfile,*) "contraction loop: ", p_list%plist(i)%pindex
      end do

      ! TODO: need some sort of shift that keeps track of which index
      ! has been used....

      ! Match external pairs, either to external ops on the same
      ! operator, or to external ops on the second operator
      do while (ncre1 + nann1 /= 0 .or. ncre2 + nann2 /= 0)
      if (ncre1 + nann1 >= ncre2 + nann2 .and. ncre1 + nann1 /= 0) then
         if (ncre1 >= nann1) then
            ! Loop through creations first
            i1 = 1
            i2 = 2
            n1 = nann1
            n2 = nann2
            n3 = ncre1
         else
            ! Loop through annihilation first
            i1 = 2
            i2 = 1
            n1 = ncre1
            n2 = ncre2
            n3 = nann1
         end if

         do i=1, 3
            shift = 1
            do j=1, ne1(i,i1)
               ! Need sp-1-nloop to shift index so they aren't the same
               ! after a loop
               ii = 1+(4*(i-1)) + (sp-1-nloop)
               ! Need extra shift if this is annihilator
               if (i1 == 2) ii = ii + ne1(i,1)
               p_list%plist(sp)%pindex(shift) = ind(ii)
               p_list%plist(sp)%ops(1) = 1
               shift = shift + 1

               if (i1 == 2) then
                  nann1 = nann1 - 1
               else
                  ncre1 = ncre1 - 1
               end if

               ! Look for matching operator
               if (n1>0) then
                  ! Look on the same (first) operator
                  do k=1, 3
                     do l=1, ne1(k,i2)
                        ii = 1+(4*(k-1)) + (sp-1-nloop)
                        ! Need to shift extra places if this is
                        ! annihilation ops
                        if (i1 == 1) ii = ii + ne1(k,1)

                        p_list%plist(sp)%pindex(shift) = ind(ii)

                        if (i1 == 2) then
                            ncre1 = ncre1 - 1
                        else
                            nann1 = nann1 - 1
                        end if

                        p_list%plist(sp)%linked = .false.
                        p_list%plist(sp)%ops(2) = 1

                        ! Found a pair, so increment pair list index
                        sp = sp + 1
                        exit
                     end do
                     !TODO: this is a problem, need to exit when found
                     !pair. Need some sort of logical
                  end do
               else if (n2 > 0) then
                  ! Look on the second operator
                  do k=1, 3
                     do l=1, ne2(k,i2)
                        ii = 1+(4*(k-1))+ne1(k,i1)+ne1(k,i2) +
     &                       (sp-1 - nloop)
                        ! Need to shift for annihilation
                        if (i1 == 1) ii = ii + ne2(k,1)

                        p_list%plist(sp)%pindex(shift) = ind(ii)

                        if (i1 == 2) then
                            ncre2 = ncre2 - 1
                        else
                            nann2 = nann2 - 1
                        end if

                        p_list%plist(sp)%ops(2) = 2

                        ! We need to link the external indices on two
                        ! different operators with a contraction index
                        p_list%plist(sp)%linked = .true.
                        do m=1, 3
                           do n=1, nc(m,i2)
         ii = 1+(4*(m-1))+ne1(m,i1)+ne1(m,i2)+ne2(m,i1)+ne2(m,i2)+nloop
                              ! Need to shift for annihilation
                              if (i1 == 1) ii = ii + nc(m,1)

                              p_list%plist(sp)%link = ind(ii)

                              if (i1 == 2) then
                                  ncre3 = ncre3 - 1
                              else
                                  nann3 = nann3 - 1
                              end if
                              exit
                           end do
                        end do

                        ! Found a pair + made a link, so increment pair list index
                        sp = sp + 1
                        exit
                     end do
                  end do
               else
                  ! No ops of opposite type to match...
                  call line_error("Particle number not conserving",item)
               end if

               ! A pair loop has been found so exit
               exit
            end do
            ! Can't match creation or annihilation so exit
            !write(item%logfile,*) "NLOOP: ", nloop, "SP: ", sp-1
            if (n3 == 0) exit
         end do

      ! If the second operator has more ops, then search from there
      else if (ncre2 + nann2 /= 0)then
         write(item%logfile,*) "HELLO"

         if (ncre1 >= nann1) then
            ! Loop through creations first
            i1 = 1
            i2 = 2
            n1 = nann2
            n2 = ncre1
            n3 = ncre2
         else
            ! Loop through annihilation first
            i1 = 2
            i2 = 1
            n1 = ncre2
            n2 = nann1
            n3 = nann2
         end if

         do i=1, 3
            shift = 1
            ! Search second operator
            do j=1, ne2(i,i1)
               ii = 1+(4*(i-1))+ne1(i,i1)+ne1(i,i2)

               ! Need extra shift if annihilator
               if (i1 == 2) ii = ii + ne2(i,1)

               p_list%plist(sp)%pindex(shift) = ind(ii)
               p_list%plist(sp)%ops(1) = 2
               shift = shift + 1

               if (i1 == 2) then
                  nann2 = nann2 - 1
               else
                  ncre2 = ncre2 - 1
               end if

               ! Look for matching operator
               if (n1>0) then
                  ! Look on the same (second) operator
                  do k=1, 3
                     do l=1, ne2(k,i2)
                        ii = 1+(4*(k-1))+ne1(i,i1)+ne1(i,i2) +
     &                       (sp-1 - nloop)

                        ! Need extra shift if annihilator
                        if (i1 == 1) ii = ii + ne2(k,1)

                        p_list%plist(sp)%pindex(shift) = ind(ii)

                        if (i1 == 2) then
                            ncre2 = ncre2 - 1
                        else
                            nann2 = nann2 - 1
                        end if

                        p_list%plist(sp)%linked = .false.
                        p_list%plist(sp)%ops(2) = 2

                        ! Found a pair, so increment pair list index
                        sp = sp + 1
                        exit
                     end do
                  end do
               else if (n2 > 0) then
                  ! Look on the first operator
                  do k=1, 3
                     do l=1, ne1(k,i2)
                        ii = 1+(4*(k-1)) + (sp-1 - nloop)

                        ! Need to shift for annihilation
                        if (i1 == 1) ii = ii + ne1(k,1)

                        p_list%plist(sp)%pindex(shift) = ind(ii)

                        if (i1 == 2) then
                            ncre1 = ncre1 - 1
                        else
                            nann1 = nann1 - 1
                        end if

                        p_list%plist(sp)%ops(2) = 1

                        ! We need to link the external indices on two
                        ! different operators with a contraction index
                        p_list%plist(sp)%linked = .true.
                        do m=1, 3
                           do n=1, nc(m,i2)
         ii = 1+(4*(m-1))+ne1(m,i1)+ne1(m,i2)+ne2(m,i1)+ne2(m,i2)+nloop
                              ! Need to shift for annihilation
                              if (i1 == 1) ii = ii + nc(m,1)

                              p_list%plist(sp)%link = ind(ii)

                              if (i1 == 2) then
                                  ncre3 = ncre3 - 1
                              else
                                  nann3 = nann3 - 1
                              end if
                              exit
                           end do
                        end do

                        ! Found a pair, so increment pair list index
                        sp = sp + 1
                        exit
                     end do
                  end do
               else
                  ! No ops of opposite type to match...
                  call line_error("Particle number not conserving",item)
               end if

               ! A pair loop has been found so exit
               exit
            end do
            ! Can't match creation or annihilation so exit
            if (n3 == 0) exit
         end do
      end if
      end do   ! do while loop



      ! We now have list of pairs + which ops they belong to and any
      ! contraction indices linking external indices on different
      ! operators
      ! Now we must assign them to ITF index string, in the correct
      ! positions




      ! First arrange the result index made from external indices only
      ! TODO: Construct a canonical order for index, so intermediate
      ! indices match in the result and use
      ! TODO: Do not construct canonical order for integrals - need to
      ! pick them out

      ! IDEA: Create two pair_lists for op1 and op2
      ! Create list with just letters which belong to that op AND no
      ! links (they now belong to each op)
      ! Now have two lists to order with only the indices that belong to
      ! that op
      ! Swap between pairs
      ! Swap between pair_lists, using first index as a guide
      ! Thereby get canonical order but maintain index pairs
      ! Don't do this for integrals - need to be picked out in python
      ! to determine if they are K or J

      ! Create pair lists for each tensor
      shift = 1
      do i = 1, sp-1
         if (p_list%plist(i)%ops(1) == 1) then
            ! Add to t1_list + think about link
            t1_list%plist(shift)%pindex(1) = p_list%plist(i)%pindex(1)

            if (p_list%plist(i)%linked) then
             t1_list%plist(shift)%pindex(2) = p_list%plist(i)%link
            else
             t1_list%plist(shift)%pindex(2) = p_list%plist(i)%pindex(2)
            end if

            t1_list%plist(shift)%linked = .false.
            t1_list%plist(shift)%ops = 1
            shift = shift + 1
         else if (p_list%plist(i)%ops(2) == 1) then
            ! Add to t1_list + think about link
            t1_list%plist(shift)%pindex(2) = p_list%plist(i)%pindex(2)

            if (p_list%plist(i)%linked) then
             t1_list%plist(shift)%pindex(1) = p_list%plist(i)%link
            else
             t1_list%plist(shift)%pindex(1) = p_list%plist(i)%pindex(1)
            end if

            t1_list%plist(shift)%linked = .false.
            t1_list%plist(shift)%ops = 1
            shift = shift + 1
         end if
      end do

      shift = 1
      do i = 1, sp-1
         if (p_list%plist(i)%ops(1) == 2) then
            ! Add to t1_list + think about link
            t2_list%plist(shift)%pindex(1) = p_list%plist(i)%pindex(1)

            if (p_list%plist(i)%linked) then
             t2_list%plist(shift)%pindex(2) = p_list%plist(i)%link
            else
             t2_list%plist(shift)%pindex(2) = p_list%plist(i)%pindex(2)
            end if

            t2_list%plist(shift)%linked = .false.
            t2_list%plist(shift)%ops = 1
            shift = shift + 1
         else if (p_list%plist(i)%ops(2) == 2) then
            ! Add to t1_list + think about link
            t2_list%plist(shift)%pindex(2) = p_list%plist(i)%pindex(2)

            if (p_list%plist(i)%linked) then
             t2_list%plist(shift)%pindex(1) = p_list%plist(i)%link
            else
             t2_list%plist(shift)%pindex(1) = p_list%plist(i)%pindex(1)
            end if

            t2_list%plist(shift)%linked = .false.
            t2_list%plist(shift)%ops = 1
            shift = shift + 1
         end if
      end do

      ! Create pair list for result tensor, only need external index
      shift = 1
      do i = nloop+1, sp-1
         r_list%plist(shift) = p_list%plist(i)
         shift = shift + 1
      end do



      ! Swap between pairs to canonical order
      ! TODO: will not work with pqrstu...
      do i = 1, r1/2
       if (.not. item%int(1)) then
       if (t1_list%plist(i)%pindex(1) > t1_list%plist(i)%pindex(2)) then
            tmp = t1_list%plist(i)%pindex(1)
            t1_list%plist(i)%pindex(1) = t1_list%plist(i)%pindex(2)
            t1_list%plist(i)%pindex(2) = tmp
         end if
       end if
      end do

      do i = 1, r2/2
       if (.not. item%int(2)) then
       if (t2_list%plist(i)%pindex(1) > t2_list%plist(i)%pindex(2)) then
            tmp = t2_list%plist(i)%pindex(1)
            t2_list%plist(i)%pindex(1) = t2_list%plist(i)%pindex(2)
            t2_list%plist(i)%pindex(2) = tmp
         end if
       end if
      end do

      do i = 1, r3/2
       if (.not. item%int(3)) then
       if (r_list%plist(i)%pindex(1) > r_list%plist(i)%pindex(2)) then
            tmp = r_list%plist(i)%pindex(1)
            r_list%plist(i)%pindex(1) = r_list%plist(i)%pindex(2)
            r_list%plist(i)%pindex(2) = tmp
         end if
       end if
      end do



      ! Sort index pairs into order with bubble sort
      ! If only 1 pair, then this is be skipped
      ! TODO: Dont need all of tmp_list
!      do i = 1, r1/2
!         write(item%logfile,*) "T1 LIST: ", t1_list%plist(i)%pindex
!      end do
!      do i = 1, r2/2
!         write(item%logfile,*) "T2 LIST: ", t1_list%plist(i)%pindex
!      end do

      ! Sort out the first tensor
      if (r1 > 2 .and. .not.item%int(1)) then
      sort = .true.
      do while (sort)
      sort = .false.
      do i = 1, r1/2
         if (i == (r1/2)) exit
         if (t1_list%plist(i+1)%pindex(1) < t1_list%plist(i)%pindex(1))
     &      then
            tmp_list%plist(1) = t1_list%plist(i)
            t1_list%plist(i) = t1_list%plist(i+1)
            t1_list%plist(i+1) = tmp_list%plist(1)
            sort = .true.
         end if
      end do
      end do
      end if

      ! Sort out the second tensor
      if (r2 > 2 .and. .not.item%int(2)) then
      sort = .true.
      do while (sort)
      sort = .false.
      do i = 1, r2/2
         if (i == (r2/2)) exit
         if (t2_list%plist(i+1)%pindex(1) < t2_list%plist(i)%pindex(1))
     &      then
            tmp_list%plist(1) = t2_list%plist(i)
            t2_list%plist(i) = t2_list%plist(i+1)
            t2_list%plist(i+1) = tmp_list%plist(1)
            sort = .true.
         end if
      end do
      end do
      end if

      ! Sort out the second tensor
      if (r3 > 2 .and. .not.item%int(3)) then
      sort = .true.
      do while (sort)
      sort = .false.
      do i = 1, r3/2
         if (i == (r3/2)) exit
         if (r_list%plist(i+1)%pindex(1) < r_list%plist(i)%pindex(1))
     &      then
            tmp_list%plist(1) = r_list%plist(i)
            r_list%plist(i) = r_list%plist(i+1)
            r_list%plist(i+1) = tmp_list%plist(1)
            sort = .true.
         end if
      end do
      end do
      end if

      ! Insert ordered lists into ITF index strings
      ! TODO: change these variables c1, just using now for convince
      c1 = '        '
      c2 = '        '
      c3 = '        '

      do i = 1, r1/2
         c1(i:i) = t1_list%plist(i)%pindex(1)
         c1(i+(r1/2):i+(r1/2)) = t1_list%plist(i)%pindex(2)
      end do
      do i = 1, r2/2
         c2(i:i) = t2_list%plist(i)%pindex(1)
         c2(i+(r2/2):i+(r2/2)) = t2_list%plist(i)%pindex(2)
      end do
      do i = 1, r3/2
         c3(i:i) = r_list%plist(i)%pindex(1)
         c3(i+(r3/2):i+(r3/2)) = r_list%plist(i)%pindex(2)
      end do





!      do i = 1, r1/2
!         write(item%logfile,*) "NEW T1 LIST: ", t1_list%plist(i)%pindex
!      end do
!      do i = 1, r2/2
!         write(item%logfile,*) "NEW T2 LIST: ", t2_list%plist(i)%pindex
!      end do






!      do i = nloop+1, sp-1
!         write(item%logfile,*) "externals: ", p_list%plist(i)%pindex
!         if (p_list%plist(i)%linked) then
!            write(item%logfile,*) "link: ", p_list%plist(i)%link
!         end if
!         write(item%logfile,*) "ops: ", p_list%plist(i)%ops
!      end do


      ! Create ITF index strings from pair lists
      ! Result tensor first
      ! (sp-1-nloop)*2 is the rank of the result tensor, ie. the number
      ! of external pairs multiplied by 2

!      c1 = '        '
!      do i = nloop+1, sp-1
!         c1(i:i) = p_list%plist(i)%pindex(1)
!         c1(i+(sp-1-nloop):i+(sp-1-nloop)) = p_list%plist(i)%pindex(2)
!      end do
!
!      ! Insert externals into the two tensors
!      c2 = '        '
!      c3 = '        '
!      ext1 = 1
!      ext2 = 1
!      do i = nloop+1, sp-1
!         if (p_list%plist(i)%ops(1)==1) then
!            ! Index belongs on the first tensor
!            c2(i:i) = p_list%plist(i)%pindex(1)
!            if (p_list%plist(i)%linked) then
!               c2(i+(r1/2):i+(r1/2)) = p_list%plist(i)%link
!            else
!               c2(i+(r1/2):i+(r1/2)) = p_list%plist(i)%pindex(2)
!            end if
!            ext1 = ext1 + 1
!
!            if (p_list%plist(i)%ops(2)==2) then
!               ! Second index belongs on second tensor
!               c3(i:i) = p_list%plist(i)%pindex(2)
!               if (p_list%plist(i)%linked) then
!                  c3(i+(r2/2):i+(r2/2)) = p_list%plist(i)%link
!               else
!                  c3(i+(r2/2):i+(r2/2)) = p_list%plist(i)%pindex(1)
!               end if
!               ext2 = ext2 + 1
!            end if
!         else
!            ! Index belongs on the second tensor
!            c3(i:i) = p_list%plist(i)%pindex(1)
!            if (p_list%plist(i)%linked) then
!               c3(i+(r2/2):i+(r2/2)) = p_list%plist(i)%link
!            else
!               c3(i+(r2/2):i+(r2/2)) = p_list%plist(i)%pindex(2)
!            end if
!            ext2 = ext2 + 1
!
!            if (p_list%plist(i)%ops(2)==1) then
!               c2(i:i) = p_list%plist(i)%pindex(2)
!               if (p_list%plist(i)%linked) then
!                  c2(i+(r1/2):i+(r1/2)) = p_list%plist(i)%link
!               else
!                  c2(i+(r1/2):i+(r1/2)) = p_list%plist(i)%pindex(1)
!               end if
!               ext1 = ext1 + 1
!            end if
!         end if
!      end do
!
!      ! Insert contraction indices into the two tensors
!      do i = 1, nloop
!         c2(i+ext1:i+ext1) = p_list%plist(i)%pindex(1)
!         c2(i+ext1+(r1/2):i+ext1+(r1/2)) = p_list%plist(i)%pindex(2)
!
!         c3(i+ext2:i+ext2) = p_list%plist(i)%pindex(2)
!         c3(i+ext2+(r2/2):i+ext2+(r2/2)) = p_list%plist(i)%pindex(1)
!      end do




      write(item%logfile,*) "=============================="
      write(item%logfile,*) "RESULT: ", trimal(c3)
      write(item%logfile,*) "T1: ", trimal(c1)
      write(item%logfile,*) "T2: ", trimal(c2)
      write(item%logfile,*) "=============================="


      deallocate(p_list%plist)
      deallocate(t1_list%plist)
      deallocate(t2_list%plist)
      deallocate(r_list%plist)
      deallocate(tmp_list%plist)


!======================================================================
      ! End new idea for index (Jan 2019)
!======================================================================






      ! Construct final index strings. For different tensors, there are
      ! different orders in which to place the operators. We pick out
      ! these cases and set a convention) array to give the order.
      ! Order of letters in x_array: a, p, i, x
      !write(item%logfile,*) "T1 array", t1_array
      !write(item%logfile,*) "C array", c_array
      !write(item%logfile,*) "T2 array", t2_array

      ! Operator 1
      select case(idx_type(1))
      case(ham)
         ! Hamiltonian/integral convention
         ! Check if K[ijab] -> K[abij]
         e1 = 0
         do i = 1, contr_info%nj_op1
           call count_index(i,
     &        contr_info%occ_op1(1:,1:,i),
     &        contr_info%rst_op1(1:,1:,1:,1:,1:,i),
     &        contr_info%ngas,contr_info%nspin,e1)
         end do

         if (e1(1,1)==2 .and. e1(2,2)==2) then
            conv = (/ 5, 6, 7, 1, 2, 3 /)
         else
            conv = (/ 3, 1, 2, 5, 6, 7 /)
         end if

      case default
         conv = (/ 1, 2, 3, 5, 6, 7 /)
      end select
      call index_convention(t1_array,c_array,item%idx1,conv,conv)

      ! Operator 2
      ! c_array annhilations correspond to t2 creations and vice versa
      select case(idx_type(2))
      case(ham)
         e2 = 0
         do i = 1, contr_info%nj_op2
           call count_index(i,
     &        contr_info%occ_op2(1:,1:,i),
     &        contr_info%rst_op2(1:,1:,1:,1:,1:,i),
     &        contr_info%ngas,contr_info%nspin,e2)
         end do

         if (e2(1,1)==2 .and. e2(2,2)==2) then
            conv =  (/ 5, 6, 7, 1, 2, 3 /)
            conv2 = (/ 1, 2, 3, 5, 6, 7 /)
         else
            conv =  (/ 3, 1, 2, 5, 6, 7 /)
            conv2 = (/ 5, 6, 7, 3, 1, 2 /)
         end if

      case default
         conv =  (/ 1, 2, 3, 5, 6, 7 /)
         conv2 = (/ 5, 6, 7, 1, 2, 3 /)
      end select
      call index_convention(t2_array,c_array,item%idx2,conv,conv2)

      !write(item%logfile,*) "INDEX1 ", item%idx1
      !write(item%logfile,*) "INDEX2 ", item%idx2
      !call correct_index(item%idx1)
      !call correct_index(item%idx2)
      !write(item%logfile,*) "INDEX1 ", item%idx1
      !write(item%logfile,*) "INDEX2 ", item%idx2

!      conv = (/ 1, 2, 3, 5, 6, 7 /)
!      call index_convention(e1_array,e2_array,tmp,conv,conv)
!
!      write(item%logfile,*) "TMP ", tmp
!      write(item%logfile,*) "1 ", item%idx1
!      write(item%logfile,*) "2 ", item%idx2
!      do i = 1, len(tmp)
!         do j = 1, len(item%idx2)
!            if (tmp(i:i) == item%idx2(j:j)) then
!               write(item%logfile,*) "HELLO2 ", j
!               ntmp(j:j) = item%idx2(j:j)
!            end if
!         end do
!      end do
!
!      do i = 1, len(tmp)
!         do j = 1, len(item%idx1)
!            if (tmp(i:i) == item%idx1(j:j)) then
!               write(item%logfile,*) "HELLO1 ", j
!            end if
!         end do
!      end do
!
!      write(item%logfile,*) "E1 ", e1_array
!      write(item%logfile,*) "E2 ", e2_array






      ! Result, problems occur rank 2 tensors...
      conv = (/ 1, 2, 3, 5, 6, 7 /)
      if (len(trim(item%idx1)) < len(trim(item%idx2))) then
         call index_convention(e2_array,e1_array,item%idx3,conv,conv)
      else
         call index_convention(e1_array,e2_array,item%idx3,conv,conv)
      end if


      ! Need to correct index and check indicies are correctly paired
      ! together


      !call correct_index_pairs(item)





      ! Need to corrrect index for R[ai]; the a and i must be in the
      ! same slot
      ! This should be considered as a hack....
      if (len(trim(item%idx3)) == 2) then
      if (len(trim(item%idx1)) + len(trim(item%idx2))/=4) then
      if (.not. item%int(1)) then

         j = 0
         k = 0

         do i = 1, len(trim(item%idx1))
            if (item%idx3(1:1) == item%idx1(i:i)) then
               j = i
            end if
            if (item%idx3(2:2) == item%idx1(i:i)) then
               k = i
            end if
         end do

         do i = 1, len(trim(item%idx2))
            if (item%idx3(1:1) == item%idx2(i:i)) then
               j = i
            end if
            if (item%idx3(2:2) == item%idx2(i:i)) then
               k = i
            end if
         end do

         if (mod(j+k,2)>0) then
            if (item%inter(1) .and. len(trim(item%idx2))==4) then
               item%idx2 = t_index(item%idx2)
            else if (item%inter(2) .and. len(trim(item%idx1))==4) then
               item%idx1 = t_index(item%idx1)
            else
               item%idx1 = t_index(item%idx1)
            end if
         end if

      end if
      end if
      end if




      return
      end

!*----------------------------------------------------------------------*
!      subroutine correct_index_pairs(item)
!*----------------------------------------------------------------------*
!!
!*----------------------------------------------------------------------*
!
!      use itf_utils
!      implicit none
!      include 'opdim.h'
!      include 'def_contraction.h'
!      include 'def_itf_contr.h'
!
!      type(itf_contr), intent(inout) ::
!     &     item
!
!      character(len=1) ::
!     &     p1(2),p2(2),p3(2)
!      character(len=:), allocatable ::
!     &     pair_list(:,:)
!
!      integer :: length, pairs, rank1, rank2, rank3, i,j
!      logical :: found
!
!      rank3 = (len(trim(item%idx3))/2)
!      if (rank3 == 0) return
!      rank1 = (len(trim(item%idx1))/2)
!      rank2 = (len(trim(item%idx2))/2)
!
!      length = 1
!      pairs = (len(trim(item%idx1)) + len(trim(item%idx2))) / 2
!
!
!      allocate(character(length) :: pair_list(pairs,2))
!
!      pair_list = ' '
!
!      do i = 1, rank1
!         pair_list(i,2) = item%idx1(i:i)
!         pair_list(i,1) = item%idx1(rank1+i:rank1+i)
!      end do
!
!      do i = 1, rank2
!         pair_list(i+rank1,2) = item%idx2(i:i)
!         pair_list(i+rank1,1) = item%idx2(rank2+i:rank2+i)
!      end do
!
!!      write(item%logfile,*) "PAIR LIST "
!!      do i = 1, 2
!!         write(item%logfile,*) (pair_list(j,i), j=1,pairs)
!!      end do
!
!      p1(1) = item%idx3(1:1)
!      p1(2) = item%idx3(rank3+1:rank3+1)
!      p2 = p1
!
!      !write(item%logfile,*) "PAIR ", p1
!
!      found = .false.
!      ! TODO: check top row as well
!      do i = 1, pairs
!         if (p1(1) == pair_list(i,2)) then
!            p1(1) = pair_list(i,1)
!            !write(item%logfile,*) "FOUND ", p1(1)
!            do j = 1, pairs
!               !write(item%logfile,*) pair_list(j,2)
!               if (p1(1) == pair_list(j,2) .and.
!     &             p1(2) == pair_list(j,1)) then
!                  p1(1) = pair_list(j,1)
!               end if
!               if (p1(1)==p1(2)) then
!                  found = .true.
!                  exit
!               end if
!            end do
!            exit
!         end if
!      end do
!
!      if (.not. found) then
!      p1 = p2
!      do i = 1, pairs
!         if (p1(1) == pair_list(i,2)) then
!            p1(1) = pair_list(i,1)
!            !write(item%logfile,*) "FOUND ", p1(1)
!            do j = 1, pairs
!               !write(item%logfile,*) pair_list(j,1)
!               if (p1(1) == pair_list(j,1) .and.
!     &             p1(2) == pair_list(j,2)) then
!                  p1(1) = pair_list(j,2)
!               end if
!               if (p1(1)==p1(2)) then
!                  found = .true.
!                  exit
!               end if
!            end do
!            exit
!         end if
!      end do
!      end if
!
!      if (.not. found) then
!      p1 = p2
!      do i = 1, pairs
!         if (p1(1) == pair_list(i,1)) then
!            p1(1) = pair_list(i,2)
!            !write(item%logfile,*) "FOUND ", p1(1)
!            do j = 1, pairs
!               !write(item%logfile,*) pair_list(j,2)
!               !write(item%logfile,*) p1
!               if (p1(1) == pair_list(j,2) .and.
!     &             p1(2) == pair_list(j,1)) then
!                  p1(1) = pair_list(j,1)
!               end if
!               if (p1(1)==p1(2)) then
!                  found = .true.
!                  exit
!               end if
!            end do
!            exit
!         end if
!      end do
!      end if
!
!      if (.not. found) then
!      p1 = p2
!      do i = 1, pairs
!         if (p1(1) == pair_list(i,1)) then
!            p1(1) = pair_list(i,2)
!            !write(item%logfile,*) "FOUND ", p1(1)
!            do j = 1, pairs
!               !write(item%logfile,*) pair_list(j,1)
!               !write(item%logfile,*) p1
!               if (p1(1) == pair_list(j,1) .and.
!     &             p1(2) == pair_list(j,2)) then
!                  p1(1) = pair_list(j,2)
!               end if
!               if (p1(1)==p1(2)) then
!                  found = .true.
!                  exit
!               end if
!            end do
!            exit
!         end if
!      end do
!      end if
!
!      !if (found) write(item%logfile,*) "GOOD ", p1
!      !if (.not. found) write(item%logfile,*) "Transpose ", p1
!      !write(item%logfile,*) item%idx3, item%idx1, item%idx2
!      !write(item%logfile,*)
!      if (.not. found) then
!      if (rank1 > rank2) then
!         !write(item%logfile,*) item%idx3, t_index(item%idx1), item%idx2
!         item%idx1 = t_index(item%idx1)
!      else if (rank2 > rank1) then
!         !write(item%logfile,*) item%idx3, item%idx1, t_index(item%idx2)
!         item%idx2 = t_index(item%idx2)
!      else
!         !write(item%logfile,*) item%idx3,item%idx1,t_index(item%idx2)
!         if (item%int(1)) then
!            item%idx1 = t_index(item%idx1)
!         else
!            item%idx2 = t_index(item%idx2)
!         end if
!      end if
!      end if
!
!
!      deallocate(pair_list)
!
!
!      return
!      end


*----------------------------------------------------------------------*
      subroutine index_convention(arr1,arr2,item_idx,conv,conv2)
*----------------------------------------------------------------------*
!     Arrange index letters into correct order according to convention
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=8), dimension(8), intent(in) ::
     &     arr1,        ! Index array
     &     arr2
      integer, intent(in) ::
     &     conv(6),
     &     conv2(6)     ! TODO: Make this optional? Could get it to work...
      character(len=index_len), intent(inout) ::
     &     item_idx

      item_idx=trimal(arr1(conv(1)))//trimal(arr2(conv2(1)))//
     &         trimal(arr1(conv(2)))//trimal(arr2(conv2(2)))//
     &         trimal(arr1(conv(3)))//trimal(arr2(conv2(3)))//
     &         trimal(arr1(conv(4)))//trimal(arr2(conv2(4)))//
     &         trimal(arr1(conv(5)))//trimal(arr2(conv2(5)))//
     &         trimal(arr1(conv(6)))//trimal(arr2(conv2(6)))

      return
      end

*----------------------------------------------------------------------*
      subroutine correct_index(item_idx)
*----------------------------------------------------------------------*
!     Arrange index letters into correct order according to convention
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=index_len), intent(inout) ::
     &     item_idx

      character(len=4) ::
     &     tmp

      tmp = item_idx
      if (len(trim(item_idx)) == 4) then

!      if (item_idx(1:1) > item_idx(2:2) .and.
!     &    item_idx(2:2) /= 'b') then
!
!         item_idx(1:1) = item_idx(2:2)
!         item_idx(2:2) = tmp(1:1)
!      end if
!      if (item_idx(3:3) < item_idx(4:4) .and.
!     &    item_idx(3:3) /= 'i') then
!         item_idx(3:3) = item_idx(4:4)
!         item_idx(4:4) = tmp(3:3)
!      end if

      if (item_idx(1:1) == 'b') then
         item_idx(1:1) = item_idx(2:2)
         item_idx(2:2) = tmp(1:1)
      end if

      if (item_idx(3:3) == 'j') then
         item_idx(3:3) = item_idx(4:4)
         item_idx(4:4) = tmp(3:3)
      end if

      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine permute_tensors(contr_info,perm_case,lulog)
*----------------------------------------------------------------------*
!     Find permutation case
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(binary_contr), intent(in) ::
     &     contr_info   ! Information about binary contraction
      integer, intent(inout) ::
     &     perm_case
      integer, intent(in) ::
     &     lulog

      integer ::
     &     e1(4,2),      ! Occupations of external index 1
     &     e2(4,2),      ! Occupations of external index 2
     &     c(4,2)
      integer ::
     &     i,
     &     sum_c1,sum_c2,sum_a1,sum_a2

      ! Check if not antisym over different verticies
      ! Check for tensor products
      ! Not going to antisymm intermediates...
      ! The intermediates can have eeaa structure
      e1=0
      e2=0
      c=0

      ! Get occuation info
      do i = 1, contr_info%n_cnt
        call count_index(i,
     &     contr_info%occ_cnt(1:,1:,i),
     &     contr_info%rst_cnt(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,c)
      end do
      do i = 1, contr_info%nj_op1
        call count_index(i,
     &     contr_info%occ_ex1(1:,1:,i),
     &     contr_info%rst_ex1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e1)
      end do
      do i = 1, contr_info%nj_op2
        call count_index(i,
     &     contr_info%occ_ex2(1:,1:,i),
     &     contr_info%rst_ex2(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e2)
      end do

      ! C: |i|a|p|x|
      ! A: |i|a|p|x|

      if ((sum(sum(e1,dim=1))+sum(sum(e2,dim=1)))==2) then
         return
      else if ((sum(sum(e1,dim=1))+sum(sum(e2,dim=1)))==0) then
         return
      else if ((sum(sum(e1,dim=1))+sum(sum(e2,dim=1)))==6) then
         return
      end if

      perm_case = 0

      if (e1(2,1)+e2(2,1)==2 .and. e1(1,2)+e2(1,2)==2 .or.
     &    e1(3,1)+e2(3,1)==2 .and. e1(1,2)+e2(1,2)==2 .or.
     &    e1(2,1)+e2(2,1)==2 .and. e1(3,2)+e2(3,2)==2) then

         sum_c1=0
         sum_c2=0
         sum_a1=0
         sum_a2=0

         do i=1, 4
            ! Sum creation ops
            sum_c1=sum_c1+e1(i,1)
            sum_c2=sum_c2+e2(i,1)

            ! Summ annhilation ops
            sum_a1=sum_a1+e1(i,2)
            sum_a2=sum_a2+e2(i,2)
         end do

         ! If sum==2, then both indicies come from same operator, therefore
         ! it doesn't need symmetrised
         if (sum_c1/=2 .and. sum_c2/=2) then
            if (sum_c1+sum_c2==2) then
               perm_case = 1
            end if
         end if

         if (sum_a1/=2 .and. sum_a2/=2) then
            if (sum_a1+sum_a2==2) then
               if (perm_case == 1) then
                  ! P(ab), then P(abij)
                  ! When we have (1-P(ab))(1-P(ij)) then we only need to
                  ! generate 1 + P(ab)
                  perm_case = 2
               else
                  perm_case = 1
               end if
            end if
         end if

      else
         return
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine index_to_groups(cov,contv,index,half_rank)
*----------------------------------------------------------------------*
!     Assign index to covariant and contravarient groups
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     half_rank         ! Rank of tensor divided by 2
      character(len=1), intent(inout) ::
     &     cov(half_rank),      ! Covarient indicies
     &     contv(half_rank)    ! Contravarient indicies
      character(len=*), intent(in) ::
     &     index       ! ITF tensor index

      integer ::
     &     i            ! Loop index

      do i = 1, half_rank
         cov(i) = index(i:i)
         contv(i) = index(i+half_rank:i+half_rank)
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine assign_spin(item)
*----------------------------------------------------------------------*
!     Assign spin to tensors, then sum remaining contraction indices
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &     item

      character(len=1) ::
     &     t1a(2),   ! First group of first tensor
     &     t1b(2),   ! Second group of first tensor
     &     t2a(2),
     &     t2b(2),
     &     r1a(2),
     &     r1b(2)
      integer ::
     &     i,j,      ! Loop index
     &     sum_c1,sum_c2,sum_a1,sum_a2,
     &     s1a(2),
     &     s1b(2),
     &     s2a(2),
     &     s2b(2),
     &     second_idx,third_idx,fourth_idx,
     &     zero_a,zero_b,
     &     result_spin(4) ! Spin case of result
      logical ::
     &     logi ! Debug delete

      ! There are certain cases that don't need to be spin summed
      ! because the spin summed line is the same as the non-spin summed
      ! line.
      if (item%rank1 + item%rank2 + item%rank3 == 4) then
         ! This covers all cases where we have two rank-2 tensors
         ! + one rank-0 tensor somewhere on a line

         ! If the line involves an intermediate, then we must add the
         ! spin name to the intermedite name. For all these cases, this
         ! is spimple
         if (item%inter(1)) then
            item%label_t1 = trim(item%label_t1)//'aa'
         else if (item%inter(2)) then
            item%label_t2 = trim(item%label_t2)//'aa'
         end if

         ! The spin case is the line itself and wont contain any other
         ! terms, so mark the start and end
         call print_itf_line(item,.false.,.false.)

         return
      else if (item%rank1 + item%rank2 + item%rank3 == 6) then
         ! This covers all cases where we have three rank-2 tensors

         if (item%inter(1)) then
            item%label_t1 = trim(item%label_t1)//'aa'
         else if (item%inter(2)) then
            item%label_t2 = trim(item%label_t2)//'aa'
         end if

         call print_itf_line(item,.false.,.false.)
         return
      else if (item%rank1 + item%rank2 + item%rank3 == 0) then
         ! Scalar contributions
         call print_itf_line(item,.false.,.false.)
         return
      else if (item%rank3==6 .or. item%rank1==6 .or. item%rank2==6) then
         ! Don't care about higher rank contractions for now...
         call print_itf_line(item,.false.,.false.)
         return
      else if (item%rank3==4 .and. item%rank1==2
     &         .and. item%rank2==2) then
         ! Don't care about tensor products now
         call print_itf_line(item,.false.,.false.)
         return
      end if

      s1a=0
      s1b=0
      s2a=0
      s2b=0
      t1a=''
      t1b=''
      t2a=''
      t2b=''
      r1a=''
      r1b=''

      ! Start by splitting the index into covarient and contravarient
      ! groups
      if (item%rank1 >= item%rank2) then
         call index_to_groups(t1a,t1b,item%idx1,item%rank1/2)
         call index_to_groups(t2a,t2b,item%idx2,item%rank2/2)
      else
         ! Swapping tensors, so the largest rank tensor is now in t1a
         ! and t1b
         call index_to_groups(t2a,t2b,item%idx1,item%rank1/2)
         call index_to_groups(t1a,t1b,item%idx2,item%rank2/2)
         item%swapped = .true.
      end if

      ! Get result index
      call index_to_groups(r1a,r1b,item%idx3,item%rank3/2)

      ! Each result has a specific spin case, ie. for rank-2 we just
      ! need the aa case, for rank-4 we need abab. So find the external
      ! indicies in the contraction tensors and set their spins. For
      ! certain intermedites we may need different cases, so set
      ! result_spin accordingly.
      if (item%inter(3)) then
         result_spin = item%spin_case
      else
         ! abab
         result_spin = (/ 1, 2, 1, 2 /)
      end if

!      write(item%logfile,*) "RESULT_SPIN: ", result_spin

      do i=1, 2
         do j=1, 2
            ! Assign spin of first tensor
            if (r1a(j)==t1a(i) .and. r1a(j)/='') then
               s1a(i) = result_spin(j)
            else if (r1a(j)==t1b(i) .and. r1a(j)/='') then
               s1b(i) = result_spin(j)
            end if

            if (r1b(j)==t1a(i) .and. r1b(j)/='') then
               s1a(i) = result_spin(j+2)
            else if (r1b(j)==t1b(i) .and. r1b(j)/='') then
               s1b(i) = result_spin(j+2)
            end if

            ! Assign spin of second tensor
            if (r1a(j)==t2a(i) .and. r1a(j)/='') then
               s2a(i) = result_spin(j)
            else if (r1a(j)==t2b(i) .and. r1a(j)/='') then
               s2b(i) = result_spin(j)
            end if

            if (r1b(j)==t2a(i) .and. r1b(j)/='') then
               s2a(i) = result_spin(j+2)
            else if (r1b(j)==t2b(i) .and. r1b(j)/='') then
               s2b(i) = result_spin(j+2)
            end if
         end do
      end do

      ! Index assigned from result
      ! Debug
      !write(item%logfile,*) "s1b: ", s1b
      !write(item%logfile,*) "s1a: ", s1a
      !write(item%logfile,*)
      !write(item%logfile,*) "s2b: ", s2b
      !write(item%logfile,*) "s2a: ", s2a
      !write(item%logfile,*)
      !write(item%logfile,*)

      ! The sum will start with the index group with the most unassigned
      ! contraction indcies (0s)
      zero_a=0
      zero_b=0
      do i=1, size(s1a)
         if (s1a(i)==0) then
            zero_a=zero_a+1
         end if
         if (s1b(i)==0) then
            zero_b=zero_b+1
         end if
      end do

      if (zero_a>=zero_b) then
         logi=.true.
         call spin_sum(s1a,s1b,s2a,s2b,t1a,t1b,t2a,t2b,zero_a,
     &                 zero_b,logi,item)
      else
         logi=.false.
         call spin_sum(s1b,s1a,s2a,s2b,t1b,t1a,t2a,t2b,zero_b,
     &                 zero_a,logi,item)
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine spin_sum(s1a,s1b,s2a,s2b,t1a,t1b,t2a,t2b,zero_a,
     &                    zero_b,logi,item)
*----------------------------------------------------------------------*
!     Spin sum
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(inout) ::
     &     s1a(2),
     &     s1b(2),
     &     s2a(2),
     &     s2b(2)
      character(len=1), intent(in) ::
     &     t1a(2),   ! First group of first tensor
     &     t1b(2),   ! Second group of first tensor
     &     t2a(2),
     &     t2b(2)
      integer, intent(in) ::
     &     zero_a,zero_b
      logical, intent(in) ::
     &     logi      ! Debug delete
      type(itf_contr), intent(inout) ::
     &     item

      integer ::
     &     i,j,k,l,m,n,o,p,q,r,s,
     &     second_idx,third_idx,fourth_idx,
     &     first_idx
      logical ::
     &     eloop     ! Check if at least one spin case is printed out

      first_idx=0
      second_idx=0
      third_idx=0
      fourth_idx=0

      eloop=.false.

      ! Check for first unassigned index
      if (any(s1a==0) .or. first_idx>0) then
      do i=1, size(s1a)
       if (s1a(i)==0 .or. first_idx>0 .and. i==first_idx) then
        if (first_idx==0) then
         first_idx=i
        end if
        do j=1, 2
         s1a(first_idx)=j
         call find_idx(s2a,s2b,t1a,t2a,t2b,first_idx,j)

         ! Check for second index
         if (any(s1a==0) .or. second_idx>0) then
          do k=i, size(s1a)
           if (s1a(k)==0 .or. second_idx>0 .and. k==second_idx) then
            if (second_idx==0) then
             second_idx=k
            end if
            do l=1, 2
             s1a(second_idx)=l
             call find_idx(s2a,s2b,t1a,t2a,t2b,second_idx,l)

             ! Check for third contraction index
             if (any(s1b==0) .or. third_idx>0) then
              do o=1, size(s1b)
               if (s1b(o)==0 .or. third_idx>0 .and. o==third_idx) then
                if (third_idx==0) then
                 third_idx=o
                end if
                do p=1, 2
                 s1b(third_idx)=p
                 call find_idx(s2a,s2b,t1b,t2a,t2b,third_idx,p)

                 ! Check for fourth contraction index
                 if (any(s1b==0) .or. fourth_idx>0) then
                  do q=1, size(s1b)
                   if (s1b(q)==0.or.fourth_idx>0.and.q==fourth_idx)
     &             then
                    if (fourth_idx==0) then
                     fourth_idx=q
                    end if
                    do r=1, 2
                     s1b(fourth_idx)=r
                     call find_idx(s2a,s2b,t1b,t2a,t2b,
     &                             fourth_idx,r)
                     ! Print result, sum over four indicies
                     call print_spin_case(s1b,s1a,s2a,
     &                                    s2b,logi,eloop,item)
                    end do
                   end if
                  end do
                 else if (.not. any(s1b==0) .and. fourth_idx==0) then
                  ! Print result, sum over three indicies
                  call print_spin_case(s1b,s1a,s2a,s2b,logi,
     &                                 eloop,item)
                 end if

                end do
               end if
              end do
             else if (.not. any(s1b==0) .and. third_idx==0) then
              ! Print result, sum over two indicies
              call print_spin_case(s1a,s1b,s2a,s2b,.not.logi,
     &                             eloop,item)
             end if
            end do
           end if

          end do ! Check if second index is 0
         else if (zero_a==1 .and. zero_b==1) then
          ! Case where one contraction over s1a and s1b
          do o=1, size(s1b)
           if (s1b(o)==0 .or. third_idx>0 .and. o==third_idx) then
            if (third_idx==0) then
             third_idx=o
            end if
            do p=1, 2
             s1b(third_idx)=p
             call find_idx(s2a,s2b,t1b,t2a,t2b,third_idx,p)
             call print_spin_case(s1a,s1b,s2a,s2b,.not.logi,
     &                            eloop,item)
            end do
           end if
          end do
         else if (.not. any(s1a==0) .and. second_idx==0) then
          ! Print result, sum over one indicies
          call print_spin_case(s1a,s1b,s2a,s2b,.not.logi,
     &                         eloop,item)
         end if
        end do ! Loop over a/b for first index
       end if ! Check for first 0 index
      end do
      else if (.not. any(s1a==0) .and. first_idx==0) then
       ! Print result, sum over no index
       call print_spin_case(s1a,s1b,s2a,s2b,.not.logi,
     &                      eloop,item)
      end if


      ! loop check for 0
      !  if 0
      !  |loop 1, 2
      !  | loop check for 0
      !  |  if 0
      !  |   loop 1, 2
      !  |  else
      !  |   print
      !  else
      !   print

      if (.not. eloop) then
         call line_error("Didn't print out spin case", item)
      end if

      ! Mark the end of the spin summed block, if we will print a
      ! permutation of this line next (ie. permute==1) then we will not
      ! end the block just yet. This will save load/drop calls for the
      ! same tensors
      if (eloop .and. item%permute /= 1 .and. item%print_line) then
         if (.not. item%symm) write(item%logfile,*) "END"
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine print_spin_case(s1a,s1b,s2a,s2b,logi,eloop,item)
*----------------------------------------------------------------------*
!     Print spin case
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &     s1a(2),
     &     s1b(2),
     &     s2a(2),
     &     s2b(2)
      logical, intent(in) ::
     &     logi      ! Debug delete
      logical, intent(inout) ::
     &     eloop     ! Check if at least one spin case is printed out
      type(itf_contr), intent(inout) ::
     &     item

      logical ::
     &     s1,       ! True if tensor 1 is mixed spin
     &     s2        ! True if tensor 2 is mixed spin
      integer ::
     &     r1,r2,i,shift,
     &     ishift
      character(len=index_len) ::
     &     spin_name


       s1=.false.
       s2=.false.

      ! Determine rank; may have swapped t1 and t2 to move larger tensor
      ! to t1
      if (item%swapped) then
         r1 = item%rank2
         r2 = item%rank1
      else
         r1 = item%rank1
         r2 = item%rank2
      end if

      !write(item%logfile,*) "s1b: ", s1b
      !write(item%logfile,*) "s1a: ", s1a
      !write(item%logfile,*)
      !write(item%logfile,*) "s2b: ", s2b
      !write(item%logfile,*) "s2a: ", s2a
      !write(item%logfile,*)
      !write(item%logfile,*)


       ! Pick out specific spin cases here
       if (sum(s1a)==sum(s1b) .and.
     &     sum(s2a)==sum(s2b)) then
         if (modulo(sum(s1a)+sum(s1b),2)==0 .and.
     &       modulo(sum(s2a)+sum(s2b),2)==0) then

            ! TODO: Doesn't work for rank 6 tensors yet...
            if (r1==2 .or. r1==0) then
               s1=.false.
            else if (s1a(1)/=s1a(2)) then
               s1=.false.
            else if (s1a(1)==s1a(2) .and. s1a(1)/=0) then
               ! Pure spin
               s1=.true.
            end if

            if (r2==2 .or. r2==0) then
               s2=.false.
            else if (s2a(1)/=s2a(2)) then
               s2=.false.
            else if (s2a(1)==s2a(2) .and. s2a(1)/=0) then
               s2=.true.
            end if

            !DEBUG
            !if (logi) then
            !   write(item%logfile,*) "s1b: ", s1b
            !   write(item%logfile,*) "s1a: ", s1a
            !   write(item%logfile,*)
            !   write(item%logfile,*) "s2b: ", s2b
            !   write(item%logfile,*) "s2a: ", s2a
            !   write(item%logfile,*)
            !   write(item%logfile,*)
            !else
            !   write(item%logfile,*) "s1b: ", s1a
            !   write(item%logfile,*) "s1a: ", s1b
            !   write(item%logfile,*)
            !   write(item%logfile,*) "s2b: ", s2b
            !   write(item%logfile,*) "s2a: ", s2a
            !   write(item%logfile,*)
            !   write(item%logfile,*)
            !end if


            ! Change name of intermediates involved in a result line
            spin_name = '        '
            if (item%inter(1)) then
               if (item%swapped) then
                  do i = 1, item%rank1/2
                     if (s2a(i)==1) then
                        spin_name(i:i) = 'a'
                     else if (s2a(i)==2) then
                        spin_name(i:i) = 'b'
                     end if

                     if (s2b(i)==1) then
                        spin_name(i+item%rank1/2:i+item%rank1/2) = 'a'
                     else if (s2b(i)==2) then
                        spin_name(i+item%rank1/2:i+item%rank1/2) = 'b'
                     end if
                  end do
               else
                  do i = 1, item%rank1/2
                     if (s1a(i)==1) then
                        spin_name(i:i) = 'a'
                     else if (s1a(i)==2) then
                        spin_name(i:i) = 'b'
                     end if

                     if (s1b(i)==1) then
                        spin_name(i+item%rank1/2:i+item%rank1/2) = 'a'
                     else if (s1b(i)==2) then
                        spin_name(i+item%rank1/2:i+item%rank1/2) = 'b'
                     end if
                  end do
               end if
               item%inter1 = spin_name

            else if (item%inter(2)) then
               if (item%swapped) then
                  do i = 1, item%rank2/2
                     if (s1a(i)==1) then
                        spin_name(i:i) = 'a'
                     else if (s1a(i)==2) then
                        spin_name(i:i) = 'b'
                     end if

                     if (s1b(i)==1) then
                        spin_name(i+item%rank2/2:i+item%rank2/2) = 'a'
                     else if (s1b(i)==2) then
                        spin_name(i+item%rank2/2:i+item%rank2/2) = 'b'
                     end if
                  end do
               else
                  do i = 1, item%rank2/2
                     if (s2a(i)==1) then
                        spin_name(i:i) = 'a'
                     else if (s2a(i)==2) then
                        spin_name(i:i) = 'b'
                     end if

                     if (s2b(i)==1) then
                        spin_name(i+item%rank2/2:i+item%rank2/2) = 'a'
                     else if (s2b(i)==2) then
                        spin_name(i+item%rank2/2:i+item%rank2/2) = 'b'
                     end if
                  end do
               end if
               item%inter2 = spin_name
            end if

!            write(item%logfile,*) "SPIN1: ", item%inter1
!            write(item%logfile,*) "SPIN2: ", item%inter2

            ! Check if we have to deal with intermediates
            ! TODO: what if there are two intermediates on one line??
            if (associated(item%inter_spins)) then

               ! Number of intermediates
               ishift = 0

               ! Check if the first tensor is an intermediate
               if (item%inter(1)) then
                  ishift = ishift + 1
                  ! Number of spin cases
                  shift = item%inter_spins(ishift)%ncase + 1
                  !write(item%logfile,*)
                  !write(item%logfile,*) "intermeiate", item%label_t1

                  ! For now, if intermeidiate is part of a permutation
                  ! line, then we add a P to its name
                  if (item%permute == 2) then
                  item%inter_spins(ishift)%name=trim(item%label_t1)//'P'
                  else
                     item%inter_spins(ishift)%name=item%label_t1
                  end if

                  if (item%swapped) then
                     ! t1 and t2 were swapped in summation
                     do i=1, 2
                        item%inter_spins(ishift)%cases(i,shift)=s2a(i)
                        item%inter_spins(ishift)%cases(i+2,shift)=s2b(i)
                     end do
                  else
                     do i=1, 2
                        item%inter_spins(ishift)%cases(i,shift)=s1a(i)
                        item%inter_spins(ishift)%cases(i+2,shift)=s1b(i)
                     end do
                  end if
               end if

               ! Check if the second tensor is an intermediate
               if (item%inter(2)) then
                  ishift = ishift + 1
                  ! Number of spin cases
                  shift = item%inter_spins(ishift)%ncase + 1
                  !write(item%logfile,*)
                  !write(item%logfile,*) "intermediate", item%label_t2
                  if (item%permute == 2) then
                     item%inter_spins(ishift)%name=item%label_t2//'P'
                  else
                     item%inter_spins(ishift)%name=item%label_t2
                  end if

                  if (item%swapped) then
                     ! t1 and t2 were swapped in summation
                     do i=1, 2
                        item%inter_spins(ishift)%cases(i,shift)=s1a(i)
                        item%inter_spins(ishift)%cases(i+2,shift)=s1b(i)
                     end do
                  else
                     do i=1, 2
                        item%inter_spins(ishift)%cases(i,shift)=s2a(i)
                        item%inter_spins(ishift)%cases(i+2,shift)=s2b(i)
                     end do
                  end if
               end if

!               do i=1, 4
!                  write(item%logfile,*) "SC ",
!     &            item%inter_spins(1)%cases(i,shift)
!               end do

!               write(item%logfile,*) "Name: ", item%inter_spins(1)%name
               ! Update number of spin cases for each different
               ! intermediate
               do i = 1, ishift
                  item%inter_spins(ishift)%ncase =
     &                              item%inter_spins(ishift)%ncase + 1
               end do

               item%ninter = ishift
            end if


            ! Print the spin summed line
            if (item%print_line) then
               ! Mark the start of the spin summed block, if the current
               ! line is a permutation of the previous line (ie. permute >1)
               ! then we should not begin a new block
               if (eloop==.false. .and. .not. (item%permute > 1)) then
                  write(item%logfile,*) 'BEGIN'
               end if

               if (item%swapped) then
                  ! t1 and t2 were swapped in summation
                  call print_itf_line(item,s2,s1)
               else
                  call print_itf_line(item,s1,s2)
               end if
            end if

            eloop=.true.
         end if
       end if


       return
       end


*----------------------------------------------------------------------*
      subroutine find_idx(s2a,s2b,t1,t2a,t2b,idx,j)
*----------------------------------------------------------------------*
!     Find corresponding index in second tenor
*----------------------------------------------------------------------*

      implicit none

      integer, intent(inout) ::
     &     s2a(2),
     &     s2b(2)
      integer, intent(in) ::
     &     idx,       ! Position of contraction index
     &     j        ! Spin value of contraction index
      character(len=1), intent(in) ::
     &     t1(2),   ! First group of first tensor
     &     t2a(2),
     &     t2b(2)
      integer ::
     &     m

      do m=1, size(t2a)
         if (t2a(m)==t1(idx)) then
            s2a(m)=j
         else if (t2b(m)==t1(idx)) then
            s2b(m)=j
         end if
      end do

      return
      end

*----------------------------------------------------------------------*
      subroutine itf_contr_init(contr_info,itf_item,perm,comm,
     &                          lulog)
*----------------------------------------------------------------------*
!     Initialise ITF contraction object
*----------------------------------------------------------------------*

      use itf_utils
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h' ! For def_formular_item.h
      include 'def_contraction.h'
      include 'def_formula_item.h' ! For command parameters
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &     contr_info   ! Information about binary contraction
      type(itf_contr), intent(inout) ::
     &     itf_item     ! Object which holds information necessary to print out an ITF algo line
      integer, intent(in) ::
     &     perm,        ! Permutation information
     &     comm,        ! formula_item command
     &     lulog        ! Output file
      ! This should be in a 'constructor'

      ! Assign output file
      itf_item%logfile=lulog

      ! Assign command type
      itf_item%command=comm

      ! Assign permutation number
      if (comm/=command_cp_intm .or. comm/=command_add_intm) then
         itf_item%permute=perm
      end if

      ! Assign labels
      itf_item%label_t1=contr_info%label_op1
      if (comm/=command_cp_intm .or. comm/=command_add_intm) then
         ! Operator 2 does not exist in [ADD] or [COPY]
         itf_item%label_t2=contr_info%label_op2
      end if
      itf_item%label_res=contr_info%label_res


      ! Check if an intermediate
      !call check_inter(itf_item%label_t1,itf_item%inter(1))
      !if (comm/=command_cp_intm .or. comm/=command_add_intm) then
      !   call check_inter(itf_item%label_t2,itf_item%inter(2))
      !end if
      !call check_inter(itf_item%label_res,itf_item%inter(3))

      itf_item%inter(1) = check_inter(itf_item%label_t1)
      if (comm/=command_cp_intm .or. comm/=command_add_intm) then
         itf_item%inter(2) = check_inter(itf_item%label_t2)
      end if
      itf_item%inter(3) = check_inter(itf_item%label_res)

      ! Check if an integral
      itf_item%int(1) = check_int(itf_item%label_t1)
      if (comm/=command_cp_intm .or. comm/=command_add_intm) then
         itf_item%int(2) = check_int(itf_item%label_t2)
      end if
      itf_item%int(3) = .false.

      ! Assign factor --- use special ITF factor
      ! the ITF factor is closer to the value expected from standard
      ! diagram rules, but still some care has to be taken when translating
      ! to ITF tensors; e.g. the (HP;HP) integrals are stored by GeCCo as
      ! <aj||bi> while the standard sign for ring terms assumes <aj||ib>
      itf_item%fact=contr_info%fact_itf

      ! Account for negative sign as explained from above...
      call integral_fact(contr_info,itf_item%fact)

      ! Assign index string
      if (comm==command_cp_intm .or. comm==command_add_intm) then
         ! For [ADD] and [COPY]
         call assign_add_index(contr_info,itf_item)
      else
         ! For other contractions
         call assign_index(contr_info,itf_item)
      end if

      ! Assign rank
      call itf_rank(itf_item%idx1,itf_item%rank1)
      if (comm/=command_cp_intm .or. comm/=command_add_intm) then
         call itf_rank(itf_item%idx2,itf_item%rank2)
      end if
      call itf_rank(itf_item%idx3,itf_item%rank3)

      itf_item%inter1 = '        '
      itf_item%inter2 = '        '

      return
      end


*----------------------------------------------------------------------*
      subroutine integral_fact(contr,fact)
*----------------------------------------------------------------------*
!     Check if negative factor needs to be applied due to the different
!     storage of the <aj||bi> integrals
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'

      type(binary_contr), intent(in) ::
     &     contr   ! Information about binary contraction
      real, intent(inout) ::
     &     fact

      integer ::
     &     i,
     &     c(4,2)


      c = 0

      if (contr%label_op1 == 'H') then
         do i = 1, contr%n_cnt
           call count_index(i,
     &        contr%occ_op1(1:,1:,i),
     &        contr%rst_op1(1:,1:,1:,1:,1:,i),
     &        contr%ngas,contr%nspin,c)
         end do
      else if (contr%label_op2 == 'H') then
         do i = 1, contr%n_cnt
           call count_index(i,
     &        contr%occ_op2(1:,1:,i),
     &        contr%rst_op2(1:,1:,1:,1:,1:,i),
     &        contr%ngas,contr%nspin,c)
         end do
      end if

      if (c(1,1)==1 .and. c(2,1)==1 .and. c(1,2)==1 .and. c(2,2)==1)
     &   then
            fact = fact * -1.0
      end if

      ! Catch three external integrals
      if (c(1,1) + c(1,2)==3)
     &   then
            fact = fact * -1.0
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine print_symmetrise(result, item)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=maxlen_bc_label), intent(inout) ::
     &     result
      type(itf_contr), intent(in) ::
     &     item

      character(len=70) ::
     &     line
      character(len=index_len) ::
     &     tindex
      character(len=maxlen_bc_label) ::
     &     new

      new = rename_tensor(result, item%rank3)

      line = '.'//trim(new)//'['//trim(item%idx3)//'] += '//
     &       trim(item%label_res)//'['//trim(item%idx3)//']'
      write(item%logfile,*) trim(line)

      tindex = ' '
      tindex(1:1) = item%idx3(2:2)
      tindex(2:2) = item%idx3(1:1)
      tindex(3:3) = item%idx3(4:4)
      tindex(4:4) = item%idx3(3:3)

      line = '.'//trim(new)//'['//trim(item%idx3)//'] += '//
     &       trim(item%label_res)//'['//trimal(tindex)//']'
      write(item%logfile,*) trim(line)

      ! Mark end of spin block
      write(item%logfile,*) 'END '

      return
      end

*----------------------------------------------------------------------*
      subroutine itf_rank(idx,nrank)
*----------------------------------------------------------------------*
!     Calculate rank of tensor
*----------------------------------------------------------------------*

      ! TODO: Possibly set rank in assign_index - so don't need this...
      !       Also think about assign_add_index
      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=index_len), intent(in) ::
     &     idx
      integer, intent(inout) ::
     &     nrank

      nrank=len(trimal(idx))

      return
      end
