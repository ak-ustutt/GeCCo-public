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
!     Rename tensor according to taste
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
!     Return array with number of operators of each type
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     iocc(ngastp,2), ngas, nspin,
     &     irstr(2,ngas,2,2,nspin), idx
      integer, intent(inout) ::
     &     nops(4,2)                     ! Matrix of index info

      ! Creation operators
      ! Hole, abcd
      nops(1,1)=nops(1,1) + iocc(2,1)
      ! Particle, ijkl
      nops(2,1)=nops(2,1) + iocc(1,1)
      ! Valence, pqrs
      nops(3,1)=nops(3,1) + iocc(3,1)
      ! Explicit, x
      nops(4,1)=nops(4,1) + iocc(4,1)

      ! Annihilation operators as above
      nops(1,2)=nops(1,2) + iocc(2,2)
      nops(2,2)=nops(2,2) + iocc(1,2)
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
     &     contr_info      ! Information about binary contraction
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

      ! Initialise permutation factors to 0 == no permutation
      perm_case = 0

      ! Determine if result needs permuting
      inter = check_inter(contr_info%label_res)

      if (.not.inter) then
         call permute_tensors(contr_info,perm_case,itflog)
      end if

      ! If the perm_array doesn't contain any zeros, then we should
      ! introduce an intermediate which collects half of the different
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
               ! Loop over permutation cases and send separately to
               ! assign_spin. For most cases this is just one, however
               ! for (1-Pij)(1-Pab), we need to generate one of these
               ! permutations before symmetrising
               call itf_contr_init(contr_info,itf_item,i,command,itflog)

               if (i == 2) then
                  ! Need to transpose by tensors after permutation, to
                  ! avoid symmetry problem when using (1 + Pabij)
                  itf_item%idx1 = t_index(itf_item%idx1)
                  itf_item%idx2 = t_index(itf_item%idx2)
               end if

               call assign_spin(itf_item)
            end do

            ! If created a perm intermediate, print the symmetrised lines
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


      ! Initialise permutation factors to 0 == no permutation
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
      ! need to create this new intermediate. This requires the
      ! transpose
      if (scan('P', label)) then
         itf_item%idx3 = t_index(itf_item%idx3)
         if (itf_item%rank2 > 2) then
            ! Don't need to permute T[ai] etc.
            itf_item%idx2 = t_index(itf_item%idx2)
         end if
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
            ! Because we take the transpose of the spin orbital eqns twice,
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
      if (abs(abs(item%fact) - 1.0d+0) > 1.0d-15) then
            write(sfact,*) item%fact

            if (item%fact < 0.0d+0) then
               do i = 1, len(sfact)
                  if (sfact(i:i) == '-') then
                     ! Remove leading negative sign
                     sfact(i:i) = ''
                     exit
                  end if
               end do
            end if

            do i=1, len(sfact)
               ! Remove leading zero
               if (sfact(i:i) == '0') then
                  sfact(i:i)=''
               else if (sfact(i:i) == '.') then
                  exit
               end if
            end do

            do i=len(sfact), 1, -1
               ! Remove trailing zeros
               if (sfact(i:i) == '0' .or. sfact(i:i) == ' ') then
                  sfact(i:i) = ''
               else
                  exit
               end if
            end do
            sfact_star=trimal(sfact)//'*'
      end if

      ! Determine what the contraction operator looks like
      equal_op='  '
      if (item%fact < 0.0d+0) then
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
!     Simple ITF index assignment for lines that only contain a result
!     and a tensor, ie. COPY and ADD lines. Don't need to pair indices.
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &     contr_info   ! Information about binary contraction
      type(itf_contr), intent(inout) ::
     &     item         ! ITF binary contraction

      integer ::
     &     o1(4,2),     ! Operator numbers of the first tensor (T1)
     &     i            ! Loop index
      character, dimension(4) ::
     &     hol=(/ 'i','j','k','l' /),
     &     par=(/ 'a','b','c','d' /)
      character, dimension(8) ::
     &     val=(/ 'p','q','r','s','t','u','v','w' /)
      character(len=index_len) ::
     &     c1, c2, c3,
     &     a1, a2, a3
      character(len=index_len), dimension(8) ::
     &     o1_array     ! Index of operator 1

      o1=0

      ! Get occupation info
      do i = 1, contr_info%nj_op1
        call count_index(i,
     &     contr_info%occ_op1(1:,1:,i),
     &     contr_info%rst_op1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,o1)
      end do

      ! Assign ranks of tensors
      call itf_rank(o1, o1, item%rank1, .true.)
      item%rank3 = item%rank1

      c1='        '
      c2='        '
      c3='        '
      a1='        '
      a2='        '
      a3='        '

      ! Assign o1 (external indices of t1)
      do i=1, o1(1,1)
          c1(i:)=par(i)
      end do
      o1_array(1)=c1
      do i=1, o1(3,1)
          c2(i:)=val(i)
      end do
      o1_array(2)=c2
      do i=1, o1(2,1)
          c3(i:)=hol(i)
      end do
      o1_array(3)=c3

      ! Need to to be shifted to not match assignment of creations above
      do i=1, o1(1,2)
          a1(i:)=par(i+o1(1,1))
      end do
      o1_array(5)=a1
      do i=1, o1(3,2)
          a2(i:)=val(i+o1(3,1))
      end do
      o1_array(6)=a2
      do i=1, o1(2,2)
          a3(i:)=hol(i+o1(2,1))
      end do
      o1_array(7)=a3

      item%idx1=trimal(o1_array(1))//trimal(o1_array(2))//
     &          trimal(o1_array(3))//trimal(o1_array(5))//
     &          trimal(o1_array(6))//trimal(o1_array(7))

      item%idx3=item%idx1

      return
      end


*----------------------------------------------------------------------*
      subroutine assign_index(contr_info,item)
*----------------------------------------------------------------------*
!     Assign an ITF index string to each tensor in a line
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &   contr_info     ! Information about binary contraction
      type(itf_contr), intent(inout) ::
     &   item           ! ITF binary contraction

      integer ::
     &   c(4,2),        ! Operator numbers of contraction index
     &   e1(4,2),       ! Operator numbers of external index 1
     &   e2(4,2),       ! Operator numbers of external index 2
     &   e3(4,2),       ! Operator numbers of result index
     &   e4(4,2),       ! Operator numbers of result index
     &   e5(4,2),       ! Operator numbers of result index
     &   i, j, k, l,    ! Loop index
     &   ii,            ! Letter index
     &   nloop,         ! Number of contraction loops
     &   i1, i2,        ! Search creation or annihilation first
     &   shift,         ! List shift
     &   sp             ! Pair list shift
      character, dimension(16) ::
     &   ind=(/ 'a','b','c','d','i','j','k','l','p','q','r','s','t',
     &          'u','v','w' /)   ! Letters for index string
      character(len=index_len) ::
     &     s1, s2, s3   ! Tmp ITF index strings
      real(8) ::
     &   factor         ! Factor from equivalent lines
      logical ::
     &   found,         ! True if found pairing index
     &   sort           ! Used in bubble sort
      type(pair_list) ::
     &   p_list,        ! Complete list of pairs in binary contraction
     &   t1_list,       ! List of pairs in first tensor (T1)
     &   t2_list,       ! List of pairs in second tensor (T2)
     &   r_list         ! List of pairs in result tensor
      integer, dimension(4) ::
     &   t_shift,       ! Index shift for external indices
     &   c_shift        ! Index shift for contraction indices
      integer, dimension(2) ::
     &   cops,          ! Number of creation/annihilation contraction operators
     &   e1ops,         ! Number of C/A external T1 operators
     &   e2ops,         ! Number of C/A external T2 operators
     &   e3ops          ! Total number of external operators

      c=0
      e1=0
      e2=0
      e3=0

      e4=0
      e5=0

      ! Get occupation info
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

      do i = 1, contr_info%nj_op1
        call count_index(i,
     &     contr_info%occ_op1(1:,1:,i),
     &     contr_info%rst_op1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e4)
      end do
      e5 = e3

      ! Figure out factor from equivalent lines
      factor = 1.0d+0
      do j = 1, 2
         do i = 1, 4
            if (c(i,j) == 0) cycle
            if (mod(c(i,j),2) == 0) then
               factor = factor * (1.0d+0/real(c(i,j),8))
            end if
         end do
      end do
      item%fact = item%fact * factor

      ! Multiply factor by -1.0 due to permutation
      if (item%permute == 2) then
         item%fact = item%fact * -1.0d+0
      end if

      ! Find out number of creation/annihilation operators per operator
      cops = sum(c, dim=1)
      e1ops = sum(e1, dim=1)
      e2ops = sum(e2, dim=1)
      e3ops = sum(e3, dim=1)

      ! Set ranks of tensors
      call itf_rank(e1, c, item%rank1, .false.)
      call itf_rank(e2, c, item%rank2, .false.)
      call itf_rank(e3, c, item%rank3, .true.)

      ! Allocate pair list according to ranks of tensors
      allocate(p_list%plist(item%rank1+item%rank2))
      allocate(t1_list%plist(item%rank1))
      allocate(t2_list%plist(item%rank2))
      allocate(r_list%plist(item%rank3))

      ! To start, we pair off contraction loops
      ! Pair list index (shift pair)
      sp = 1

      ! Keep searching until all paired loops are found (ie. the number
      ! of creation or annihilation operators becomes 0)

      ! Set letter shift values for contraction indices
      do i = 1, 4
         c_shift(i) = e1(i,1) + e1(i,2) + e2(i,1) + e2(i,2)
      end do

      do while (cops(1) /= 0 .and. cops(2) /= 0)
         ! Need to start the search with the largest number of operators;
         ! check to start the search with the creation or annihilation
         ! operators first.
         ! i1 and i2 allow us to index the creation or annihilation
         ! operators
         if (cops(1) >= cops(2)) then
            ! Loop through creation first
            i1 = 1
            i2 = 2
         else
            ! Loop through annihilation first
            i1 = 2
            i2 = 1
         end if

         found = .false.

         ! Loop over P/H/V/X
         do i = 1, 4
            do j = 1, c(i,i1)
               ii = 1+(4*(i-1)) + c_shift(i)

               ! Assign index letter to pair list
               p_list%plist(sp)%pindex(i1) = ind(ii)

               ! Mark which operator this index belongs to
               ! Contraction index always w.r.t the first operator
               p_list%plist(sp)%ops(i1) = 1

               ! Increase index letter shift
               c_shift(i) = c_shift(i) + 1

               ! Decrease number of annihilation or creation operators
               c(i,i1) = c(i,i1) - 1

               ! Look for matching operator
               do k = 1, 4
                  do l = 1, c(k,i2)
                     ii = 1+(4*(k-1)) + c_shift(k)

                     p_list%plist(sp)%pindex(i2) = ind(ii)
                     p_list%plist(sp)%linked = .false.
                     ! Ultimately contraction indices belong on both
                     ! tensors, but marking these for both the first
                     ! and second tensor helps to pick them out in the
                     ! code below
                     p_list%plist(sp)%ops(i2) = 2
                     c_shift(k) = c_shift(k) + 1
                     c(k,i2) = c(k,i2) - 1

                     ! Found a pair, so increment pair list index
                     sp = sp + 1
                     found = .true.
                     exit
                  end do
                  ! A pair has been found, so exit
                  if (found) exit
               end do
               if (found) exit
            end do
            ! Can't pair creation or annihilation so exit
            if (found) exit
         end do
         ! Check how many operators are left to assign
         cops = sum(c,dim=1)
      end do

      ! Set number of contraction loops
      nloop = sp-1

      ! Match external pairs, either to external ops on the same
      ! operator, or to external ops on the second operator
      t_shift = 0
      do while (sum(e1ops) /= 0 .or. sum(e2ops) /= 0)
         if (sum(e1ops) /= 0) then
            ! Search on first tensor
            call find_pairs(p_list, sp, 1, e1, e2, c, t_shift, c_shift,
     &                      e1ops, e2ops, item)

         else if (sum(e2ops) /= 0)then
            ! Search on second tensor
            ! Note the number of creation/annihilation ops has been
            ! switched
            call find_pairs(p_list, sp, 2, e2, e1, c, t_shift, c_shift,
     &                      e2ops, e1ops, item)
         end if

         ! Check for more operators
         e1ops = sum(e1, dim=1)
         e2ops = sum(e2, dim=1)
      end do


      ! We now have list of pairs and which ops they belong to + any
      ! contraction indices linking external indices on different
      ! operators.
      ! Now they must be assigned to an ITF index string for each tensor,
      ! in the correct positions

      ! Create a pair list for T1, T2 and the result tensor. This will
      ! allow manipulation of indices on different tensor without
      ! interfering with each other
      call make_pair_list(p_list, t1_list, 1, sp-1)
      call make_pair_list(p_list, t2_list, 2, sp-1)

      ! Create pair list for result tensor, only need external index
      shift = 1
      do i = nloop+1, sp-1
         r_list%plist(shift) = p_list%plist(i)
         shift = shift + 1
      end do


      if (item%permute == 2) then
         ! Need to swap annihilation operators between tensors:
         ! T1_{ac}^{ik} T2_{cb}^{kj} -> T1_{ac}^{jk} T2_{cb}^{ki}

         t1_list%plist(nloop+1)%pindex(2) =
     &                               p_list%plist(nloop+2)%pindex(2)
         t2_list%plist(nloop+1)%pindex(2) =
     &                               p_list%plist(nloop+1)%pindex(2)
         ! Final permutations are made in command_to_itf
      end if


      ! Swap between creation and annihilation operators to make sure
      ! external ops are before internal. This is important for the
      ! amplitudes and fock tensors, as these require a specific order of
      ! indices. Integrals and intermediates are ignored by this
      ! subroutine
      ! TODO: will not work with pqrstu...
      call swap_index(t1_list, item%rank1, item%int(1), item%inter(1))
      call swap_index(t2_list, item%rank2, item%int(2), item%inter(2))
      call swap_index(r_list, item%rank3, item%int(3), item%inter(3))


      ! Sort index pairs into order with bubble sort. If only 1 pair,
      ! then this is skipped. This is important to assure consistent use
      ! of indices for the declaration and use of an intermediate (the
      ! slot structure must be the same when it is constructed and when
      ! it is used)
      call swap_pairs(t1_list, item%rank1, item%int(1), e4)
      e4 = 0
      call swap_pairs(t2_list, item%rank2, item%int(2), e4)
      call swap_pairs(r_list, item%rank3, item%int(3), e5)


      ! Insert ordered lists into ITF index strings
      s1 = '        '
      s2 = '        '
      s3 = '        '

      do i = 1, item%rank1/2
         s1(i:i) = t1_list%plist(i)%pindex(1)
         s1(i+(item%rank1/2):i+(item%rank1/2)) =
     &                                        t1_list%plist(i)%pindex(2)
      end do
      do i = 1, item%rank2/2
         s2(i:i) = t2_list%plist(i)%pindex(1)
         s2(i+(item%rank2/2):i+(item%rank2/2)) =
     &                                        t2_list%plist(i)%pindex(2)
      end do
      do i = 1, item%rank3/2
         s3(i:i) = r_list%plist(i)%pindex(1)
         s3(i+(item%rank3/2):i+(item%rank3/2)) =
     &                                         r_list%plist(i)%pindex(2)
      end do

      ! Assign an ITF index string to each tensor
      item%idx1 = trim(s1)
      item%idx2 = trim(s2)
      item%idx3 = trim(s3)

      ! Release memory for pair lists
      deallocate(p_list%plist)
      deallocate(t1_list%plist)
      deallocate(t2_list%plist)
      deallocate(r_list%plist)

      return
      end


*----------------------------------------------------------------------*
      subroutine find_pairs(list, sp, tensor, e1, e2, c, t_shift,
     &                      c_shift, e1ops, e2ops, item)
*----------------------------------------------------------------------*
!     Find external pairs in a binary contraction. Assign a contraction
!     index if they are on different tensors.
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(pair_list), intent(inout) ::
     &   list              ! Complete pair list for binary contraction
      integer, intent(inout) ::
     &   sp               ! Pair list shift index
      integer, intent(inout) ::
     &   tensor,           ! Label T1 or T2
     &   e1(4,2),          ! External index occupations for a tensor
     &   e2(4,2),          ! External index occupations the other tensor
     &   c(4,2)            ! Contraction index occupations
      integer, dimension(3), intent(inout) ::
     &   t_shift,          ! External index letter shift
     &   c_shift           ! Creation index letter shift
      integer, dimension(2), intent(in) ::
     &   e1ops,
     &   e2ops
      type(itf_contr), intent(in) ::
     &   item              ! ITF binary contraction info

      character, dimension(16) ::
     &   ind=(/ 'a','b','c','d','i','j','k','l','p','q','r','s','t',
     &          'u','v','w' /)       ! Index letters
      logical ::
     &   found
      integer ::
     &   opp_tensor,       ! Label opposite tensor to 'tensor'
     &   i1, i2,           ! Label annihilation/creation index
     &   n1, n2,           ! Label occupations
     &   start,            ! Labels P/H/V/X to look for pairs
     &   ii,               ! Letter index
     &   i,j,k,l,m,n,      ! Loop indices
     &   nops(4)


      ! Decide which tensors we are dealing with (T1 or T2)
      if (tensor == 1) then
         opp_tensor = 2
      else
         opp_tensor = 1
      end if

      ! i1 and i2 label creation/annihilation
      ! operators, pindex always has a creation in position 1, so the
      ! indices must be placed in their correct position.
      if (e1ops(1) >= e1ops(2)) then
         ! Loop through creations first
         i1 = 1
         i2 = 2
         n1 = e1ops(2)
         n2 = e2ops(2)
      else
         ! Loop through annihilation first
         i1 = 2
         i2 = 1
         n1 = e1ops(1)
         n2 = e2ops(1)
      end if

      found = .false.

      ! Change i = 2, to start loop at hole index, create pair, then
      ! change back to i = 1 on next call, continue as normal till every
      ! index gone
      ! Have i = start, start set to be 1, 2, 3, 4 based on what the
      ! largest number of possible loops
      nops(1) = e1(1,i1) + e1(2,i2) + e1(3,i2) + e1(4,i2)
      nops(2) = e1(2,i1) + e1(1,i2) + e1(3,i2) + e1(4,i2)
      nops(3) = e1(3,i1) + e1(1,i2) + e1(2,i2) + e1(4,i2)
      nops(4) = e1(4,i1) + e1(1,i2) + e1(2,i2) + e1(3,i2)

      ! By default, always start pairing from particle index
      start = 1
      do i = 2, size(nops)
         if (nops(start) < nops(i)) then
            start = i
         end if
      end do

      !write(11,*) " NOPS ", nops
      !write(11,*) "start ", start


      ! Start main search loop
      !do i = 1, 4
      do i = start, 4
         do j = 1, e1(i,i1)
            ! Search for first operator
            ii = 1+(4*(i-1)) + t_shift(i)
            list%plist(sp)%pindex(i1) = ind(ii)
            list%plist(sp)%ops(i1) = tensor

            t_shift(i) = t_shift(i) + 1
            e1(i,i1) = e1(i,i1) - 1

            ! Look for matching operator
            if (n1 > 0) then
               ! Look on the same tensor
               do k = 1, 4
                  do l = 1, e1(k,i2)
                     ii = 1+(4*(k-1)) + t_shift(k)
                     list%plist(sp)%pindex(i2) = ind(ii)
                     list%plist(sp)%linked = .false.
                     list%plist(sp)%ops(i2) = tensor

                     t_shift(k) = t_shift(k) + 1
                     e1(k,i2) = e1(k,i2) - 1

                     ! Found a pair, so increment pair list index
                     sp = sp + 1
                     found = .true.
                     exit
                  end do
                  if (found) exit
               end do

            else if (n2 > 0) then
               ! Look on the opposite tensor
               do k = 1, 4
                  do l = 1, e2(k,i2)
                     ii = 1+(4*(k-1))+t_shift(k)
                     list%plist(sp)%pindex(i2) = ind(ii)
                     list%plist(sp)%ops(i2) = opp_tensor
                     list%plist(sp)%linked = .true.

                     t_shift(k) = t_shift(k) + 1
                     e2(k,i2) = e2(k,i2) - 1

                     ! We need to link the external indices on two
                     ! different operators with a contraction index
                     do m=1, 4
                        do n=1, c(m,i2)
                            ii = 1+(4*(m-1)) + c_shift(m)

                           list%plist(sp)%link = ind(ii)

                           c_shift(m) = c_shift(m) + 1
                           c(m,i2) = c(m,i2) - 1

                           sp = sp + 1
                           found = .true.
                           exit
                        end do
                        if (found) exit
                     end do
                     if (found) exit
                  end do
                  if (found) exit
               end do
            else
               ! No ops of opposite type to match...
               call line_error("Particle number not conserving",item)
            end if
            ! A pair loop has been found so exit
            if (found) exit
         end do
         if (found) exit
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine make_pair_list(p_list, t_list, place, npairs)
*----------------------------------------------------------------------*
!     Construct a pair list for a specific tensor from the pair list of
!     the complete binary contraction
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(pair_list), intent(in) ::
     &   p_list            ! Complete pair list for binary contraction
      type(pair_list), intent(inout) ::
     &   t_list            ! Pair list for specific tensor
      integer, intent(in) ::
     &   place,            ! Tensor place in binary contraction: {1,2,3}
     &   npairs            ! Number of pairs in binary contraction

      integer ::
     &   s,                ! Index to help placement of pairs
     &   i                 ! Loop index

      s = 1
      do i = 1, npairs
         if (p_list%plist(i)%ops(1) == place) then
            t_list%plist(s)%pindex(1) = p_list%plist(i)%pindex(1)

            if (p_list%plist(i)%linked) then
               t_list%plist(s)%pindex(2) = p_list%plist(i)%link
            else
               t_list%plist(s)%pindex(2) = p_list%plist(i)%pindex(2)
            end if

            t_list%plist(s)%linked = .false.
            t_list%plist(s)%ops = place
            s = s + 1
         else if (p_list%plist(i)%ops(2) == place) then
            t_list%plist(s)%pindex(2) = p_list%plist(i)%pindex(2)

            if (p_list%plist(i)%linked) then
               t_list%plist(s)%pindex(1) = p_list%plist(i)%link
            else
               t_list%plist(s)%pindex(1) = p_list%plist(i)%pindex(1)
            end if

            t_list%plist(s)%linked = .false.
            t_list%plist(s)%ops = place
            s = s + 1
         end if
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine swap_index(list, rank, integral, inter)
*----------------------------------------------------------------------*
!     Swap creation/annihilation operators around in pair index to make
!     sure that the external slots come before the internal. This is
!     needed because tensors need to be defined as T[eecc] and f[ec].
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(pair_list), intent(inout) ::
     &   list              ! Complete pair list for a tensor
      integer, intent(in) ::
     &   rank              ! Rank of tensor
      logical, intent(in) ::
     &   integral,
     &   inter

      character(len=1) ::
     &   tmp               ! Temporary holder for index
      integer ::
     &   i                 ! Loop index


      ! If the tensor is a rank-4 integral or an intermediate,
      ! creation/annihilation order is important, so these are skipped
      if (.not. integral .or. integral .and. rank == 2) then
         if (.not. inter) then
            do i = 1, rank/2
               if (list%plist(i)%pindex(1) >
     &                                     list%plist(i)%pindex(2)) then
                  tmp = list%plist(i)%pindex(1)
                  list%plist(i)%pindex(1) = list%plist(i)%pindex(2)
                  list%plist(i)%pindex(2) = tmp
               end if
            end do
         end if
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine swap_pairs(list, rank, integral, e4)
*----------------------------------------------------------------------*
!     Swap pairs of indices within a pair list using bubble sort.
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(pair_list), intent(inout) ::
     &   list              ! Complete pair list for a tensor
      integer, intent(in) ::
     &   rank              ! Rank of tensor
      logical, intent(in) ::
     &   integral
      integer, intent(in) ::
     &   e4(4,2)

      type(pair_list) ::
     &   tmp_list          ! Temporary holder for pair list
      logical ::
     &   sort              ! True is list has been sorted
      integer ::
     &   i                 ! Loop index

      allocate(tmp_list%plist(1))

      ! Skip if there is only one pair, or if this is an integral
      if (e4(1,2) > e4(1,1)) then
      if (rank > 2 .and. .not. integral) then
         sort = .true.
         do while (sort)
            sort = .false.
            do i = 1, rank/2
               if (i == rank/2) exit

               if (list%plist(i+1)%pindex(2) <
     &                                     list%plist(i)%pindex(2)) then
                  tmp_list%plist(1) = list%plist(i)
                  list%plist(i) = list%plist(i+1)
                  list%plist(i+1) = tmp_list%plist(1)
                  sort = .true.
               end if
            end do
         end do
      end if
      else
      ! The case where there is an intermediate with particle
      ! annihilation operators
      if (rank > 2 .and. .not. integral) then
         sort = .true.
         do while (sort)
            sort = .false.
            do i = 1, rank/2
               if (i == rank/2) exit

               if (list%plist(i+1)%pindex(1) <
     &                                     list%plist(i)%pindex(1)) then
                  tmp_list%plist(1) = list%plist(i)
                  list%plist(i) = list%plist(i+1)
                  list%plist(i+1) = tmp_list%plist(1)
                  sort = .true.
               end if
            end do
         end do
      end if
      end if

      ! If there are more particle annihilation operators, such as in an
      ! intermediate, then the pairs should be sorted according to these

      deallocate(tmp_list%plist)

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

      ! Check if not antisym over different vertices
      ! Check for tensor products
      ! Not going to antisymm intermediates...
      ! The intermediates can have eeaa structure
      e1=0
      e2=0
      c=0

      ! Get occupation info
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

      if (e1(1,1)+e2(1,1)==2 .and. e1(2,2)+e2(2,2)==2 .or.
     &    e1(3,1)+e2(3,1)==2 .and. e1(2,2)+e2(2,2)==2 .or.
     &    e1(1,1)+e2(1,1)==2 .and. e1(3,2)+e2(3,2)==2) then

         sum_c1=0
         sum_c2=0
         sum_a1=0
         sum_a2=0

         do i=1, 4
            ! Sum creation ops
            sum_c1=sum_c1+e1(i,1)
            sum_c2=sum_c2+e2(i,1)

            ! Sum annihilation ops
            sum_a1=sum_a1+e1(i,2)
            sum_a2=sum_a2+e2(i,2)
         end do

         ! If sum==2, then both indices come from same operator, therefore
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
         call spin_sum(s1a,s1b,s2a,s2b,t1a,t1b,t2a,t2b,zero_a,
     &                 zero_b,item)
      else
         call spin_sum(s1b,s1a,s2a,s2b,t1b,t1a,t2a,t2b,zero_b,
     &                 zero_a,item)
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine spin_sum(s1a,s1b,s2a,s2b,t1a,t1b,t2a,t2b,zero_a,
     &                    zero_b,item)
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
     &                                    s2b,eloop,item)
                    end do
                   end if
                  end do
                 else if (.not. any(s1b==0) .and. fourth_idx==0) then
                  ! Print result, sum over three indicies
                  call print_spin_case(s1b,s1a,s2a,s2b,
     &                                 eloop,item)
                 end if

                end do
               end if
              end do
             else if (.not. any(s1b==0) .and. third_idx==0) then
              ! Print result, sum over two indicies
              call print_spin_case(s1a,s1b,s2a,s2b,
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
             call print_spin_case(s1a,s1b,s2a,s2b,
     &                            eloop,item)
            end do
           end if
          end do
         else if (.not. any(s1a==0) .and. second_idx==0) then
          ! Print result, sum over one indicies
          call print_spin_case(s1a,s1b,s2a,s2b,
     &                         eloop,item)
         end if
        end do ! Loop over a/b for first index
       end if ! Check for first 0 index
      end do
      else if (.not. any(s1a==0) .and. first_idx==0) then
       ! Print result, sum over no index
       call print_spin_case(s1a,s1b,s2a,s2b,
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
      subroutine print_spin_case(s1a,s1b,s2a,s2b,eloop,item)
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
     &    sum(s2a)==sum(s2b)) then
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

            if (item%inter(1)) then
               if (item%swapped) then
                  call inter_spin_name(s2a,s2b,item%rank1/2,item%inter1)
               else
                  call inter_spin_name(s1a,s1b,item%rank1/2,item%inter1)
               end if
            end if

            if (item%inter(2)) then
               if (item%swapped) then
                  call inter_spin_name(s1a,s1b,item%rank2/2,item%inter2)
               else
                  call inter_spin_name(s2a,s2b,item%rank2/2,item%inter2)
               end if
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
                  !write(item%logfile,*) "intermediate", item%label_t1

                  ! For now, if intermediate is part of a permutation
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
      subroutine inter_spin_name(sa,sb,h_rank,label)
*----------------------------------------------------------------------*
!     Add spin name to intermediate (ie. STIN001 -> STIN001abab)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   sa(2),      ! Spin info
     &   sb(2),      ! Spin info
     &   h_rank      ! Rank/2
      character(len=index_len), intent(inout) ::
     &   label       ! Spin name of intermediate

      character(len=index_len) ::
     &   spin_name
      integer ::
     &   i

      spin_name = ''

      do i = 1, h_rank
         if (sa(i)==1) then
            spin_name(i:i) = 'a'
         else if (sa(i)==2) then
            spin_name(i:i) = 'b'
         end if

         if (sb(i)==1) then
            spin_name(i+h_rank:i+h_rank) = 'a'
         else if (sb(i)==2) then
            spin_name(i+h_rank:i+h_rank) = 'b'
         end if
      end do

      label = spin_name

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

      ! Assign index string. Tensor ranks are also set here
      if (comm==command_cp_intm .or. comm==command_add_intm) then
         ! For [ADD] and [COPY]
         call assign_add_index(contr_info,itf_item)
      else
         ! For other contractions
         call assign_index(contr_info,itf_item)
      end if

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
      real(8), intent(inout) ::
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
            fact = fact * -1.0d+0
      end if

      ! Catch three internal integrals
      if (c(2,1) + c(2,2)==3)
     &   then
            fact = fact * -1.0d+0
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
      write(item%logfile,*) 'END'

      return
      end


*----------------------------------------------------------------------*
      subroutine itf_rank(ops1, ops2, rank, flag)
*----------------------------------------------------------------------*
!     Calculate rank of tensor
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   ops1(4,2),
     &   ops2(4,2)
      integer, intent(inout) ::
     &   rank
      logical, intent(in) ::
     &   flag

      if (flag) then
         ! Only use the first occupation array to calculate the rank
         rank = sum(sum(ops1, dim=1))
      else
         rank = sum(sum(ops1, dim=1)) + sum(sum(ops2, dim=1))
      end if

      return
      end
