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

      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &    string
      integer, intent(in) ::
     &    rank          ! Rank of tensor
      character(len=MAXLEN_BC_LABEL) ::
     &    rename_tensor

      if (trim(string).eq.'O2g' .or. trim(string).eq.'O2') then
          rename_tensor='R'
      else if (trim(string).eq.'O1') then
          rename_tensor='R'
      else if (trim(string).eq.'T2g' .or. trim(string).eq.'T2') then
          rename_tensor='T'
      else if (trim(string).eq.'T1') then
          rename_tensor='T'
      else if (trim(string).eq.'T3') then
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

      character(len=MAXLEN_BC_LABEL), intent(in) ::
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

      character(len=MAXLEN_BC_LABEL), intent(in) ::
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
      recursive function c_index(idx, n) result(cidx)
*----------------------------------------------------------------------*
!     Cycle covarient tensor index: abcijk => cabijk
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   n        ! Number of cycles
      character(len=INDEX_LEN), intent(in) ::
     &   idx      ! ITF index string

      character(len=INDEX_LEN) ::
     &   cidx
      character(len=INDEX_LEN) ::
     &   tmp

      !TODO: only works for rank 6
      if (n == 0) then
         cidx = idx
      else
         tmp = idx
         tmp(1:1) = idx(3:3)
         tmp(2:2) = idx(1:1)
         tmp(3:3) = idx(2:2)
         cidx = c_index(tmp, n-1)
      end if

      end function c_index


*----------------------------------------------------------------------*
      pure function f_index(index, hrank, upper)
*----------------------------------------------------------------------*
!     Flip co/contravarient tensor index: abcijk => acbijk
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(*), intent(in) ::
     &   index          ! ITF index string
      integer, intent(in) ::
     &   hrank          ! Half rank
      logical, optional, intent(in) ::
     &   upper          ! Flip contravarient index
      character(hrank*2) ::
     &   f_index        ! Transpose of ITF index string

      logical ::
     &   contra

      if (present(upper)) then
         contra = upper
      else
         contra = .false.
      end if

      f_index=index

      if (hrank==1 .or. hrank==0) return

      if (contra) then
         f_index(hrank-1+hrank:hrank-1+hrank)=
     &                                index(hrank+hrank:hrank+hrank)
         f_index(hrank+hrank:hrank+hrank)=
     &                                index(hrank-1+hrank:hrank-1+hrank)
      else
         f_index(hrank-1:hrank-1)=index(hrank:hrank)
         f_index(hrank:hrank)=index(hrank-1:hrank-1)
      end if

      end function


      end module itf_utils


*----------------------------------------------------------------------*
      subroutine print_spin(spins, rank, label, logfile)
*----------------------------------------------------------------------*
!     Print spin of a tensor
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &   spins(2,*),
     &   rank,
     &   logfile
      character(len=*), intent(in) ::
     &   label

      integer ::
     &   i

      write(logfile,*) "SPIN OF ", trim(label)
      write(logfile,*) "================================="
      write(logfile,*) (spins(2,i), i=1, rank/2)
      write(logfile,*) (spins(1,i), i=1, rank/2)
      write(logfile,*) "================================="

      return
      end


*----------------------------------------------------------------------*
      subroutine print_plist(p_list, size, label, logfile)
*----------------------------------------------------------------------*
!     Print spin of a tensor
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(pair_list), intent(in) ::
     &   p_list
      integer, intent(in) ::
     &   size,
     &   logfile
      character(len=*), intent(in) ::
     &   label

      integer ::
     &   i, j

      write(logfile,*) "P_list: ", label
      write(logfile,*) "=============================="
      do i = 1, size
      write(logfile,*) p_list%plist(i)%pindex(1), p_list%plist(i)%ops(1)
      write(logfile,*) p_list%plist(i)%pindex(2), p_list%plist(i)%ops(2)
      write(logfile,*) p_list%plist(i)%link
      write(logfile,'(3i3)') (p_list%plist(i)%nval(j), j=1, 3)
      write(logfile,*) "---------------------------"
      end do
      write(logfile,*) "=============================="

      return
      end


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
      write(item%logfile,*) "Result: ", item%label_res, trim(item%idx3)
      write(item%logfile,*) "Tensor1: ", item%label_t1, trim(item%idx1)
      write(item%logfile,*) "Tensor2: ", item%label_t2, trim(item%idx2)
      write(item%logfile,*)

      return
      end


*----------------------------------------------------------------------*
      subroutine count_index2(idx,iocc,irstr,ngas,nspin,nops)
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

      integer :: i,j

      ! Creation operators
      ! Particle, ijkl
      nops(1,1)=nops(1,1) + iocc(1,1)
      ! Hole, abcd
      nops(2,1)=nops(2,1) + iocc(2,1)
      ! Valence, pqrs
      nops(3,1)=nops(3,1) + iocc(3,1)
      ! Explicit, x
      nops(4,1)=nops(4,1) + iocc(4,1)

      ! Annihilation operators as above
      nops(1,2)=nops(1,2) + iocc(1,2)
      nops(2,2)=nops(2,2) + iocc(2,2)
      nops(3,2)=nops(3,2) + iocc(3,2)
      nops(4,2)=nops(4,2) + iocc(4,2)

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

      integer :: i,j

!      ! Creation operators
!      ! Hole, abcd
!      nops(1,1)=nops(1,1) + iocc(2,1)
!      ! Particle, ijkl
!      nops(2,1)=nops(2,1) + iocc(1,1)
!      ! Valence, pqrs
!      nops(3,1)=nops(3,1) + iocc(3,1)
!      ! Explicit, x
!      nops(4,1)=nops(4,1) + iocc(4,1)
!
!      ! Annihilation operators as above
!      nops(1,2)=nops(1,2) + iocc(2,2)
!      nops(2,2)=nops(2,2) + iocc(1,2)
!      nops(3,2)=nops(3,2) + iocc(3,2)
!      nops(4,2)=nops(4,2) + iocc(4,2)

      ! Creation operators
      ! Hole, abcd
      nops(1,1)=nops(1,1) + iocc(2,1)
      ! Valence, pqrs
      nops(2,1)=nops(2,1) + iocc(3,1)
      ! Particle, ijkl
      nops(3,1)=nops(3,1) + iocc(1,1)
      ! Explicit, x
      nops(4,1)=nops(4,1) + iocc(4,1)

      ! Annihilation operators as above
      nops(1,2)=nops(1,2) + iocc(2,2)
      nops(2,2)=nops(2,2) + iocc(3,2)
      nops(3,2)=nops(3,2) + iocc(1,2)
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
     &     item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
!     &    perm_array(4),   ! Info of permutation factors
     &    perm_case,   ! Info of permutation factors
     &    i, j                ! Loop index
      logical ::
     &    inter,           ! True if result is an intermediate
     &    found
      character(len=MAXLEN_BC_LABEL) ::
     &    old_name,
     &    old_inter
      character(len=INDEX_LEN) ::
     &    old_idx

      ! Initialise permutation factors:
      ! 0 == no permutation
      ! 1 == (1-Pxy)
      ! 2 == (1-Pxy)(1-Pvw) = (1+Pxy) in spatial orbitals
      perm_case = 0
      do i = 1, ngastp
         if(contr_info%perm(i)) perm_case = perm_case + 1
      end do

      ! Determine if result needs permuting
      !inter = check_inter(contr_info%label_res)

      !if (.not.inter) then
      !   call permute_tensors(contr_info,perm_case,itflog)
      !end if

      ! If the perm_array doesn't contain any zeros, then we should
      ! introduce an intermediate which collects half of the different
      ! permutation cases, then do:
      ! .R[abij] += I[abij]
      ! .R[abij] += I[baji]
      ! Save old name and replace it with a new one
      if (perm_case > 0) then
         old_name = contr_info%label_res
         contr_info%label_res = "ITIN"
         item%symm = .true.
      end if

      ! Pick out specific commands, form the itf_contr object, spin sum
      ! and print out contraction line
      if (command==command_add_intm .or. command==command_cp_intm) then
         ! For [ADD] and [COPY] cases
         call itf_contr_init(contr_info,item,0,command,itflog)
         call print_itf_line(item,.false.,.false.)
      else
         ! For other binary contractions
         if (perm_case == 0) then
            ! No permutations
            call itf_contr_init(contr_info,item,0,command,itflog)
            call assign_spin(item)
         else
            do i=1, perm_case
               ! Loop over permutation cases and send separately to
               ! assign_spin. For most cases this is just one, however
               ! for (1-Pij)(1-Pab), we need to generate one of these
               ! permutations before symmetrising
               call itf_contr_init(contr_info,item,i,command,itflog)

               if (i == 2) then
                  ! Need to transpose by tensors after permutation, to
                  ! avoid symmetry problem when using (1 + Pabij)
                  ! If the intermediate has three internal/external
                  ! indicies, then permute the covarient index
                  if (item%inter(1)) then
                     found = .false.
                     do j = 1, ngastp
                        if (item%nops1(j) > 2) then
                        item%idx1=f_index(item%idx1,item%rank1/2,.true.)
                        found = .true.
                        exit
                        end if
                     end do
                     if (.not.found) then
                        item%idx1=f_index(item%idx1,item%rank1/2)
                     end if
                  else
                     item%idx1=f_index(item%idx1,item%rank1/2)
                  end if
                  item%idx2=f_index(item%idx2,item%rank2/2)

                  ! Whenever we tranpose a tensor, we intoroduce a sign
                  ! chage
                  ! No sign due to the tranpose of idx1 which defines an
                  ! intermeidate. Extra signs to to tranpose of tensors
                  ! which define intermediates are included in the
                  ! intermediate line
                  if (item%permute==2) then
                     if (item%rank2>2) item%fact = item%fact * -1.0d+0
                  end if
                  !write(11,*) "index flip fact: ", item%fact

                  ! Sometimes, the permuation intermediate goes into the
                  ! same intermediate as the non-permuation one.
                  ! This means there will be a factor of 2, so need to
                  ! get rid of this
                  if (trim(old_inter)==trim(item%label_t1) .and.
     &                trim(old_idx)==trim(item%idx1)) then
                     exit
                  end if
               end if

               old_inter = trim(item%label_t1)
               old_idx = trim(item%idx1)
               call assign_spin(item)
            end do

            ! If created a perm intermediate, print the symmetrised lines
            call print_symmetrise(old_name,item)
         end if
      end if

      ! Deallocate memroy used when construcitng item
      call itf_deinit(item)

      return
      end

*----------------------------------------------------------------------*
      subroutine intermediate_spin_info(contr_info,itflog,command,
     &                               spin_inters,n_inter,permute)
*----------------------------------------------------------------------*
!     Search for intermediates, doesn't print out line
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
     &    i,j,k
      logical ::
     &    summed


      call itf_contr_init(contr_info,itf_item,permute,
     &                    command,itflog)

      if (permute == 2) then
         ! Need to transpose by tensors after permutation, to
         ! avoid symmetry problem when using (1 + Pabij)
         ! When we transpose tensors, we get a sign change, however
         ! here, we are just searching for intermediates and don't care
         ! about the final sign - that will come when the intermediate
         ! is printed out
         itf_item%idx1 = f_index(itf_item%idx1,itf_item%rank1/2)
         itf_item%idx2 = f_index(itf_item%idx2,itf_item%rank2/2)
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

                     ! TODO: This code is repeated in print_itf.f
                     do j = 1, itf_item%rank3/2
                        itf_item%i_spin%spin(1,j) =
     &                                       spin_inters(i)%cases(j,k)
                        itf_item%i_spin%spin(2,j) =
     &                             spin_inters(i)%cases(j+INDEX_LEN/2,k)
                     end do

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

      ! TODO: Ignore tensor product cases for now...
      if (itf_item%rank3 /= 4 .and. itf_item%rank1 /= 2 .and.
     &    itf_item%rank2 /= 2) then
         if (itf_item%ninter == 0) call line_error("Couldn't find
     &                                    intermediate", itf_item)
      end if


      ! Copy information back to array in print_itf()
      do i = 1, itf_item%ninter
         spin_inters(i+n_inter)=itf_item%inter_spins(i)
         !write(itf_item%logfile,*) "INER SPIN: ",itf_item%inter_spins(i)
      end do

      ! Overall number of intermediates used to index spin_inters
      n_inter = n_inter + itf_item%ninter

      deallocate(itf_item%inter_spins)

      ! Deallocate memroy used when construcitng itf_item
      call itf_deinit(itf_item)

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
      character(len=MAXLEN_BC_LABEL) ::
     &    old_name


      ! Initialise permutation factors to 0 == no permutation
      perm_case = 0
      do i = 1, ngastp
         if(contr_info%perm(i)) perm_case = perm_case + 1
      end do

      if (perm_case == 0) then
         ! No permutations
         call intermediate_spin_info(contr_info,itflog,command,
     &                            spin_inters,n_inter,1)
      else
         do i=1, perm_case
            call intermediate_spin_info(contr_info,itflog,command,
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
!               item%inter_spins(j)%cases(:,j) = (/ 1, 0, 1, 0 /)
               item%inter_spins(j)%cases(:,j) = 0
               item%inter_spins(j)%cases(1,j) = 1
               item%inter_spins(j)%cases(1+INDEX_LEN/2,j) = 1
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
               !item%inter_spins(j)%cases(:,j) = (/ 0, 0, 0, 0 /)
               item%inter_spins(j)%cases(:,j) = 0
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
      subroutine intermediate_to_itf(contr_info,itflog,command,
     &                               label,spin_case)
*----------------------------------------------------------------------*
!     Actually print out the intermediate lines
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
     &     spin_case(INDEX_LEN)
      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &     label

      type(itf_contr) ::
     &     item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &    i, j, k, l
      character(len=INDEX_LEN) ::
     &    spin_name
      logical ::
     &   found


      call itf_contr_init(contr_info,item,0,
     &                    command,itflog)

      ! Set overall spin case of result
      !write(item%logfile,*) "spin_case ", spin_case
      do i = 1, item%rank3/2
         item%i_spin%spin(1,i) = spin_case(i)
         item%i_spin%spin(2,i) = spin_case(i+INDEX_LEN/2)
      end do

      ! Change intermediate name to reflect spin case
      spin_name = ''
      j = 1
      do k = 1, 2
         do i = 1, item%rank3/2
            if (item%i_spin%spin(k,i)==1) then
               spin_name(j:j) = 'a'
               j = j + 1
            else if (item%i_spin%spin(k,i)==2) then
               spin_name(j:j) = 'b'
               j = j + 1
            end if
         end do
      end do

      ! If an intermediate arises as the result of a permutation, we
      ! need to create this new intermediate. This requires the
      ! transpose
      ! TODO: use a flag in itf_contr
      if (scan('P', label)) then

         found = .false.
         do j = 1, ngastp
            if (item%nops3(j) > 2) then
               ! If there are 3 external/internal indcies, don't need to
               ! permute index, but do need to permute spin_name
               do k = 1, len(spin_name)
                  if (spin_name(k:k)=='a') then
                     spin_name(k:k)='b'
                  else if (spin_name(k:k)=='b') then
                     spin_name(k:k)='a'
                  end if
               end do
               item%idx3=f_index(item%idx3,item%rank3/2,.true.)
               found = .true.
               exit
            end if
         end do

         do j = 1, ngastp
            ! Need to catach three internal integrals which result from a
            ! permuation and mark them. process.py will turn KP -> J
            if (item%nops1(j) > 2) then
               item%label_t1 = 'KP'
            end if
         end do

         if (.not.found) then
            item%idx3=f_index(item%idx3,item%rank3/2)
         end if

         item%idx1 = f_index(item%idx1,item%rank1/2,.true.)
         item%idx2 = f_index(item%idx2, item%rank2/2)


         ! This factor is from (1-Pij)
         item%fact = item%fact * -1.0d+0

         ! These factors come from permutations in order to match up
         ! indcices in final intermediate
         if (item%rank1>2) item%fact = item%fact * -1.0d+0
         if (item%rank2>2) item%fact = item%fact * -1.0d+0


      end if

      item%label_res = trim(item%label_res)//trim(spin_name)
      call assign_spin(item)

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

      character(len=MAXLEN_BC_LABEL) ::
     &     nres, nt1, nt2          ! Name of tensors involved in the contraction
      character(len=5) ::
     &     s_int                 ! Intermdiate tensor number
      character(len=264) ::
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

      ! Add intermediate spin strings to names
      if (item%inter(1)) nt1 = trim(nt1)//trim(item%inter1)
      if (item%inter(2)) nt2 = trim(nt2)//trim(item%inter2)

      ! Change tensor to spatial orbital quantity, unless it is an
      ! intermediate
      call spatial_string(st1,item%idx1,nt1,s1,item%inter(1),item%rank1,
     &                    1,item%binary,item%three(1),item%logfile)
      call spatial_string(st2,item%idx2,nt2,s2,item%inter(2),item%rank2,
     &                    2,item%binary,item%three(2),item%logfile)

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
      write(item%logfile,'(a)') trim(itf_line)

      return
      end


*----------------------------------------------------------------------*
      subroutine spatial_string(st,idx,nt,spin,inter,rank,tensor,binary,
     &                          three,lulog)
*----------------------------------------------------------------------*
!     Construct spatial tensor representation
*----------------------------------------------------------------------*

      use itf_utils

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=264), intent(inout) ::
     &   st          ! Name of spin summed tensors + index
      character(len=INDEX_LEN), intent(in) ::
     &   idx         ! Index of tensor
      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &   nt          ! Name of tensors involved in the contraction
      logical, intent(in) ::
     &   spin,       ! True if pure spin
     &   inter,      ! True if an intermediate
     &   binary,     ! True if a binary contraction
     &   three       ! True if three operators of the same type
      integer, intent(in) ::
     &   rank,       ! Rank of tensor
     &   tensor,     ! T1 or T2
     &   lulog       ! Logfile

      integer ::
     &   hrank       ! Half rank

      hrank = rank / 2

      if (spin .and. .not.inter) then
         ! Pure spin
         select case (rank)
            case (4)
               st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &            trimal(nt)//'['//f_index(idx,hrank)//']'//')'
            case (6)
               st='('//trimal(nt)//'['//trim(idx)//']'//' + '//
     &            trimal(nt)//'['//trim(c_index(idx,1))//']'//' + '//
     &            trimal(nt)//'['//trim(c_index(idx,2))//']'//' - '//
     &            trimal(nt)//'['//f_index(idx,hrank)//']'//' - '//
     &            trimal(nt)//'['//
     &            f_index(c_index(idx,1),hrank)//']'//' - '//
     &            trimal(nt)//'['//
     &            f_index(c_index(idx,2),hrank)//']'//')'
            case default
               write(lulog,*) "ERROR: Couldn't determine rank"
         end select
      else
         if (tensor==2 .and. .not.binary) then
            ! Don't need second operator for [ADD] or [COPY]
            st=''
            return
         else
         select case (rank)
            case (0)
               st=trimal(nt)//'['//trim(idx)//']'
            case (2)
               st=trimal(nt)//'['//trim(idx)//']'
            case (4)
               st=trimal(nt)//'['//trim(idx)//']'
            case (6)
               st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &              trimal(nt)//'['//f_index(idx,hrank)//']'//')'
            case default
               write(lulog,*) "ERROR: Couldn't determine rank"
         end select
         end if
      end if

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
      character(len=INDEX_LEN) ::
     &     c1, c2, c3,
     &     a1, a2, a3
      character(len=INDEX_LEN), dimension(8) ::
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
!      do i=1, o1(1,1)
!          c1(i:)=par(i)
!      end do
!      o1_array(1)=c1
!      do i=1, o1(3,1)
!          c2(i:)=val(i)
!      end do
!      o1_array(2)=c2
!      do i=1, o1(2,1)
!          c3(i:)=hol(i)
!      end do
!      o1_array(3)=c3

      do i=1, o1(1,1)
          c1(i:)=par(i)
      end do
      o1_array(1)=c1
      do i=1, o1(2,1)
          c2(i:)=val(i)
      end do
      o1_array(2)=c2
      do i=1, o1(3,1)
          c3(i:)=hol(i)
      end do
      o1_array(3)=c3

      ! Need to to be shifted to not match assignment of creations above
!      do i=1, o1(1,2)
!          a1(i:)=par(i+o1(1,1))
!      end do
!      o1_array(5)=a1
!      do i=1, o1(3,2)
!          a2(i:)=val(i+o1(3,1))
!      end do
!      o1_array(6)=a2
!      do i=1, o1(2,2)
!          a3(i:)=hol(i+o1(2,1))
!      end do
!      o1_array(7)=a3


      do i=1, o1(1,2)
          a1(i:)=par(i+o1(1,1))
      end do
      o1_array(5)=a1
      do i=1, o1(2,2)
          a2(i:)=val(i+o1(2,1))
      end do
      o1_array(6)=a2
      do i=1, o1(3,2)
          a3(i:)=hol(i+o1(3,1))
      end do
      o1_array(7)=a3

      item%idx1=trimal(o1_array(1))//trimal(o1_array(2))//
     &          trimal(o1_array(3))//trimal(o1_array(5))//
     &          trimal(o1_array(6))//trimal(o1_array(7))

      item%idx3=item%idx1

      return
      end


*----------------------------------------------------------------------*
      subroutine assign_nval(ncase, index, np, p_list)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   ncase,
     &   index,
     &   np
      type(pair_list), intent(inout) ::
     &   p_list

      ! Assign 'value' to index
      select case (ncase)
         case (1)
            p_list%plist(np)%nval(index) = 1
         case (2)
            p_list%plist(np)%nval(index) = 2
         case (3)
            p_list%plist(np)%nval(index) = 3
         case (4)
            p_list%plist(np)%nval(index) = 4
      end select

      return
      end


*----------------------------------------------------------------------*
      subroutine create_index_str(idx, cnt, ex, c_shift, e_shift, rank)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(inout) ::
     &   idx
      integer, intent(in) ::
     &   cnt(4,2),        ! Operator numbers of contraction index
     &   ex(4,2),
     &   c_shift(4),
     &   rank
      integer, intent(inout) ::
     &   e_shift(4)

      integer ::
     &   shift_c,
     &   shift_a,
     &   cnt_shift(ngastp),
     &   i, j, k, l,
     &   shift_p
      character, dimension(21) ::
     &   ind=(/ 'i','j','k','l','m','n','o','a','b','c','d','e','f',
     &          'g','p','q','r','s','t','u','v' /)   ! Letters for index string
      type(pair_list) ::
     &   p_list

      allocate(p_list%plist(rank/2))

      shift_c = 1
      shift_a = rank
      cnt_shift = c_shift
      l = 1

      do i = 1, ngastp
         if (cnt(i,1)>0) then
            do j = 1, cnt(i,1)
               k = 1+(7*(i-1)) + cnt_shift(i)
               idx%str(shift_c) = ind(k)
               idx%cnt_poss(l) = shift_c
               l = l + 1
               cnt_shift(i) = cnt_shift(i) + 1
               shift_c = shift_c + 1
            end do
         end if
         if (ex(i,1)>0) then
            do j = 1, ex(i,1)
               k = 1+(7*(i-1)) + e_shift(i)
               idx%str(shift_c) = ind(k)
               e_shift(i) = e_shift(i) + 1
               shift_c = shift_c + 1
            end do
         end if

         if (cnt(i,2)>0) then
            do j = 1, cnt(i,2)
               k = 1+(7*(i-1)) + cnt_shift(i)
               idx%str(shift_a) = ind(k)
               idx%cnt_poss(l) = shift_a
               l = l + 1
               cnt_shift(i) = cnt_shift(i) + 1
               shift_a = shift_a - 1
            end do
         end if
         if (ex(i,2)>0) then
            do j = 1, ex(i,2)
               k = 1+(7*(i-1)) + e_shift(i)
               idx%str(shift_a) = ind(k)
               e_shift(i) = e_shift(i) + 1
               shift_a = shift_a - 1
            end do
         end if
      end do

      !call print_plist(p_list, rank/2, "P list", 11)

      deallocate(p_list%plist)

      return
      end


*----------------------------------------------------------------------*
      subroutine assign_new_index(contr_info,item)
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
     &   ci(4,2),        ! Operator numbers of contraction index
     &   e1(4,2),       ! Operator numbers of external index 1
     &   e2(4,2),       ! Operator numbers of external index 2
     &   e3(4,2),       ! Operator numbers of result index
     &   e4(4,2),       ! Operator numbers of result index
     &   e5(4,2),       ! Operator numbers of result index
     &   te1(4,2),      ! Test Operator numbers of external index 1
     &   te2(4,2),      ! Test Operator numbers of external index 2
     &   tc(4,2),       ! Test Operator numbers of external index 2
     &   i, j, k, l, m,    ! Loop index
     &   ii,            ! Letter index
     &   nloop,         ! Number of contraction loops
     &   i1, i2,        ! Search creation or annihilation first
     &   shift,         ! List shift
     &   shift_a,         ! List shift
     &   shift_c,        ! List shift
     &   sp             ! Pair list shift
      character, dimension(21) ::
     &   ind=(/ 'a','b','c','d','e','f','g','p','q','r','s','t','u',
     &          'v','i','j','k','l','m','n','o' /)   ! Letters for index string
      character(len=INDEX_LEN) ::
     &     s1, s2, s3   ! Tmp ITF index strings
      character(len=1) ::
     &   tmp,
     &   tmp2
      real(8) ::
     &   factor         ! Factor from equivalent lines
      logical ::
     &   found,         ! True if found pairing index
     &   sort           ! Used in bubble sort
      type(pair_list) ::
     &   p_list,        ! Complete list of pairs in binary contraction
     &   p_list2,        ! Complete list of pairs in binary contraction
     &   t1_list,       ! List of pairs in first tensor (T1)
     &   t2_list,       ! List of pairs in second tensor (T2)
     &   r_list,        ! List of pairs in result tensor
     &   d_list         ! debug list
      integer, dimension(4) ::
     &   t_shift,       ! Index shift for external indices
     &   e_shift,       ! Index shift for external indices
     &   c_shift        ! Index shift for contraction indices
      integer, dimension(2) ::
     &   cops,          ! Number of creation/annihilation contraction operators
     &   e1ops,         ! Number of C/A external T1 operators
     &   e2ops,         ! Number of C/A external T2 operators
     &   e3ops          ! Total number of external operators
      integer, dimension(4) ::
     &   tops           ! Scratch space
!      character, pointer ::
!     &     str1(:) => null(),
!     &     str2(:) => null(),
!     &     str3(:) => null()
      type(index_str) ::
     &   str1,
     &   str2,
     &   str3,
     &   t_str1,
     &   t_str2
      integer :: n_cnt, s, rank, rank1, tensor
      logical :: is_cnt, found_end, found_match, found_cnt

      c=0
      e1=0
      e2=0
      e3=0

      ! Get occupation info
      do i = 1, contr_info%n_cnt
        call count_index2(i,
     &     contr_info%occ_cnt(1:,1:,i),
     &     contr_info%rst_cnt(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,c)
      end do
      do i = 1, contr_info%nj_op1
        call count_index2(i,
     &     contr_info%occ_ex1(1:,1:,i),
     &     contr_info%rst_ex1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e1)
      end do
      do i = 1, contr_info%nj_op2
        call count_index2(i,
     &     contr_info%occ_ex2(1:,1:,i),
     &     contr_info%rst_ex2(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e2)
      end do
      do i = 1, contr_info%nj_res
        call count_index2(i,
     &     contr_info%occ_res(1:,1:,i),
     &     contr_info%rst_res(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e3)
      end do

!      ! Figure out factor from equivalent lines
!      factor = 1.0d+0
!      do j = 1, 2
!         do i = 1, 4
!            if (c(i,j) == 0) cycle
!            if (mod(c(i,j),2) == 0) then
!               factor = factor * (1.0d+0/real(c(i,j),8))
!            end if
!         end do
!      end do
!      item%fact = item%fact * factor
!
      ! Set ranks of tensors
      call itf_rank(e1, c, item%rank1, .false.)
      call itf_rank(e2, c, item%rank2, .false.)
      call itf_rank(e3, c, item%rank3, .true.)
!
!      ! Set number of indcies
!      item%nops1 = sum(e1, dim=2) + sum(c, dim=2)
!      item%nops2 = sum(e2, dim=2) + sum(c, dim=2)
!      item%nops3 = sum(e3, dim=2)

      allocate(str1%str(item%rank1))
      allocate(str2%str(item%rank2))
      allocate(str3%str(item%rank1+item%rank2))

      n_cnt = sum(sum(c, dim=1))
      allocate(str1%cnt_poss(n_cnt))
      allocate(str2%cnt_poss(n_cnt))
      allocate(str3%cnt_poss(n_cnt))

      allocate(str1%fact(n_cnt))
      allocate(str2%fact(n_cnt))
      allocate(str3%fact(n_cnt))

      allocate(p_list%plist(item%rank1/2+item%rank2/2))
      allocate(p_list2%plist(item%rank3/2))

      ! Set letter shift values for contraction indices
      do i = 1, 4
         c_shift(i) = e1(i,1) + e1(i,2) + e2(i,1) + e2(i,2)
      end do

      ! Make 'tranpose' of c array
      do i = 1, 2
         do j = 1, ngastp
            if (i==1) then
               k = 2
            else
               k = 1
            end if
            ci(j,i) = c(j,k)
         end do
      end do

      ! e_shift updates after each call
      e_shift = 0
      call create_index_str(str1,c,e1, c_shift, e_shift, item%rank1)
      call create_index_str(str2,ci,e2,c_shift, e_shift, item%rank2)

      do i = 1, item%rank1
         str3%str(i) = str1%str(i)
      end do
      do i = 1, item%rank2
         str3%str(i+item%rank1) = str2%str(i)
      end do

      write(item%logfile,*) "STR1: {", str1%str, "}"
      write(item%logfile,*) "STR2: {", str2%str, "}"
      write(item%logfile,*) "STR3: {", str1%str, "}{", str2%str, "}"

      write(item%logfile,*) "CNT POSS1: ", str1%cnt_poss
      write(item%logfile,*) "CNT POSS2: ", str2%cnt_poss

      do i = 1, n_cnt
       if (mod(item%rank1-str1%cnt_poss(i)+str2%cnt_poss(i)-1,2)==0)then
          str3%fact(i) = 1
       else
          str3%fact(i) = -1
       end if
       write(11,*) "FACTOR: ", str3%fact(i)
      end do


      ! Loop to find external pairs
      ! TODO: Start with most external indices
      if (sum(sum(e1,dim=2))>=sum(sum(e2,dim=2))) then
         t_str1 = str1
         rank1 = item%rank1
      else
         t_str1 = str2
         rank1 = item%rank2
      end if

      write(item%logfile,*)
      write(item%logfile,*) "====================================="
      write(item%logfile,*) "Looking for pairs in assign_new_index"
      write(item%logfile,*) "====================================="
      shift = 1
      do i = 1, rank1

         ! Find first none contraction index
         is_cnt = .false.
         do j = 1, n_cnt
            if (i==t_str1%cnt_poss(j)) then
               is_cnt = .true.
            end if
         end do

         tmp = t_str1%str(rank1-i+1)

         ! Need to check if previous pair has already been formed
         ! then continue
         if (shift>1 .and. .not. is_cnt) then
            do l = 1, shift-1
               do m = 1, 2
                  if (p_list2%plist(l)%pindex(m) == t_str1%str(i)) then
                     ! This pair has already been found
                     ! Mark is_cnt as true so as to skip to the next index
                     write(item%logfile,*) "Already found this index: ",
     &                           t_str1%str(i)
                     is_cnt = .true.
                     exit
                  end if
               end do
            end do
         end if


         if (.not. is_cnt) then

            write(item%logfile,*) "--------------------------------"
            write(item%logfile,*) "first ex index ", t_str1%str(i)
            write(item%logfile,*) "looking for ", tmp
            write(item%logfile,*) "--------------------------------"

            found_end = .false.
            tensor = 2
            do while (.not. found_end)

               !found_end = .false.
               ! Set up variables
               if (tensor == 2) then
                  rank = item%rank2
                  t_str2= str2
               else
                  rank = item%rank1
                  t_str2 = str1
               end if

               found_match = .false.
               do j = 1, rank
                  write(item%logfile,*) "Searching for match ",
     &                                  t_str2%str(j)," ",tmp
                  if (t_str2%str(j) == tmp) then
                     found_match = .true.
                     tmp2 = t_str2%str(rank-j+1)
                     write(item%logfile,*) t_str2%str(j),
     &                                     " is paired with ", tmp2

                     found_cnt = .false.
                     do k = 1, n_cnt
                        ! Is the pair a cnt index?
                        if ((rank-j+1) == t_str2%cnt_poss(k)) then
                           found_end = .false.
                           if (tensor == 2) then
                              tensor = 1
                           else
                              tensor = 2
                           end if
                           tmp = tmp2
                           found_cnt = .true.
                           exit
                        end if
                     end do

                     if (.not. found_cnt) then
                        ! If not, end
                        found_end = .true.
                        exit
                     end if
                     if (found_cnt) exit

                  end if
                  if (found_end) exit
               end do
               ! The index is not on the other tensor (two
               ! external indices on one tensor)
               if (.not. found_match) then
                  tmp2 = tmp
                  found_end = .true.
                  exit
               end if

            end do

            do l = 1, rank1/2
               if (t_str1%str(i)==t_str1%str(l)) then
                  ! Started with a creation operator
                  p_list2%plist(shift)%pindex(1) = t_str1%str(i)
                  p_list2%plist(shift)%pindex(2) = tmp2
                  exit
               else if (t_str1%str(i)==t_str1%str(l+rank1/2)) then
                  ! Started with an annhilation operator
                  p_list2%plist(shift)%pindex(1) = tmp2
                  p_list2%plist(shift)%pindex(2) = t_str1%str(i)
                  exit
               end if
            end do

            write(item%logfile,*) "================================"
            write(item%logfile,*) "Found an external pair:"
            write(item%logfile,*) "================================"
            write(item%logfile,*) "creation index    ",
     &                            p_list2%plist(shift)%pindex(1)
            write(item%logfile,*) "annhilation index ",
     &                            p_list2%plist(shift)%pindex(2)
            write(item%logfile,*) "================================"

            shift = shift + 1

            !TODO: whole thing in a loop, decrease number of externals
            !untill found all on both tensor (for example there are two
            !on T1 and two on T2

         end if
      end do

      write(item%logfile,*) "====================================="
      write(item%logfile,*) "Ending external pair search"
      write(item%logfile,*) "====================================="
      write(item%logfile,*)

!      ! Create pairs in order to get result index
!      do i = 1, item%rank1/2
!         p_list%plist(i)%pindex(1) = str1%str(i)
!         p_list%plist(i)%pindex(2) = str1%str(item%rank1-(i-1))
!         do j = 1, n_cnt
!            if (i==str1%cnt_poss(j)) p_list%plist(i)%linked = .true.
!         end do
!      end do
!      do i = 1, item%rank2/2
!      p_list%plist(i+item%rank1/2)%pindex(1)= str2%str(i)
!      p_list%plist(i+item%rank1/2)%pindex(2)= str2%str(item%rank2-(i-1))
!      do j = 1, n_cnt
!      if (i==str2%cnt_poss(j)) then
!         p_list%plist(i+item%rank1/2)%linked = .true.
!      end if
!      end do
!      end do
!
!      call print_plist(p_list, item%rank1/2+item%rank2/2, "P list", 11)
!
!      s = 1
!      do i = 1, item%rank1/2 + item%rank2/2
!         if (p_list%plist(i)%linked) write(11,*) "hello"
!         tmp = p_list%plist(i)%pindex(1)
!         do j = i+1, item%rank1/2 + item%rank2/2
!            if (tmp == p_list%plist(j)%pindex(2)) then
!               p_list2%plist(s)%pindex(1) = p_list%plist(j)%pindex(1)
!               p_list2%plist(s)%pindex(2) = p_list%plist(i)%pindex(2)
!               s = s + 1
!            end if
!         end do
!         tmp = p_list%plist(i)%pindex(2)
!         do j = i+1, item%rank1/2 + item%rank2/2
!            if (tmp == p_list%plist(j)%pindex(1)) then
!               p_list2%plist(s)%pindex(1) = p_list%plist(i)%pindex(1)
!               p_list2%plist(s)%pindex(2) = p_list%plist(j)%pindex(2)
!               s = s + 1
!            end if
!         end do
!      end do
!
!      call print_plist(p_list2, item%rank3/2, "P list2", 11)

      deallocate(p_list2%plist)
      deallocate(p_list%plist)

      deallocate(str3%fact)
      deallocate(str2%fact)
      deallocate(str1%fact)

      deallocate(str3%cnt_poss)
      deallocate(str2%cnt_poss)
      deallocate(str1%cnt_poss)

      deallocate(str3%str)
      deallocate(str2%str)
      deallocate(str1%str)

!      ! Find out number of creation/annihilation operators per operator
!      cops = sum(c, dim=1)
!      e1ops = sum(e1, dim=1)
!      e2ops = sum(e2, dim=1)
!      e3ops = sum(e3, dim=1)
!
!
!      if (item%rank1 == 0 .and. item%rank2 == 0) then
!         item%idx1 = ''
!         item%idx2 = ''
!         item%idx3 = ''
!         return
!      end if
!
!      ! Set number of contraction indicies, used later on
!      item%contri = sum(sum(c, dim=1))
!
!!      write(10,*) "e1 ", e1
!!      write(10,*) "e2 ", e2
!!      write(10,*) "e3 ", e3
!!      write(10,*) "c ", c
!!      write(10,*) "cops ", cops(1), cops(2)
!
!      ! Allocate pair list according to ranks of tensors
!      allocate(p_list%plist(item%rank1/2+item%rank2/2))
!      allocate(t1_list%plist(item%rank1/2))
!      allocate(t2_list%plist(item%rank2/2))
!      allocate(r_list%plist(item%rank3/2))
!
!
!      ! To start, we pair off contraction loops
!      ! Pair list index (shift pair)
!      sp = 1
!
!      ! Keep searching until all paired loops are found (ie. the number
!      ! of creation or annihilation operators becomes 0)
!
!      do while (cops(1) /= 0 .and. cops(2) /= 0)
!         ! Need to start the search with the largest number of operators;
!         ! check to start the search with the creation or annihilation
!         ! operators first.
!         ! i1 and i2 allow us to index the creation or annihilation
!         ! operators
!         if (cops(1) >= cops(2)) then
!            ! Loop through creation first
!            i1 = 1
!            i2 = 2
!         else
!            ! Loop through annihilation first
!            i1 = 2
!            i2 = 1
!         end if
!
!         found = .false.
!
!         ! Loop over P/H/V/X
!         do i = 1, 4
!            do j = 1, c(i,i1)
!               ii = 1+(7*(i-1)) + c_shift(i)
!
!               ! Assign index letter to pair list
!               p_list%plist(sp)%pindex(i1) = ind(ii)
!
!               ! Assign 'value' to index
!               call assign_nval(i, i1, sp, p_list)
!
!               ! Mark which operator this index belongs to
!               ! Contraction index always w.r.t the first operator
!               p_list%plist(sp)%ops(i1) = 1
!
!               ! Increase index letter shift
!               c_shift(i) = c_shift(i) + 1
!
!               ! Decrease number of annihilation or creation operators
!               c(i,i1) = c(i,i1) - 1
!
!               ! Look for matching operator
!               do k = 1, 4
!                  do l = 1, c(k,i2)
!                     ii = 1+(7*(k-1)) + c_shift(k)
!
!                     p_list%plist(sp)%pindex(i2) = ind(ii)
!                     call assign_nval(k, i2, sp, p_list)
!                     p_list%plist(sp)%linked = .false.
!                     ! Ultimately contraction indices belong on both
!                     ! tensors, but marking these for both the first
!                     ! and second tensor helps to pick them out in the
!                     ! code below
!                     p_list%plist(sp)%ops(i2) = 2
!                     c_shift(k) = c_shift(k) + 1
!                     c(k,i2) = c(k,i2) - 1
!
!                     ! Found a pair, so increment pair list index
!                     sp = sp + 1
!                     found = .true.
!                     exit
!                  end do
!                  ! A pair has been found, so exit
!                  if (found) exit
!               end do
!               if (found) exit
!            end do
!            ! Can't pair creation or annihilation so exit
!            if (found) exit
!         end do
!         ! Check how many operators are left to assign
!         cops = sum(c,dim=1)
!      end do
!
!      ! Set number of contraction loops
!      nloop = sp-1
!
!      ! Match external pairs, either to external ops on the same
!      ! operator, or to external ops on the second operator
!      t_shift = 0
!      do while (sum(e1ops) /= 0 .or. sum(e2ops) /= 0)
!         if (sum(e1ops) /= 0) then
!            ! Search on first tensor
!            call find_pairs(p_list, sp, 1, e1, e2, c, t_shift, c_shift,
!     &                      e1ops, e2ops, item)
!
!         else if (sum(e2ops) /= 0)then
!            ! Search on second tensor
!            ! Note the number of creation/annihilation ops has been
!            ! switched
!            call find_pairs(p_list, sp, 2, e2, e1, c, t_shift, c_shift,
!     &                      e2ops, e1ops, item)
!         end if
!
!         ! Check for more operators
!         e1ops = sum(e1, dim=1)
!         e2ops = sum(e2, dim=1)
!      end do
!
!      !call print_plist(p_list, sp-1, "P_LIST", item%logfile)
!
!      ! We now have list of pairs and which ops they belong to + any
!      ! contraction indices linking external indices on different
!      ! operators.
!      ! Now they must be assigned to an ITF index string for each tensor,
!      ! in the correct positions
!
!      ! Create a pair list for T1, T2 and the result tensor. This will
!      ! allow manipulation of indices on different tensor without
!      ! interfering with each other
!      call make_pair_list(p_list, t1_list, 1, sp-1)
!      call make_pair_list(p_list, t2_list, 2, sp-1)
!
!      !call print_plist(t1_list, item%rank1/2, "T1_LIST", item%logfile)
!      !call print_plist(t2_list, item%rank2/2, "T2_LIST", item%logfile)
!
!      ! Create pair list for result tensor, only need external index
!      shift = 1
!      do i = nloop+1, sp-1
!         r_list%plist(shift) = p_list%plist(i)
!         shift = shift + 1
!      end do
!
!      ! Only need to permute annihilation ops amongst themselves,
!      ! this is the case whenever we have (1-Pyx)(1-Pvw). For two rank
!      ! 4 tensors, this is straight forward. When there is a rank 2 and
!      ! rank 4, can't swap the rank 2, so have to swapped both
!      ! annihilations on the rank 4 tensor
!      if (item%permute == 2) then
!         ! Need to swap annihilation operators between tensors:
!         ! T1_{ac}^{ik} T2_{cb}^{kj} -> T1_{ac}^{jk} T2_{cb}^{ki}
!
!         if (item%rank1 /= 2) then
!            if (p_list%plist(nloop+2)%ops(2)==1) then
!               ! External indices belong on same tensor
!               tmp = p_list%plist(nloop+1)%pindex(2)
!               t1_list%plist(nloop+1)%pindex(2) =
!     &                                   p_list%plist(nloop+2)%pindex(2)
!               t1_list%plist(nloop+2)%pindex(2) = tmp
!            else
!               t1_list%plist(nloop+1)%pindex(2) =
!     &                                   p_list%plist(nloop+2)%pindex(2)
!            end if
!         end if
!
!         if (item%rank2 /= 2) then
!            if (p_list%plist(nloop+1)%ops(2)==2) then
!               tmp = p_list%plist(nloop+1)%pindex(2)
!               t2_list%plist(nloop+1)%pindex(2) =
!     &                                   p_list%plist(nloop+2)%pindex(2)
!               t2_list%plist(nloop+2)%pindex(2) = tmp
!            else
!               t2_list%plist(nloop+1)%pindex(2) =
!     &                                   p_list%plist(nloop+1)%pindex(2)
!            end if
!         end if
!
!      end if
!
!      ! Swap between creation and annihilation operators to make sure
!      ! external ops are before internal. This is important for the
!      ! amplitudes and fock tensors, as these require a specific order of
!      ! indices. Integrals and intermediates are ignored by this
!      ! subroutine
!      ! TODO: will not work with pqrstu...
!      call swap_index(t1_list, item%rank1, item%int(1), item%inter(1))
!      call swap_index(t2_list, item%rank2, item%int(2), item%inter(2))
!      call swap_index(r_list, item%rank3, item%int(3), item%inter(3))
!
!      ! Sort index pairs into order with bubble sort. If only 1 pair,
!      ! then this is skipped. This is important to assure consistent use
!      ! of indices for the declaration and use of an intermediate (the
!      ! slot structure must be the same when it is constructed and when
!      ! it is used)
!      call swap_pairs(t1_list,item%rank1,item%int(1),item%inter(1),e4)
!      e4 = 0
!      call swap_pairs(t2_list,item%rank2,item%int(2),item%inter(2),e4)
!      call swap_pairs(r_list,item%rank3,item%int(3),item%inter(3),e5)
!
!      ! Insert ordered lists into ITF index strings
!      s1 = '        '
!      s2 = '        '
!      s3 = '        '
!
!      do i = 1, item%rank1/2
!         s1(i:i) = t1_list%plist(i)%pindex(1)
!         s1(i+(item%rank1/2):i+(item%rank1/2)) =
!     &                                        t1_list%plist(i)%pindex(2)
!      end do
!      do i = 1, item%rank2/2
!         s2(i:i) = t2_list%plist(i)%pindex(1)
!         s2(i+(item%rank2/2):i+(item%rank2/2)) =
!     &                                        t2_list%plist(i)%pindex(2)
!      end do
!      do i = 1, item%rank3/2
!         s3(i:i) = r_list%plist(i)%pindex(1)
!         s3(i+(item%rank3/2):i+(item%rank3/2)) =
!     &                                         r_list%plist(i)%pindex(2)
!      end do
!
!      ! Assign an ITF index string to each tensor
!      item%idx1 = trim(s1)
!      item%idx2 = trim(s2)
!      item%idx3 = trim(s3)
!
!      ! Release memory for pair lists
!      deallocate(p_list%plist)
!      deallocate(t1_list%plist)
!      deallocate(t2_list%plist)
!      deallocate(r_list%plist)

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
     &   te1(4,2),      ! Test Operator numbers of external index 1
     &   te2(4,2),      ! Test Operator numbers of external index 2
     &   tc(4,2),       ! Test Operator numbers of external index 2
     &   i, j, k, l,    ! Loop index
     &   ii,            ! Letter index
     &   nloop,         ! Number of contraction loops
     &   i1, i2,        ! Search creation or annihilation first
     &   shift,         ! List shift
     &   sp             ! Pair list shift
!      character, dimension(21) ::
!     &   ind=(/ 'a','b','c','d','e','f','g','i','j','k','l','m','n',
!     &          'o','p','q','r','s','t','u','v' /)   ! Letters for index string
      character, dimension(21) ::
     &   ind=(/ 'a','b','c','d','e','f','g','p','q','r','s','t','u',
     &          'v','i','j','k','l','m','n','o' /)   ! Letters for index string
      character(len=INDEX_LEN) ::
     &     s1, s2, s3   ! Tmp ITF index strings
      character(len=1) ::
     &   tmp
      real(8) ::
     &   factor         ! Factor from equivalent lines
      logical ::
     &   found,         ! True if found pairing index
     &   sort           ! Used in bubble sort
      type(pair_list) ::
     &   p_list,        ! Complete list of pairs in binary contraction
     &   t1_list,       ! List of pairs in first tensor (T1)
     &   t2_list,       ! List of pairs in second tensor (T2)
     &   r_list,        ! List of pairs in result tensor
     &   d_list         ! debug list
      integer, dimension(4) ::
     &   t_shift,       ! Index shift for external indices
     &   c_shift        ! Index shift for contraction indices
      integer, dimension(2) ::
     &   cops,          ! Number of creation/annihilation contraction operators
     &   e1ops,         ! Number of C/A external T1 operators
     &   e2ops,         ! Number of C/A external T2 operators
     &   e3ops          ! Total number of external operators
      integer, dimension(4) ::
     &   tops           ! Scratch space

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

      ! Find out number of creation/annihilation operators per operator
      cops = sum(c, dim=1)
      e1ops = sum(e1, dim=1)
      e2ops = sum(e2, dim=1)
      e3ops = sum(e3, dim=1)

      ! Set ranks of tensors
      call itf_rank(e1, c, item%rank1, .false.)
      call itf_rank(e2, c, item%rank2, .false.)
      call itf_rank(e3, c, item%rank3, .true.)

      ! Set number of indcies
      item%nops1 = sum(e1, dim=2) + sum(c, dim=2)
      item%nops2 = sum(e2, dim=2) + sum(c, dim=2)
      item%nops3 = sum(e3, dim=2)

      ! Set if it has three operators of the same type
      ! TODO: probably not needed anymore
      tops = sum(e1, dim=2)
      do i = 1, ngastp
         if (tops(i)==3) then
            item%three(1) = .true.
         end if
      end do

      tops = sum(e1, dim=2)
      do i = 1, ngastp
         if (tops(i)==3) then
            item%three(2) = .true.
         end if
      end do

      if (item%rank1 == 0 .and. item%rank2 == 0) then
         item%idx1 = ''
         item%idx2 = ''
         item%idx3 = ''
         return
      end if

      ! Set number of contraction indicies, used later on
      item%contri = sum(sum(c, dim=1))

!      write(10,*) "e1 ", e1
!      write(10,*) "e2 ", e2
!      write(10,*) "e3 ", e3
!      write(10,*) "c ", c
!      write(10,*) "cops ", cops(1), cops(2)

      ! Allocate pair list according to ranks of tensors
      allocate(p_list%plist(item%rank1/2+item%rank2/2))
      allocate(t1_list%plist(item%rank1/2))
      allocate(t2_list%plist(item%rank2/2))
      allocate(r_list%plist(item%rank3/2))

!      ! New code ======================================================
!      allocate(d_list%plist(item%rank2/2))
!      sp=1
!      do i = 1, 4
!         c_shift(i) = e1(i,1) + e1(i,2) + e2(i,1) + e2(i,2)
!      end do
!      t_shift=0
!      te1 = e1
!      te2 = e2
!      tc = c
!      e1ops = sum(te1, dim=1)
!      e2ops = sum(te2, dim=1)
!
!      !write(item%logfile,*) "e1 ", te1
!      !write(item%logfile,*) "e2 ", te2
!      !write(item%logfile,*) "c ", tc
!      ! Assume exciation operator always comes second
!      if (trim(item%label_t2)=='T2g' .and. item%rank2>2) then
!
!         do while (sum(e2ops) /= 0)
!            ! Search on second tensor
!            ! Note the number of creation/annihilation ops has been
!            ! switched
!            call find_T_pairs(d_list,sp,2, te2, te1,tc,t_shift, c_shift,
!     &                      e2ops, e1ops, item)
!
!            ! Check for more operators
!            e2ops = sum(te2, dim=1)
!         end do
!
!         !call print_plist(d_list, item%rank2/2, "d_list", item%logfile)
!      end if
!
!
!      e1ops = sum(e1, dim=1)
!      e2ops = sum(e2, dim=1)
!      deallocate(d_list%plist)
!      ! End new code ==================================================

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
               ii = 1+(7*(i-1)) + c_shift(i)

               ! Assign index letter to pair list
               p_list%plist(sp)%pindex(i1) = ind(ii)

               ! Assign 'value' to index
               call assign_nval(i, i1, sp, p_list)

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
                     ii = 1+(7*(k-1)) + c_shift(k)

                     p_list%plist(sp)%pindex(i2) = ind(ii)
                     call assign_nval(k, i2, sp, p_list)
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

      !call print_plist(p_list, sp-1, "P_LIST", item%logfile)

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

      !call print_plist(t1_list, item%rank1/2, "T1_LIST", item%logfile)
      !call print_plist(t2_list, item%rank2/2, "T2_LIST", item%logfile)

      ! Create pair list for result tensor, only need external index
      shift = 1
      do i = nloop+1, sp-1
         r_list%plist(shift) = p_list%plist(i)
         shift = shift + 1
      end do

      ! Only need to permute annihilation ops amongst themselves,
      ! this is the case whenever we have (1-Pyx)(1-Pvw). For two rank
      ! 4 tensors, this is straight forward. When there is a rank 2 and
      ! rank 4, can't swap the rank 2, so have to swapped both
      ! annihilations on the rank 4 tensor
      if (item%permute == 2) then
         ! Need to swap annihilation operators between tensors:
         ! T1_{ac}^{ik} T2_{cb}^{kj} -> T1_{ac}^{jk} T2_{cb}^{ki}

         if (item%rank1 /= 2) then
            if (p_list%plist(nloop+2)%ops(2)==1) then
               ! External indices belong on same tensor
               tmp = p_list%plist(nloop+1)%pindex(2)
               t1_list%plist(nloop+1)%pindex(2) =
     &                                   p_list%plist(nloop+2)%pindex(2)
               t1_list%plist(nloop+2)%pindex(2) = tmp
            else
               t1_list%plist(nloop+1)%pindex(2) =
     &                                   p_list%plist(nloop+2)%pindex(2)
            end if
         end if

         if (item%rank2 /= 2) then
            if (p_list%plist(nloop+1)%ops(2)==2) then
               tmp = p_list%plist(nloop+1)%pindex(2)
               t2_list%plist(nloop+1)%pindex(2) =
     &                                   p_list%plist(nloop+2)%pindex(2)
               t2_list%plist(nloop+2)%pindex(2) = tmp
            else
               t2_list%plist(nloop+1)%pindex(2) =
     &                                   p_list%plist(nloop+1)%pindex(2)
            end if
         end if

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
      call swap_pairs(t1_list,item%rank1,item%int(1),item%inter(1),e4)
      e4 = 0
      call swap_pairs(t2_list,item%rank2,item%int(2),item%inter(2),e4)
      call swap_pairs(r_list,item%rank3,item%int(3),item%inter(3),e5)

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


!*----------------------------------------------------------------------*
!      subroutine find_T_pairs(list, sp, tensor, e1, e2, c, t_shift,
!     &                        c_shift, e1ops, e2ops, item)
!*----------------------------------------------------------------------*
!!     Find external pairs in a binary contraction. Assign a contraction
!!     index if they are on different tensors.
!*----------------------------------------------------------------------*
!
!      implicit none
!      include 'opdim.h'
!      include 'def_contraction.h'
!      include 'def_itf_contr.h'
!
!      type(pair_list), intent(inout) ::
!     &   list              ! Complete pair list for binary contraction
!      integer, intent(inout) ::
!     &   sp               ! Pair list shift index
!      integer, intent(inout) ::
!     &   tensor,           ! Label T1 or T2
!     &   e1(4,2),          ! External index occupations for a tensor
!     &   e2(4,2),          ! External index occupations the other tensor
!     &   c(4,2)            ! Contraction index occupations
!      integer, dimension(3), intent(inout) ::
!     &   t_shift,          ! External index letter shift
!     &   c_shift           ! Creation index letter shift
!      integer, dimension(2), intent(in) ::
!     &   e1ops,
!     &   e2ops
!      type(itf_contr), intent(in) ::
!     &   item              ! ITF binary contraction info
!
!      character, dimension(21) ::
!     &   ind=(/ 'a','b','c','d','e','f','g','p','q','r','s','t','u',
!     &          'v','i','j','k','l','m','n','o' /)   ! Letters for index string
!      logical ::
!     &   found
!      integer ::
!     &   opp_tensor,       ! Label opposite tensor to 'tensor'
!     &   i1, i2,           ! Label annihilation/creation index
!     &   n1, n2,           ! Label occupations
!     &   start,            ! Labels P/H/V/X to look for pairs
!     &   ii,               ! Letter index
!     &   i,j,k,l,m,n,z,    ! Loop indices
!     &   nops(4)
!
!
!      ! Decide which tensors we are dealing with (T1 or T2)
!      if (tensor == 1) then
!         opp_tensor = 2
!      else
!         opp_tensor = 1
!      end if
!
!      ! i1 and i2 label creation/annihilation
!      ! operators, pindex always has a creation in position 1, so the
!      ! indices must be placed in their correct position.
!      if (e1ops(1) >= e1ops(2)) then
!         ! Loop through creations first
!         i1 = 1
!         i2 = 2
!         n1 = e1ops(2)
!         n2 = e2ops(2)
!      else
!         ! Loop through annihilation first
!         i1 = 2
!         i2 = 1
!         n1 = e1ops(1)
!         n2 = e2ops(1)
!      end if
!
!      found = .false.
!
!      ! Change i = 2, to start loop at hole index, create pair, then
!      ! change back to i = 1 on next call, continue as normal till every
!      ! index gone
!      ! Have i = start, start set to be 1, 2, 3, 4 based on what the
!      ! largest number of possible loops
!      nops(1) = e1(1,i1) + e1(2,i2) + e1(3,i2) + e1(4,i2)
!      nops(2) = e1(2,i1) + e1(1,i2) + e1(3,i2) + e1(4,i2)
!      nops(3) = e1(3,i1) + e1(1,i2) + e1(2,i2) + e1(4,i2)
!      nops(4) = e1(4,i1) + e1(1,i2) + e1(2,i2) + e1(3,i2)
!
!      ! By default, always start pairing from particle index
!      start = 1
!      do i = 2, size(nops)
!         if (nops(start) < nops(i)) then
!            start = i
!         end if
!      end do
!
!      ! Start main search loop
!      !do i = 1, 4
!      do i = start, 4
!         do j = 1, e1(i,i1)
!            ! Search for first operator
!            ii = 1+(7*(i-1)) + t_shift(i)
!            list%plist(sp)%pindex(i1) = ind(ii)
!            call assign_nval(i, i1, sp, list)
!            list%plist(sp)%ops(i1) = tensor
!
!            t_shift(i) = t_shift(i) + 1
!            e1(i,i1) = e1(i,i1) - 1
!
!            ! Main loop over a and c index
!            do z = 2, 3
!            ! Look for matching operator
!            if (e1(z,i2) > 0) then
!               ! Look on the same tensor (only a or c)
!               do k = 2, 3
!                  do l = 1, e1(k,i2)
!                     ii = 1+(7*(k-1)) + t_shift(k)
!                     list%plist(sp)%pindex(i2) = ind(ii)
!                     call assign_nval(k, i2, sp, list)
!                     list%plist(sp)%linked = .false.
!                     list%plist(sp)%ops(i2) = tensor
!
!                     t_shift(k) = t_shift(k) + 1
!                     e1(k,i2) = e1(k,i2) - 1
!
!                     ! Found a pair, so increment pair list index
!                     sp = sp + 1
!                     found = .true.
!                     exit
!                  end do
!                  if (found) exit
!               end do
!
!            else if (e1(z,i2) > 0) then
!               ! Look on the opposite tensor (only a or c)
!               do k = 2, 3
!                  do l = 1, e2(k,i2)
!                     ii = 1+(7*(k-1))+t_shift(k)
!                     list%plist(sp)%pindex(i2) = ind(ii)
!                     call assign_nval(k, i2, sp, list)
!                     list%plist(sp)%ops(i2) = opp_tensor
!                     list%plist(sp)%linked = .true.
!
!                     t_shift(k) = t_shift(k) + 1
!                     e2(k,i2) = e2(k,i2) - 1
!
!                     ! We need to link the external indices on two
!                     ! different operators with a contraction index
!                     do m=2, 3
!                        do n=1, c(m,i2)
!                            ii = 1+(7*(m-1)) + c_shift(m)
!
!                           list%plist(sp)%link = ind(ii)
!                           call assign_nval(m, 3, sp, list)
!
!                           c_shift(m) = c_shift(m) + 1
!                           c(m,i2) = c(m,i2) - 1
!
!                           sp = sp + 1
!                           found = .true.
!                           exit
!                        end do
!                        if (found) exit
!                     end do
!                     if (found) exit
!                  end do
!                  if (found) exit
!               end do
!            else if (c(z,i2) > 0) then
!               ! Pair with a contraction index
!               do m=2, 3
!                  do n=1, c(m,i2)
!                      ii = 1+(7*(m-1)) + c_shift(m)
!
!                     list%plist(sp)%pindex(i2) = ind(ii)
!                     call assign_nval(m, 2, sp, list)
!                     list%plist(sp)%ops(i2) = opp_tensor
!                     list%plist(sp)%linked = .true.
!                     list%plist(sp)%link = ind(ii)
!                     call assign_nval(m, 3, sp, list)
!
!                     c_shift(m) = c_shift(m) + 1
!                     c(m,i2) = c(m,i2) - 1
!
!                     sp = sp + 1
!                     found = .true.
!                     exit
!                  end do
!                  if (found) exit
!               end do
!            else
!               ! No ops of opposite type to match...
!               call line_error("Particle number not conserving",item)
!            end if
!            end do
!            ! A pair loop has been found so exit
!            if (found) exit
!         end do
!         if (found) exit
!      end do
!
!      return
!      end

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

!      character, dimension(21) ::
!     &   ind=(/ 'a','b','c','d','e','f','g','i','j','k','l','m','n',
!     &          'o','p','q','r','s','t','u','v' /)   ! Letters for index string
      character, dimension(21) ::
     &   ind=(/ 'a','b','c','d','e','f','g','p','q','r','s','t','u',
     &          'v','i','j','k','l','m','n','o' /)   ! Letters for index string
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

      ! Start main search loop
      !do i = 1, 4
      do i = start, 4
         do j = 1, e1(i,i1)
            ! Search for first operator
            ii = 1+(7*(i-1)) + t_shift(i)
            list%plist(sp)%pindex(i1) = ind(ii)
            call assign_nval(i, i1, sp, list)
            list%plist(sp)%ops(i1) = tensor

            t_shift(i) = t_shift(i) + 1
            e1(i,i1) = e1(i,i1) - 1

            ! Look for matching operator
            if (n1 > 0) then
               ! Look on the same tensor
               do k = 1, 4
                  do l = 1, e1(k,i2)
                     ii = 1+(7*(k-1)) + t_shift(k)
                     list%plist(sp)%pindex(i2) = ind(ii)
                     call assign_nval(k, i2, sp, list)
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
                     ii = 1+(7*(k-1))+t_shift(k)
                     list%plist(sp)%pindex(i2) = ind(ii)
                     call assign_nval(k, i2, sp, list)
                     list%plist(sp)%ops(i2) = opp_tensor
                     list%plist(sp)%linked = .true.

                     t_shift(k) = t_shift(k) + 1
                     e2(k,i2) = e2(k,i2) - 1

                     ! We need to link the external indices on two
                     ! different operators with a contraction index
                     do m=1, 4
                        do n=1, c(m,i2)
                            ii = 1+(7*(m-1)) + c_shift(m)

                           list%plist(sp)%link = ind(ii)
                           call assign_nval(m, 3, sp, list)

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
            t_list%plist(s)%nval(1) = p_list%plist(i)%nval(1)

            if (p_list%plist(i)%linked) then
               t_list%plist(s)%pindex(2) = p_list%plist(i)%link
               t_list%plist(s)%nval(2) = p_list%plist(i)%nval(3)
            else
               t_list%plist(s)%pindex(2) = p_list%plist(i)%pindex(2)
               t_list%plist(s)%nval(2) = p_list%plist(i)%nval(2)
            end if

            t_list%plist(s)%linked = .false.
            t_list%plist(s)%ops = place
            s = s + 1
         else if (p_list%plist(i)%ops(2) == place) then
            t_list%plist(s)%pindex(2) = p_list%plist(i)%pindex(2)
            t_list%plist(s)%nval(2) = p_list%plist(i)%nval(2)

            if (p_list%plist(i)%linked) then
               t_list%plist(s)%pindex(1) = p_list%plist(i)%link
               t_list%plist(s)%nval(1) = p_list%plist(i)%nval(3)
            else
               t_list%plist(s)%pindex(1) = p_list%plist(i)%pindex(1)
               t_list%plist(s)%nval(1) = p_list%plist(i)%nval(1)
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
     &   i,                ! Loop index
     &   itmp


      ! If the tensor is a rank-4 integral or an intermediate,
      ! creation/annihilation order is important, so these are skipped
      if (.not. integral .or. integral .and. rank == 2) then
         if (.not. inter) then
            !call print_plist(list, rank/2, "hello1", 10)
            do i = 1, rank/2
!               if (list%plist(i)%pindex(1) >
!     &                                     list%plist(i)%pindex(2)) then
               if (list%plist(i)%nval(1) >=
     &                                     list%plist(i)%nval(2)) then

               if (list%plist(i)%nval(1) ==
     &                                     list%plist(i)%nval(2)) then
               if (list%plist(i)%pindex(1) >
     &                                     list%plist(i)%pindex(2)) then
                  !call print_plist(list, rank/2, "hello2", 10)
                  tmp = list%plist(i)%pindex(1)
                  list%plist(i)%pindex(1) = list%plist(i)%pindex(2)
                  list%plist(i)%pindex(2) = tmp

                  itmp = list%plist(i)%nval(1)
                  list%plist(i)%nval(1) = list%plist(i)%nval(2)
                  list%plist(i)%nval(2) = itmp
               end if
               else
                  tmp = list%plist(i)%pindex(1)
                  list%plist(i)%pindex(1) = list%plist(i)%pindex(2)
                  list%plist(i)%pindex(2) = tmp

                  itmp = list%plist(i)%nval(1)
                  list%plist(i)%nval(1) = list%plist(i)%nval(2)
                  list%plist(i)%nval(2) = itmp
               end if

               end if
            end do
         end if
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine swap_pairs(list, rank, integral, inter, e4)
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
     &   integral,
     &   inter
      integer, intent(in) ::
     &   e4(4,2)

      type(pair_list) ::
     &   tmp_list          ! Temporary holder for pair list
      logical ::
     &   sort              ! True is list has been sorted
      integer ::
     &   i                 ! Loop index

      allocate(tmp_list%plist(1))

      !if (e4(1,2) > e4(1,1)) then
      if (e4(3,2) > e4(3,1)) then
      ! Skip if there is only one pair, or if this is an integral
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


      ! Basically just the amplitude tensors
      if (rank > 2 .and. .not. integral .and. .not. inter) then
         sort = .true.
         do while (sort)
            sort = .false.
            do i = 1, rank/2
               if (i == rank/2) exit

               if (list%plist(i+1)%nval(2) <
     &                                     list%plist(i)%nval(2)) then
                  !call print_plist(list, rank/2, "hello", 10)
                  tmp_list%plist(1) = list%plist(i)
                  list%plist(i) = list%plist(i+1)
                  list%plist(i+1) = tmp_list%plist(1)
                  sort = .true.
               end if
            end do
         end do
      end if

      deallocate(tmp_list%plist)

      return
      end


!*----------------------------------------------------------------------*
!      subroutine permute_tensors(contr_info,perm_case,lulog)
!*----------------------------------------------------------------------*
!!     Find permutation case
!*----------------------------------------------------------------------*
!
!      use itf_utils
!      implicit none
!
!      include 'opdim.h'
!      include 'def_contraction.h'
!
!      type(binary_contr), intent(in) ::
!     &     contr_info   ! Information about binary contraction
!      integer, intent(inout) ::
!     &     perm_case
!      integer, intent(in) ::
!     &     lulog
!
!      integer ::
!     &     e1(4,2),      ! Occupations of external index 1
!     &     e2(4,2),      ! Occupations of external index 2
!     &     c(4,2)
!      integer ::
!     &     i,
!     &     sum_c1,sum_c2,sum_a1,sum_a2
!      logical ::
!     &   inter
!
!      ! Check if not antisym over different vertices
!      ! Check for tensor products
!      ! Not going to antisymm intermediates...
!      ! The intermediates can have eeaa structure
!      e1=0
!      e2=0
!      c=0
!
!      ! TODO: Hack to get correct permuations for terms with inters
!      inter = check_inter(contr_info%label_res)
!
!      ! Get occupation info
!      do i = 1, contr_info%n_cnt
!        call count_index(i,
!     &     contr_info%occ_cnt(1:,1:,i),
!     &     contr_info%rst_cnt(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,c)
!      end do
!      do i = 1, contr_info%nj_op1
!        call count_index(i,
!     &     contr_info%occ_ex1(1:,1:,i),
!     &     contr_info%rst_ex1(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,e1)
!      end do
!      do i = 1, contr_info%nj_op2
!        call count_index(i,
!     &     contr_info%occ_ex2(1:,1:,i),
!     &     contr_info%rst_ex2(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,e2)
!      end do
!
!      ! C: |a|i|p|x|
!      ! A: |a|i|p|x|
!
!      if ((sum(sum(e1,dim=1))+sum(sum(e2,dim=1)))==2) then
!         return
!      else if ((sum(sum(e1,dim=1))+sum(sum(e2,dim=1)))==0) then
!         return
!      else if ((sum(sum(e1,dim=1))+sum(sum(e2,dim=1)))==6) then
!         return
!      else if (sum(sum(c,dim=1))==0) then
!         ! Get rid of tensor products
!         return
!      end if
!
!      perm_case = 0
!
!!      if (e1(1,1)+e2(1,1)==2 .and. e1(2,2)+e2(2,2)==2 .or.
!!     &    e1(3,1)+e2(3,1)==2 .and. e1(2,2)+e2(2,2)==2 .or.
!!     &    e1(1,1)+e2(1,1)==2 .and. e1(3,2)+e2(3,2)==2) then
!      if (e1(1,1)+e2(1,1)==2 .or. e1(1,2)+e2(1,2)==2 .or.
!     &    e1(2,1)+e2(2,1)==2 .or. e1(2,2)+e2(2,2)==2 .or.
!     &    e1(3,1)+e2(3,1)==2 .or. e1(3,2)+e2(3,2)==2) then
!
!         sum_c1=0
!         sum_c2=0
!         sum_a1=0
!         sum_a2=0
!
!
!         do i=1, 4
!            ! Sum creation ops
!            sum_c1=sum_c1+e1(i,1)
!            sum_c2=sum_c2+e2(i,1)
!
!            ! Sum annihilation ops
!            sum_a1=sum_a1+e1(i,2)
!            sum_a2=sum_a2+e2(i,2)
!         end do
!
!         !TODO: Need to detect case (1-P(ab))(1-P(ij)) for
!         !      K_ci^ka T_j^c T_k^b
!         ! Need to distinguish between P(ab) and P(ij)
!         ! logical(2)? If both true get perm_case = 2
!         ! Keep track of this quantaty in item when searching for
!         ! intermediates, if both become true, then we need an extra
!         ! permute
!
!         ! If sum==2, then both indices come from same operator, therefore
!         ! it doesn't need symmetrised
!         if (sum_c1/=2 .and. sum_c2/=2) then
!            if (sum_c1+sum_c2==2) then
!               if (inter) then
!                  return
!               else
!                  perm_case = 1
!               end if
!            end if
!         end if
!
!         if (sum_a1/=2 .and. sum_a2/=2) then
!            if (sum_a1+sum_a2==2) then
!               if (perm_case == 1) then
!                  ! P(ab), then P(abij)
!                  ! When we have (1-P(ab))(1-P(ij)) then we only need to
!                  ! generate 1 + P(ab)
!                  perm_case = 2
!               else
!                  perm_case = 1
!               end if
!            end if
!         end if
!
!      else
!         return
!      end if
!
!      return
!      end


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
     &   item

      integer ::
     &   i,j,              ! Loop index
     &   result_spin(4),   ! Spin case of result
     &   hr1, hr2, hr3     ! Half tensor ranks

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
      else if (item%rank1 == 0 .or. item%rank2 ==0) then
         ! Tensor multiplied by a scalar
         call print_itf_line(item,.false.,.false.)
         return
      else if (item%rank3==4 .and. item%rank1==2
     &         .and. item%rank2==2) then
         ! Don't care about tensor products now
         call print_itf_line(item,.false.,.false.)
         return
      end if


      ! Calculate half ranks for use in indexing letter index strings
      hr1 = item%rank1/2
      hr2 = item%rank2/2
      hr3 = item%rank3/2

      ! Assign spin to indicies on the result tensor
      if (item%inter(3)) then
         item%t_spin(3)%spin = item%i_spin%spin
         !write(item%logfile,*) "i_spin ", item%i_spin%spin
      else
         select case (item%rank3)
            case(0)
            case(2)
               ! aa
               item%t_spin(3)%spin(1,1) = 1
               item%t_spin(3)%spin(2,1) = 1
            case(4)
               ! abab
               item%t_spin(3)%spin(1,1) = 1
               item%t_spin(3)%spin(1,2) = 2
               item%t_spin(3)%spin(2,1) = 1
               item%t_spin(3)%spin(2,2) = 2
            case(6)
               do i=1, 3
                  ! aaaaaa
                  item%t_spin(3)%spin(1,i) = 1
                  item%t_spin(3)%spin(2,i) = 1
               end do
            case default
               call line_error("Could not determine tensor rank",item)
         end select
      end if

      !call print_spin(item%t_spin(3)%spin, item%rank3, "Result", 11)

      ! Assign spin of external indicies to T1 and T2
      do j=1, hr3
         do i=1, hr1
            ! Assign spin of first tensor
            if (item%idx3(j:j)==item%idx1(i:i)) then
               item%t_spin(1)%spin(1,i) = item%t_spin(3)%spin(1,j)
            else if (item%idx3(j:j)==item%idx1(i+hr1:i+hr1)) then
               item%t_spin(1)%spin(2,i) = item%t_spin(3)%spin(1,j)
            end if

            if (item%idx3(j+hr3:j+hr3)==item%idx1(i:i)) then
               item%t_spin(1)%spin(1,i) = item%t_spin(3)%spin(2,j)
            else if (item%idx3(j+hr3:j+hr3)==item%idx1(i+hr1:i+hr1))then
               item%t_spin(1)%spin(2,i) = item%t_spin(3)%spin(2,j)
            end if
         end do

         do i=1, hr2
            ! Assign spin of second tensor
            if (item%idx3(j:j)==item%idx2(i:i)) then
               item%t_spin(2)%spin(1,i) = item%t_spin(3)%spin(1,j)
            else if (item%idx3(j:j)==item%idx2(i+hr2:i+hr2)) then
               item%t_spin(2)%spin(2,i) = item%t_spin(3)%spin(1,j)
            end if

            if (item%idx3(j+hr3:j+hr3)==item%idx2(i:i)) then
               item%t_spin(2)%spin(1,i) = item%t_spin(3)%spin(2,j)
            else if (item%idx3(j+hr3:j+hr3)==item%idx2(i+hr2:i+hr2))then
               item%t_spin(2)%spin(2,i) = item%t_spin(3)%spin(2,j)
            end if
         end do
      end do

      !call print_spin(item%t_spin(1)%spin, item%rank1, "T1", 11)
      !call print_spin(item%t_spin(2)%spin, item%rank2, "T2", 11)

      ! Sum over the remaining contraction indicies and print out the
      ! line
      call spin_index(item)

      return
      end


*----------------------------------------------------------------------*
      subroutine spin_index(item)
*----------------------------------------------------------------------*
!     Find contraction index used in spin summation
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &     item

      type(twodarray), pointer ::
     &   poss(:,:) => null()
      integer ::
     &   i, j, k, l, m, n, shift,
     &   i1, i2, i3, i4,
     &   z1, z2, r1, r2
      character(len=INDEX_LEN) ::
     &   str1, str2
      logical ::
     &   eloop

      allocate(poss(2,item%contri))

      !do i = 1, 2
      !   do j = 1, item%contri
      !      poss(i,j)%elements = 0
      !   end do
      !end do

      ! Largest tensor goes first
      if (item%rank2 > item%rank1) then
         z1 = 2
         z2 = 1
         r1 = item%rank2/2
         r2 = item%rank1/2
         str1 = item%idx2
         str2 = item%idx1
      else
         z1 = 1
         z2 = 2
         r1 = item%rank1/2
         r2 = item%rank2/2
         str1 = item%idx1
         str2 = item%idx2
      end if

      !call print_spin(item%t_spin(1)%spin, item%rank1, "T1", 10)
      !call print_spin(item%t_spin(2)%spin, item%rank2, "T2", 10)

      ! Get position information of contraction indicies
      shift = 1
      do i=1, size(item%t_spin(z1)%spin,2)
         if (item%t_spin(z1)%spin(1,i) == 0) then
            do j=1, size(item%t_spin(z2)%spin,2)
               if (str1(i:i) == str2(j+r2:j+r2)) then
                  poss(1,shift)%elements(1) = 1
                  poss(1,shift)%elements(2) = i
                  poss(2,shift)%elements(1) = 2
                  poss(2,shift)%elements(2) = j
                  shift = shift + 1
               end if
            end do
            do j=1, size(item%t_spin(z2)%spin,2)
               if (str1(i:i) == str2(j:j)) then
                  poss(1,shift)%elements(1) = 1
                  poss(1,shift)%elements(2) = i
                  poss(2,shift)%elements(1) = 1
                  poss(2,shift)%elements(2) = j
                  shift = shift + 1
               end if
            end do
         end if
      end do

      do i=1, size(item%t_spin(z1)%spin,2)
         if (item%t_spin(z1)%spin(2,i) == 0) then
            do j=1, size(item%t_spin(z2)%spin,2)
               if (str1(i+r1:i+r1) == str2(j:j)) then
                  poss(1,shift)%elements(1) = 2
                  poss(1,shift)%elements(2) = i
                  poss(2,shift)%elements(1) = 1
                  poss(2,shift)%elements(2) = j
                  shift = shift + 1
               end if
            end do
            do j=1, size(item%t_spin(z2)%spin,2)
               if (str1(i+r1:i+r1) == str2(j+r2:j+r2)) then
                  poss(1,shift)%elements(1) = 2
                  poss(1,shift)%elements(2) = i
                  poss(2,shift)%elements(1) = 2
                  poss(2,shift)%elements(2) = j
                  shift = shift + 1
               end if
            end do
         end if
      end do

      if (shift == 1) then
         call line_error("Didn't find contraction index", item)
      end if

      ! Main spin summation loop
      shift = shift - 1
      eloop = .false.

      !write(10,*) "========================"
      !do i = 1, 2
      !   do j = 1, item%contri
      !      write(10,*) "DEBUG1 poss: ", poss(i,j)%elements
      !   end do
      !end do
      !write(10,*) "========================"

      do i = 1, 2
         i1 = poss(1,1)%elements(1)
         i2 = poss(1,1)%elements(2)
         i3 = poss(2,1)%elements(1)
         i4 = poss(2,1)%elements(2)
         item%t_spin(z1)%spin(i1, i2) = i
         item%t_spin(z2)%spin(i3, i4) = i
         if (shift <= 1) then
            call print_spin_case(item,eloop)
         end if
         if (shift > 1) then
            do j = 1, 2
               i1 = poss(1,2)%elements(1)
               i2 = poss(1,2)%elements(2)
               i3 = poss(2,2)%elements(1)
               i4 = poss(2,2)%elements(2)
               item%t_spin(z1)%spin(i1, i2) = j
               item%t_spin(z2)%spin(i3, i4) = j
               if (shift <= 2) then
                  call print_spin_case(item,eloop)
               end if
               if (shift > 2) then
                  do k = 1, 2
                     i1 = poss(1,3)%elements(1)
                     i2 = poss(1,3)%elements(2)
                     i3 = poss(2,3)%elements(1)
                     i4 = poss(2,3)%elements(2)
                     item%t_spin(z1)%spin(i1, i2) = k
                     item%t_spin(z2)%spin(i3, i4) = k
                     if (shift <= 3) then
                        call print_spin_case(item,eloop)
                     end if
                     if (shift > 3) then
                        do l = 1, 2
                           i1 = poss(1,4)%elements(1)
                           i2 = poss(1,4)%elements(2)
                           i3 = poss(2,4)%elements(1)
                           i4 = poss(2,4)%elements(2)
                           item%t_spin(z1)%spin(i1, i2) = l
                           item%t_spin(z2)%spin(i3, i4) = l
                           if (shift <= 4) then
                              call print_spin_case(item,eloop)
                           end if
                           if (shift > 4) then
                              do m = 1, 2
                                 i1 = poss(1,5)%elements(1)
                                 i2 = poss(1,5)%elements(2)
                                 i3 = poss(2,5)%elements(1)
                                 i4 = poss(2,5)%elements(2)
                                 item%t_spin(z1)%spin(i1, i2) = m
                                 item%t_spin(z2)%spin(i3, i4) = m
                                 if (shift <= 5) then
                                    call print_spin_case(item,eloop)
                                 end if
                                 if (shift > 5) then
                                    do n = 1, 2
                                       i1 = poss(1,6)%elements(1)
                                       i2 = poss(1,6)%elements(2)
                                       i3 = poss(2,6)%elements(1)
                                       i4 = poss(2,6)%elements(2)
                                       item%t_spin(z1)%spin(i1, i2) = n
                                       item%t_spin(z2)%spin(i3, i4) = n
                                       if (shift <= 6) then
                                        call print_spin_case(item,eloop)
                                       end if
                                    end do
                                 end if
                              end do
                           end if
                        end do
                     end if
                  end do
               end if
            end do
         end if
      end do

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

      deallocate(poss)

      return
      end


*----------------------------------------------------------------------*
      subroutine print_spin_case(item,eloop)
*----------------------------------------------------------------------*
!     Print spin case
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &   item
      logical, intent(inout) ::
     &   eloop     ! Check if at least one spin case is printed out

      logical ::
     &   s1,       ! True if tensor 1 is mixed spin
     &   s2,       ! True if tensor 2 is mixed spin
     &   contains1,
     &   contains2
      integer ::
     &   i,
     &   shift,   ! Index number of different spin cases for an intermediate
     &   ishift  ! Index number of intermediates on a line
      character(len=INDEX_LEN) ::
     &   spin_name
      integer ::
     &   sum1a, sum1b, sum2a, sum2b

      ! Sum the spins of the covarient and contravarient indicies
      sum1a = 0
      sum1b = 0
      sum2a = 0
      sum2b = 0
      do i = 1, size(item%t_spin(1)%spin,2)
         sum1a = sum1a + item%t_spin(1)%spin(1,i)
         sum1b = sum1b + item%t_spin(1)%spin(2,i)
      end do
      do i = 1, size(item%t_spin(2)%spin,2)
         sum2a = sum2a + item%t_spin(2)%spin(1,i)
         sum2b = sum2b + item%t_spin(2)%spin(2,i)
      end do

      ! Pick out specific spin cases here
      if (sum1a==sum1b .and. sum2a==sum2b) then
         if (modulo(sum1a+sum1b,2)==0 .and.
     &           modulo(sum2a+sum2b,2)==0) then

            ! Decide if tensor is mixed spin
            contains1 = .false.
            contains2 = .false.
            s1 = .false.
            if (item%rank1>2) then
               do i = 1, size(item%t_spin(1)%spin,2)
                  if (item%t_spin(1)%spin(1,i)==1) then
                     contains1 = .true.
                  end if
                  if (item%t_spin(1)%spin(1,i)==2) then
                     contains2 = .true.
                  end if
               end do
               if (contains1 .neqv. contains2) then
                  s1 = .true.
               end if
            end if

            contains1 = .false.
            contains2 = .false.
            s2 = .false.
            if (item%rank2>2) then
               do i = 1, size(item%t_spin(2)%spin,2)
                  if (item%t_spin(2)%spin(1,i)==1) then
                     contains1 = .true.
                  end if
                  if (item%t_spin(2)%spin(1,i)==2) then
                     contains2 = .true.
                  end if
               end do
               if (contains1 .neqv. contains2) then
                  s2 = .true.
               end if
            end if

            ! Append spin name to intermediate label (i.e. STIN001abab)
            if (item%inter(1)) then
               call inter_spin_name(item%t_spin(1)%spin,
     &                                 item%rank1/2,item%inter1)
            end if

            if (item%inter(2)) then
               call inter_spin_name(item%t_spin(2)%spin,
     &                                 item%rank2/2,item%inter2)
            end if


            ! Check if we are searching for intermediates
            if (associated(item%inter_spins)) then

               ! Number of intermediates
               ishift = 0

               if (item%inter(1)) then
                  call inter_spin_case(item%t_spin(1)%spin,item%rank1/2,
     &                                 item%label_t1,ishift,item)
               end if

               if (item%inter(2)) then
                  call inter_spin_case(item%t_spin(2)%spin,item%rank2/2,
     &                                  item%label_t2,ishift,item)
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

               call print_itf_line(item,s1,s2)
            end if

            eloop=.true.
         end if
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine inter_spin_name(spin,hrank,label)
*----------------------------------------------------------------------*
!     Add spin name to intermediate (ie. STIN001 -> STIN001abab)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   hrank,            ! Rank/2
     &   spin(2,hrank)    ! Spin info
      character(len=INDEX_LEN), intent(inout) ::
     &   label             ! Spin name of intermediate

      character(len=INDEX_LEN) ::
     &   spin_name
      integer ::
     &   i

      spin_name = ''

      do i = 1, hrank
         if (spin(1,i)==1) then
            spin_name(i:i) = 'a'
         else if (spin(1,i)==2) then
            spin_name(i:i) = 'b'
         end if

         if (spin(2,i)==1) then
            spin_name(i+hrank:i+hrank) = 'a'
         else if (spin(2,i)==2) then
            spin_name(i+hrank:i+hrank) = 'b'
         end if
      end do

      label = spin_name

      return
      end


*----------------------------------------------------------------------*
      subroutine inter_spin_case(spin,hrank,label,ishift,item)
*----------------------------------------------------------------------*
!     Determine intermediate spin cases and save info to item
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   hrank,       ! Half rank
     &   spin(2,hrank)
      integer, intent(inout) ::
     &   ishift      ! Index number of intermediates
      character(len=MAXLEN_BC_LABEL), intent(inout) ::
     &   label
      type(itf_contr), intent(inout) ::
     &   item

      integer ::
     &   i,
     &   shift

      ishift = ishift + 1

      ! Number of spin cases
      shift = item%inter_spins(ishift)%ncase + 1

      ! For now, if intermediate is part of a permutation
      ! line, then we add a P to its name
      if (item%permute == 2) then
      item%inter_spins(ishift)%name=trim(label)//'P'
      else
         item%inter_spins(ishift)%name=label
      end if

      do i=1, hrank
         item%inter_spins(ishift)%cases(i,shift)=spin(1,i)
         item%inter_spins(ishift)%cases(i+INDEX_LEN/2,shift)=spin(2,i)
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine itf_contr_init(contr_info,item,perm,comm,
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
     &     item     ! Object which holds information necessary to print out an ITF algo line
      integer, intent(in) ::
     &     perm,        ! Permutation information
     &     comm,        ! formula_item command
     &     lulog        ! Output file
      ! This should be in a 'constructor'

      ! Assign output file
      item%logfile=lulog

      ! Assign command type
      item%command=comm

      ! Check if binary contraction or not
      if (comm==command_add_intm .or. comm==command_cp_intm) then
         item%binary = .false.
      end if


      ! Assign permutation number
      if (comm/=command_cp_intm .or. comm/=command_add_intm) then
         item%permute=perm
      end if

      ! Assign labels
      item%label_t1=contr_info%label_op1
      if (comm/=command_cp_intm .or. comm/=command_add_intm) then
         ! Operator 2 does not exist in [ADD] or [COPY]
         item%label_t2=contr_info%label_op2
      end if
      item%label_res=contr_info%label_res


      ! Check if an intermediate
      item%inter(1) = check_inter(item%label_t1)
      if (comm/=command_cp_intm .or. comm/=command_add_intm) then
         item%inter(2) = check_inter(item%label_t2)
      end if
      item%inter(3) = check_inter(item%label_res)

      ! Check if an integral
      item%int(1) = check_int(item%label_t1)
      if (comm/=command_cp_intm .or. comm/=command_add_intm) then
         item%int(2) = check_int(item%label_t2)
      end if
      item%int(3) = .false.

      ! Assign factor --- use special ITF factor
      ! the ITF factor is closer to the value expected from standard
      ! diagram rules, but still some care has to be taken when translating
      ! to ITF tensors; e.g. the (HP;HP) integrals are stored by GeCCo as
      ! <aj||bi> while the standard sign for ring terms assumes <aj||ib>
      item%fact=contr_info%fact_itf
      !write(11,*) "diag fact: ", contr_info%fact
      !write(11,*) "itf fact: ", contr_info%fact_itf

      ! Account for negative sign as explained from above...
      call integral_fact(contr_info,item%fact)

      !write(11,*) "integral fact: ", item%fact

      ! Inialise number of contraction indicies
      item%contri = 0

      ! Assign index string. Tensor ranks are also set here
      if (comm==command_cp_intm .or. comm==command_add_intm) then
         ! For [ADD] and [COPY]
         call assign_add_index(contr_info,item)
      else
         ! For other contractions
         call assign_index(contr_info,item)
         call assign_new_index(contr_info,item)
      end if

      item%inter1 = ''
      item%inter2 = ''

      allocate(item%t_spin(1)%spin(2, item%rank1/2))
      allocate(item%t_spin(2)%spin(2, item%rank2/2))
      allocate(item%t_spin(3)%spin(2, item%rank3/2))

      item%t_spin(1)%spin = 0
      item%t_spin(2)%spin = 0
      item%t_spin(3)%spin = 0

      if (item%inter(3)) then
         allocate(item%i_spin%spin(2, item%rank3/2))
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine itf_deinit(item)
*----------------------------------------------------------------------*
!     Deinitialise ITF contraction object
*----------------------------------------------------------------------*

      use itf_utils
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &     item     ! Object which holds information necessary to print out an ITF algo line

      deallocate(item%t_spin(1)%spin)
      deallocate(item%t_spin(2)%spin)
      deallocate(item%t_spin(3)%spin)

      if (item%inter(3)) then
         deallocate(item%i_spin%spin)
      end if

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
     &   loop,
     &     c(4,2)


      c = 0

      if (contr%label_op1 == 'H') then
         !TODO: Why looping over n_cnt????? probably want to loop over
         !      number of blocks for each operator
         do i = 1, contr%n_cnt
!            write(10,*) "n_cnt ", contr%n_cnt
!            write(10,*) "occ_op1 ", contr%occ_op1(1:,1:,i)
!            write(10,*) "rst_op1 ", contr%rst_op1(1:,1:,1:,1:,1:,i)
           call count_index(i,
     &        contr%occ_op1(1:,1:,i),
     &        contr%rst_op1(1:,1:,1:,1:,1:,i),
     &        contr%ngas,contr%nspin,c)
         end do
      end if
      if (contr%label_op2 == 'H') then
         !TODO: Why looping over n_cnt????? probably want to loop over
         !      number of vertices for each operator
         do i = 1, contr%n_cnt
           call count_index(i,
     &        contr%occ_op2(1:,1:,i),
     &        contr%rst_op2(1:,1:,1:,1:,1:,i),
     &        contr%ngas,contr%nspin,c)
         end do
      end if

      !if (c(1,1)==1 .and. c(2,1)==1 .and. c(1,2)==1 .and. c(2,2)==1)
      if (c(1,1)==1 .and. c(3,1)==1 .and. c(1,2)==1 .and. c(3,2)==1)
     &   then
            fact = fact * -1.0d+0
            !write(11,*) "Changing the factor ", fact
      end if

      ! Catch three internal integrals
      !if (c(2,1) + c(2,2)==3) then
      if (c(3,1) + c(3,2)==3) then
         fact = fact * -1.0d+0
         !write(11,*) "3 internal"
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

      character(len=MAXLEN_BC_LABEL), intent(inout) ::
     &     result
      type(itf_contr), intent(in) ::
     &     item

      character(len=70) ::
     &     line
      character(len=INDEX_LEN) ::
     &     tindex
      character(len=MAXLEN_BC_LABEL) ::
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
