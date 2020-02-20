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
      else if (trim(string).eq.'O3') then
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
             rename_tensor='Ym1'
          else if (rank==4) then
             rename_tensor='Dm2'
          else if (rank==6) then
             rename_tensor='Dm3'
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
      pure function check_den(label)
*----------------------------------------------------------------------*
!     Check if tensor is a density matrix
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &     label

      logical ::
     &     check_den

      ! Assume these are the names of intermediates
      if (index(label, "GAM")>0) then
         check_den=.true.
      else
         check_den=.false.
      end if

      end function


*----------------------------------------------------------------------*
      pure function check_energy(label)
*----------------------------------------------------------------------*
!     Check if tensor is a density matrix
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &     label

      logical ::
     &     check_energy

      ! Assume these are the names of intermediates
      if (index(label, "ECCD")>0) then
         check_energy=.true.
      else
         check_energy=.false.
      end if

      end function


*----------------------------------------------------------------------*
      function get_itype(idx, canonical) result(itype)
*----------------------------------------------------------------------*
!     Returns index type (itype) of a given index
*----------------------------------------------------------------------*
      implicit none

      character(len=1), intent(in) ::
     &   idx
      logical, optional, intent(in) ::
     &   canonical
      integer ::
     &   itype,
     &   vals(4)
      logical ::
     &   can

      if (present(canonical)) then
         can = canonical
      else
         can = .false.
      end if

      if (can) then
         vals(1) = 2
         vals(2) = 1
         vals(3) = 3
         vals(4) = 4
      else
         vals(1) = 1
         vals(2) = 3
         vals(3) = 2
         vals(4) = 4
      end if

      if (scan("abcdefgh",idx)>0) then
         itype = vals(1)
      else if (scan("ijklmno",idx)>0) then
         itype = vals(2)
      else if (scan("pqrstuvw",idx)>0) then
         itype = vals(3)
      else if (scan("xyz",idx)>0) then
         itype = vals(4)
      end if

      end function get_itype


*----------------------------------------------------------------------*
      recursive function c_index(idx, n, reverse) result(cidx)
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
      logical, optional, intent(in) ::
     &   reverse

      character(len=INDEX_LEN) ::
     &   cidx
      character(len=INDEX_LEN) ::
     &   tmp

      logical ::
     &   rev

      if (present(reverse)) then
         rev = reverse
      else
         rev = .false.
      end if


      !TODO: only works for rank 6
      if (n == 0) then
         cidx = idx
      else
         if (rev) then
            tmp = idx
            tmp(1:1) = idx(2:2)
            tmp(2:2) = idx(3:3)
            tmp(3:3) = idx(1:1)
            cidx = c_index(tmp, n-1, rev)
         else
            tmp = idx
            tmp(1:1) = idx(3:3)
            tmp(2:2) = idx(1:1)
            tmp(3:3) = idx(2:2)
            cidx = c_index(tmp, n-1)
         end if
      end if

      end function c_index


*----------------------------------------------------------------------*
      pure function f_index(index, hrank, upper, middle, first, second)
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
     &   upper,         ! Flip contravarient index
     &   middle,        ! Flip index slots 2 and 3
     &   first,         ! Flip first pair
     &   second         ! Flip second pair
      character(hrank*2) ::
     &   f_index        ! Transpose of ITF index string

      logical ::
     &   contra,
     &   mid,
     &   f_pair,
     &   s_pair

      f_index=index
      if (hrank==1 .or. hrank==0) return

      if (present(upper)) then
         contra = upper
      else
         contra = .false.
      end if

      if (present(middle)) then
         mid = middle
      else
         mid = .false.
      end if

      if (present(first)) then
         f_pair = first
      else
         f_pair = .false.
      end if

      if (present(second)) then
         s_pair = second
      else
         s_pair = .false.
      end if

      if (contra) then
         f_index(hrank-1+hrank:hrank-1+hrank)=
     &                                index(hrank+hrank:hrank+hrank)
         f_index(hrank+hrank:hrank+hrank)=
     &                                index(hrank-1+hrank:hrank-1+hrank)
      else if (mid) then
         f_index(2:2)=index(3:3)
         f_index(3:3)=index(2:2)
      else if (f_pair) then
         f_index(1:1)=index(3:3)
         f_index(3:3)=index(1:1)
      else if (s_pair) then
         f_index(2:2)=index(4:4)
         f_index(4:4)=index(2:2)
      else
         f_index(hrank-1:hrank-1)=index(hrank:hrank)
         f_index(hrank:hrank)=index(hrank-1:hrank-1)
      end if

      end function


      end module itf_utils


*----------------------------------------------------------------------*
      subroutine init_spin_cases(spin_inters)
*----------------------------------------------------------------------*
!     Initalise spin_cases array
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(spin_cases), dimension(MAXINT), intent(inout) ::
     &     spin_inters

      integer ::
     &   i, j, k

      do i = 1, MAXINT
         spin_inters(i)%name = ''
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine print_inter_spin_cases(spin_inters, ninter, label,
     &                                  logfile)
*----------------------------------------------------------------------*
!     Print out intermediate spin cases
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(spin_cases), dimension(MAXINT), intent(in) ::
     &     spin_inters
      integer, intent(in) ::
     &   ninter,
     &   logfile
      character(len=*), intent(in) ::
     &   label

      integer ::
     &   i, j, k

      if (ninter==0) return

      write(logfile,*) "=============================="
      write(logfile,*) "TITLE: ", trim(label)

      do i = 1, ninter
         write(logfile,*) "=============================="
         write(logfile,*) "NCASE: ", spin_inters(i)%ncase
         write(logfile,*) "NAME: ", trim(spin_inters(i)%name)
         write(logfile,*) "SYMM_RES: ", spin_inters(i)%symm_res
         do j = 1, spin_inters(i)%ncase
            write(logfile,*) "Spin:"
            write(logfile,'(4i2)')
     &                        (spin_inters(i)%cases(k,j),k=1,INDEX_LEN)
            write(logfile,*) "------------------------------"
         end do

         write(logfile,'(8i2)')
     &                         (spin_inters(i)%itype(j),j=1,INDEX_LEN)
         write(logfile,*) "------------------------------"
         write(logfile,*) "=============================="
      end do
      write(logfile,*)

      return
      end


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
      subroutine print_itf_contr(item, label)
*----------------------------------------------------------------------*
!     Print line of ITF code
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(in) ::
     &   item
      character(len=*), intent(in) ::
     &   label

      integer ::
     &   i,j, f

      write(item%logfile,*) "ift_contr: ", label
      write(item%logfile,*) "=============================="
      write(item%logfile,*) "R:  ", item%label_res, "  ", item%idx3
      write(item%logfile,*) "T1: ", item%label_t1, "  ", item%idx1
      write(item%logfile,*) "T2: ", item%label_t2, "  ", item%idx2
      write(item%logfile,*) "---------------------------"
      write(item%logfile,*) "itype:"
      write(item%logfile,'(8i3)') (item%itype(j), j=1, INDEX_LEN)
      write(item%logfile,*) "factor ", item%fact
      if (item%inter(3)) then
         write(item%logfile,*) "i_spin:"
         do i = 1, 2
            write(item%logfile,*)
     &                        (item%i_spin%spin(i,j),j=1,item%rank3/2)
            write(item%logfile,*) "------------------------------"
         end do
      end if
      write(item%logfile,*) "=============================="

      return
      end

*----------------------------------------------------------------------*
      subroutine count_index2(iocc, nops)
*----------------------------------------------------------------------*
!     Return array with number of operators of each type
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     iocc(ngastp,2)
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
      subroutine command_to_itf2(contr_info, itin, itflog, command,
     &                           inter_itype)
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
     &   contr_info      ! Information about binary contraction
      logical, intent(in) ::
     &   itin              ! Print ITIN lines or not
      integer, intent(in) ::
     &   itflog,         ! Output file
     &   command,        ! Type of formula item command, ie. contraction, copy etc.
     &   inter_itype     ! Store itypes of intermediates between lines

      type(itf_contr2) ::
     &   item,      ! ITF contraction object; holds all info about the ITF algo line
     &   pitem      ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &   perm_case,   ! Info of permutation factors
     &   i, j, l, k,                ! Loop index
     &   ntest = 00
      logical ::
     &   inter,           ! True if result is an intermediate
     &   upper,
     &   symmetric,
     &   intpp,
     &   pline
      character(len=MAXLEN_BC_LABEL) ::
     &   old_name,
     &   un_perm_name,
     &   old_inter
      character(len=INDEX_LEN) ::
     &   old_idx,
     &   un_perm_idx

!      ! Being a special block which the python processor will pull out
!      ! into its own code block
!      intpp = .false.
!      if (contr_info%label_res=='INTpp') then
!         intpp = .true.
!         write(itflog,'(a)') "BEGIN_INTPP"
!      end if
!
!
      intpp = .false.

      ! Mark begining of spin summed block
      write(itflog,'(a5)') 'BEGIN'


      ! 1. Initalise itf_contr
      call itf_contr_init2(contr_info,item,1,itin,command,itflog,
     &                     inter_itype)


      ! 2. Assign index / Determine sign
      call assign_new_index2(item)


      ! 3. Determine if we need a permutation line and create new item
      !    for it
      !    Permutation line required for:
      !        a) Symmetric residuals with two pairs of external indicies
      !        a) Non-symmetric residuals with one pair of external indicies
      pline = .false.
      perm_case = 0
      do i = 1, ngastp
         if(contr_info%perm(i)) perm_case = perm_case + 1
      end do

      ! 3.5. Decide whether to symmetrise after every term
      ! TODO: Symmetrise at the end for some results
      ! TODO: keep old name in item
      old_name = item%label_res
      item%intpp = .false.
      call prepare_symmetrise2(contr_info, command, perm_case,
     &                         old_name, item)


      if (item%symmetric .and. perm_case==2 .or.
     &    .not. item%symmetric .and. perm_case==1) then
         ! Don't need perm line for intermediates
         if (.not. item%inter(3)) then
         pline = .true.

         !TODO: just copy item??
         call itf_contr_init2(contr_info,pitem,1,itin,command,itflog,
     &                        inter_itype)
         call assign_new_index2(pitem)

         ! TODO: For now, so perm lines have same name as non-perm lines
         pitem%label_res = item%label_res

         ! Permute indicies and update factor
         ! -1 factor for permuation line
         call create_permutation2(pitem, contr_info%perm)
         !call print_itf_contr2(pitem)
         end if
      end if


      ! 4. Spin sum
      call assign_spin2(item)
      if (pline) call assign_spin2(pitem)

      ! 5. Loop over spin cases and print out each line
      call print_spin_cases(item)
      if (pline) then
         call print_spin_cases(pitem)
      end if

      ! 7. Print symmetrisation term
      if (.not. item%inter(3)) then
         call print_symmetrise2(old_name, item)
      end if

      ! 8. If an intermediate, set inter_itype for use in next line
      ! where the intermediate is created
      call set_itype(item, inter_itype)


      if (ntest>0) then
      !if (.true.) then
         call print_itf_contr2(item)
         if (pline) call print_itf_contr2(pitem)
      end if


      ! Mark end of spin block
      write(itflog,'(a)') "END"

      ! Deallocate memroy used when construcitng item
      call itf_deinit2(item)

      !if (intpp) then
      !   write(itflog,'(a)') "END_INTPP"
      !end if

      return
      end



*----------------------------------------------------------------------*
      subroutine command_to_itf(contr_info, itin, itflog, command)
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
     &   contr_info      ! Information about binary contraction
      logical, intent(in) ::
     &   itin              ! Print ITIN lines or not
      integer, intent(in) ::
     &   itflog,         ! Output file
     &   command         ! Type of formula item command, ie. contraction, copy etc.

      type(itf_contr) ::
     &   item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &   perm_case,   ! Info of permutation factors
     &   i, j, l, k                ! Loop index
      logical ::
     &   inter,           ! True if result is an intermediate
     &   found,
     &   upper,
     &   symmetric,
     &   intpp
      character(len=MAXLEN_BC_LABEL) ::
     &   old_name,
     &   un_perm_name,
     &   old_inter
      character(len=INDEX_LEN) ::
     &   old_idx,
     &   un_perm_idx

      ! Being a special block which the python processor will pull out
      ! into its own code block
      intpp = .false.
      if (contr_info%label_res=='INTpp') then
         intpp = .true.
         write(itflog,'(a)') "BEGIN_INTPP"
      end if


      ! Initialise permutation factors:
      ! 0 == no permutation
      ! 1 == (1-Pxy)
      ! 2 == (1-Pxy)(1-Pvw) = (1+Pxy) in spatial orbitals
      perm_case = 0
      do i = 1, ngastp
         if(contr_info%perm(i)) perm_case = perm_case + 1
      end do

      ! If symmetrising after every term, rename residual ITIN
      call prepare_symmetrise(contr_info, itin, intpp, symmetric,
     &                            command, old_name)


      ! Mark begining of spin summed block
      write(itflog,'(a5)') 'BEGIN'

      ! Pick out specific commands, form the itf_contr object, spin sum
      ! and print out contraction line
      if (command==command_add_intm .or. command==command_cp_intm) then
         ! For [ADD] and [COPY] cases
         call itf_contr_init(contr_info,item,0,itin,command,itflog)
         call print_itf_line(item,.false.,.false.)
      else
         ! For other binary contractions
         if (perm_case == 0) then
            ! No permutations
            call itf_contr_init(contr_info,item,0,itin,command,itflog)
            call assign_spin(item)
         else
            if (symmetric) then
!               call itf_contr_init(contr_info,item,1,itin,command,
!     &                                itflog)
!               call assign_spin(item)
!
!               ! check_permutation()
!               if (perm_case==2) then
!
!                  ! create_permutation_line()
!
!               end if

               do i=1, perm_case
                  ! Loop over permutation cases and send separately to
                  ! assign_spin. For most cases this is just one, however
                  ! for (1-Pij)(1-Pab), we need to generate one of these
                  ! permutations before symmetrising
                  call itf_contr_init(contr_info,item,i,itin,command,
     &                                itflog)


                  if (i == 2) then

                     if (item%rank1==6) then
                        !item%idx1=c_index(item%idx1,1,.true.)
                        item%idx1=c_index(item%idx1,1)
                     else
                        item%idx1=f_index(item%idx1,item%rank1/2)
                        !if (item%rank1/=0 .and. item%rank2==0) then
                        !   item%idx1=f_index(item%idx1,item%rank1/2,.true.)
                        !end if
                     end if

                     item%idx2=f_index(item%idx2,item%rank2/2)

                     ! Whenever we tranpose a tensor, we intoroduce a sign
                     ! chage
                     ! No sign due to the tranpose of idx1 which defines an
                     ! intermeidate. Extra signs to to tranpose of tensors
                     ! which define intermediates are included in the
                     ! intermediate line
                     if (item%permute==2) then
                        if (item%rank2>2) item%fact = item%fact*-1.0d+0
                        !write(item%logfile,*)"index flip fact: ", item%fact
                     end if

                  end if

                  call assign_spin(item)

                  !if (i == 1 .and. perm_case>1) then
                  !   call create_permutation(item, contr_info%perm)
                  !end if
               end do
            else
               un_perm_name=''
               un_perm_idx=''
               do i=1, perm_case+1
                  ! Loop over permutation cases and send separately to
                  ! assign_spin. For most cases this is just one, however
                  ! for (1-Pij)(1-Pab), we need to generate one of these
                  ! permutations before symmetrising
                  call itf_contr_init(contr_info,item,i,itin,
     &                                command,itflog)

                  ! Don't print permutation iter if we already have it
                  if (i==2 .and. item%inter(1)) then
                     if (un_perm_idx==item%idx1 .and.
     &                   un_perm_name==item%label_t1) then
                        item%print_line = .false.
                     end if
                  else if (i==2 .and. item%inter(2)) then
                     if (un_perm_idx==item%idx2 .and.
     &                   un_perm_name==item%label_t2) then
                        item%print_line = .false.
                     end if
                  end if

                  call assign_spin(item)

                  if (i==1 .and. item%inter(1)) then
                     un_perm_idx = item%idx1
                     un_perm_name = item%label_t1
                  else if (i==1 .and. item%inter(2)) then
                     un_perm_idx = item%idx2
                     un_perm_name = item%label_t2
                  end if

               end do

               if (item%inter(1)) then
                  if (un_perm_idx==item%idx1 .and.
     &                un_perm_name==item%label_t1) then
                     write(itflog,'(a)') "END"
                  end if
               else if (item%inter(2)) then
                  if (un_perm_idx==item%idx2 .and.
     &                un_perm_name==item%label_t2) then
                     write(itflog,'(a)') "END"
                  end if
               end if
            end if

         end if
      end if


      ! If created a perm intermediate, print the symmetrised lines
      if (symmetric .and. itin) then
         call print_symmetrise(old_name,item)
      end if


      ! Mark end of spin block
      if (item%print_line) write(itflog,'(a)') "END"

      ! Deallocate memroy used when construcitng item
      call itf_deinit(item)

      if (intpp) then
         write(itflog,'(a)') "END_INTPP"
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine set_itype(item, itype)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
     &   item        ! ITF contraction object; holds all info about the ITF algo line
      integer, intent(inout) ::
     &   itype(INDEX_LEN)

      integer ::
     &   i

      ! If the current result is not an intermeidte, don't need an itype
      ! for the next line
      if (.not. item%inter(3)) then
         itype = 0
         return
      end if

      ! Get current itype of intermediate result
      !do i = 1, item%rank3
      !   itype(i) = get_itype(item%idx3(i:i))
      !end do

      if (item%rank3==2) then
         itype(1) = get_itype(item%idx3(1:1))
         itype(2) = get_itype(item%idx3(2:2))
      else if (item%rank3==4) then
         itype(1) = get_itype(item%idx3(1:1))
         itype(2) = get_itype(item%idx3(2:2))
         itype(3) = get_itype(item%idx3(4:4))
         itype(4) = get_itype(item%idx3(3:3))
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine create_permutation(item, perm)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &   item        ! ITF contraction object; holds all info about the ITF algo line
      logical, intent(in) ::
     &   perm(ngastp)

      integer ::
     &   ex_itype,
     &   i, j,
     &   shift
      character(len=1) ::
     &   ex_ind(2)
      logical ::
     &   found1,
     &   found2
      character(len=INDEX_LEN) ::
     &   tmp1, tmp2

      tmp1 = item%idx1
      tmp2 = item%idx2

      ! Identify and swap external indicies
      do i = 1, ngastp
         if (perm(i)) then
            ex_itype = i
            exit
         end if
      end do


      shift = 1
      do i = 1, item%rank3
         if (get_itype(item%idx3(i:i),.true.)/=ex_itype) cycle

         do j = 1, item%rank1
            if (item%idx3(i:i)==item%idx1(j:j)) then
               ex_ind(shift) = item%idx1(j:j)
               shift = shift + 1
            end if
         end do

         do j = 1, item%rank2
            if (item%idx3(i:i)==item%idx2(j:j)) then
               ex_ind(shift) = item%idx2(j:j)
               shift = shift + 1
            end if
         end do

      end do

      if (shift < 3) write(item%logfile,*) "ERROR"

      found1 = .false.
      found2 = .false.
      do i = 1, item%rank1
         if (tmp1(i:i)==ex_ind(1)) then
            tmp1(i:i) = ex_ind(2)
            found1 = .true.
            exit
         end if
      end do
      if (.not. found1) then
         do i = 1, item%rank1
            if (tmp1(i:i)==ex_ind(2)) then
               tmp1(i:i) = ex_ind(1)
               found1 = .true.
               exit
            end if
         end do
      end if

      do i = 1, item%rank2
         if (tmp2(i:i)==ex_ind(1)) then
            tmp2(i:i) = ex_ind(2)
            found2 = .true.
         end if
      end do
      if (.not. found2) then
         do i = 1, item%rank2
            if (tmp2(i:i)==ex_ind(2)) then
               tmp2(i:i) = ex_ind(1)
               found2 = .true.
            end if
         end do
      end if

      write(item%logfile,*) "ex ind ", ex_ind
      write(item%logfile,*) "ex type ", ex_itype
      write(item%logfile,*) "perm ", perm

      ! Permute indicies to get correct pairing
      tmp1 = f_index(tmp1, item%rank1/2)
      tmp2 = f_index(tmp2, item%rank2/2)

      write(item%logfile,*) "idx1 ", tmp1
      write(item%logfile,*) "idx2 ", tmp2

      return
      end


*----------------------------------------------------------------------*
      subroutine create_permutation2(item, perm)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
     &   item        ! ITF contraction object; holds all info about the ITF algo line
      logical, intent(in) ::
     &   perm(ngastp)

      integer ::
     &   ex_itype,
     &   i, j,
     &   shift
      character(len=1) ::
     &   ex_ind(2)
      logical ::
     &   found1,
     &   found2
      character(len=INDEX_LEN) ::
     &   tmp1, tmp2

      tmp1 = item%idx1
      tmp2 = item%idx2

      ! Identify and swap external indicies
      do i = 1, ngastp
         if (perm(i)) then
            ex_itype = i
            exit
         end if
      end do


      shift = 1
      do i = 1, item%rank3
         if (get_itype(item%idx3(i:i),.true.)/=ex_itype) cycle

         do j = 1, item%rank1
            if (item%idx3(i:i)==item%idx1(j:j)) then
               ex_ind(shift) = item%idx1(j:j)
               shift = shift + 1
            end if
         end do

         do j = 1, item%rank2
            if (item%idx3(i:i)==item%idx2(j:j)) then
               ex_ind(shift) = item%idx2(j:j)
               shift = shift + 1
            end if
         end do

      end do

      if (shift < 3) write(item%logfile,*) "ERROR"

      do i = 1, item%rank1
         if (item%idx1(i:i)==ex_ind(1)) then
            tmp1(i:i) = ex_ind(2)
            exit
         end if
      end do
      do i = 1, item%rank1
         if (item%idx1(i:i)==ex_ind(2)) then
            tmp1(i:i) = ex_ind(1)
            exit
         end if
      end do

      do i = 1, item%rank2
         if (item%idx2(i:i)==ex_ind(1)) then
            tmp2(i:i) = ex_ind(2)
            exit
         end if
      end do
      do i = 1, item%rank2
         if (item%idx2(i:i)==ex_ind(2)) then
            tmp2(i:i) = ex_ind(1)
            exit
         end if
      end do

      !write(item%logfile,*) "ex ind ", ex_ind
      !write(item%logfile,*) "ex type ", ex_itype
      !write(item%logfile,*) "perm ", perm

      ! Update orginal index
      item%idx1 = tmp1
      item%idx2 = tmp2

      ! Negative sign from (1-P_xy)
      item%fact = item%fact * -1.0d+0

      return
      end


*----------------------------------------------------------------------*
      subroutine convert_to_abab_block(item, t_spin, new_idx1, new_idx2,
     &                                 new_fact)
*----------------------------------------------------------------------*
!     Convert integrals and amplitudes to abab spn blocks, this may also
!     intoduce a sign change
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
     &   item
      type(spin_info2), intent(in) ::
     &      t_spin(3)
      character(len=INDEX_LEN), intent(inout) ::
     &     new_idx1, new_idx2
      real(8), intent(inout) ::
     &     new_fact

      integer ::
     &   i
      character(len=1) ::
     &   tmp
      integer ::
     &   ops(4,2)

      ! TODO: Don't bother for eeee or cccc integrals - a hack for now
      ! to allow the summation of K:eeee terms
      ops = item%e1 + item%c
      if (ops(2,2)==2 .and. ops(2,1)==2 .or.
     &    ops(1,2)==2 .and. ops(1,1)==2) then
         if (.not. item%inter(1) .and. .not. item%inter(3)) then
            return
         end if
      end if

      ! TODO: only work for rank 4
      if (item%rank1>2 .and. .not. item%inter(1)) then
         if (t_spin(1)%spin(1,1)>
     &       t_spin(1)%spin(1,2)) then
            tmp = new_idx1(2:2)
            new_idx1(2:2) = new_idx1(1:1)
            new_idx1(1:1) = tmp
            new_fact = new_fact * -1.0d+0
         end if
         if (t_spin(1)%spin(2,1)>
     &       t_spin(1)%spin(2,2)) then
            tmp = new_idx1(3:3)
            new_idx1(3:3) = new_idx1(4:4)
            new_idx1(4:4) = tmp
            new_fact = new_fact * -1.0d+0
         end if
      end if

      if (item%rank2>2 .and. .not. item%inter(2)) then
         if (t_spin(2)%spin(1,1)>
     &       t_spin(2)%spin(1,2)) then
            tmp = new_idx2(2:2)
            new_idx2(2:2) = new_idx2(1:1)
            new_idx2(1:1) = tmp
            new_fact = new_fact * -1.0d+0
         end if
         if (t_spin(2)%spin(2,1)>
     &       t_spin(2)%spin(2,2)) then
            tmp = new_idx2(3:3)
            new_idx2(3:3) = new_idx2(4:4)
            new_idx2(4:4) = tmp
            new_fact = new_fact * -1.0d+0
         end if
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine prepare_symmetrise(contr_info, itin, intpp,
     &                                  symmetric, command, old_name)
*----------------------------------------------------------------------*
!     If a symmetric residual and symmetrising after every term, rename
!     result tensor to ITIN. If not symmetrising after every term,
!     rename tensor to G (and symmetrise at end)
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'

      type(binary_contr), intent(inout) ::
     &   contr_info      ! Information about binary contraction
      logical, intent(in) ::
     &   itin,             ! Print ITIN lines or not
     &   intpp
      logical, intent(inout) ::
     &   symmetric
      integer, intent(in) ::
     &   command         ! Type of formula item command, ie. contraction, copy etc.
      character(len=MAXLEN_BC_LABEL), intent(inout) ::
     &   old_name

      ! Check if result is a symmetric matrix, if not, then no
      ! permuational symmetry and not extra factors
      call check_symmetric(contr_info, command, symmetric)

      ! If a symmetric residual, symmetrise after every term. Introduce
      ! ITIN intermeidate to collect terms
      ! .R[abij] += I[abij]
      ! .R[abij] += I[baji]
      if (symmetric .and. itin) then
         old_name = contr_info%label_res
         contr_info%label_res = "ITIN"
      end if

      ! If not symmetrising after every term, rename residual to G
      if (symmetric .and. .not. itin) then
         old_name = contr_info%label_res

         if (intpp) then
            contr_info%label_res = "INTpp"
         else
            contr_info%label_res = "G"
         end if
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine prepare_symmetrise2(contr_info, command, perm_case,
     &                               old_name, item)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(inout) ::
     &   contr_info      ! Information about binary contraction
      integer, intent(in) ::
     &   command         ! Type of formula item command, ie. contraction, copy etc.
      character(len=MAXLEN_BC_LABEL), intent(inout) ::
     &   old_name
      integer, intent(in) ::
     &   perm_case    ! Info of permutation factors
      type(itf_contr2), intent(inout) ::
     &   item

      logical ::
     &   symmetric

      ! Return if INTPP block
      if (item%intpp) return

      ! Return if an intermediate or result rank less than 2
      if (item%inter(3) .or. item%rank3<=2) return

      ! Check permuational symmetry, if permutational symmetry, then
      ! multiply by 0.5 (we will symmetrise each result)
      if (perm_case == 0) item%fact = item%fact * 0.5d+0

      ! Introduce ITIN intermeidate to collect terms
      ! .R[abij] += I[abij]
      ! .R[abij] += I[baji]
      old_name = item%label_res
      item%label_res = "ITIN"


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
     &     item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &    i,j,k,
     &    type3(2,INDEX_LEN),
     &    t_shift
      logical ::
     &    summed,
     &    symmetric

      do i = 1, MAXINT
         if (contr_info%label_res == spin_inters(i)%name) then
            item%itype = spin_inters(i)%itype
            item%symm_res = spin_inters(i)%symm_res
         !write(itflog,*) "multiple intermediate ", spin_inters(i)%itype
            exit
         end if
      end do

      call itf_contr_init(contr_info,item,permute,.false.,command,
     &                    itflog)

      !write(10,*) "intermediate_spin_info ", item%label_res
      !write(10,*) "intermediate_spin_info ", item%idx3

      !write(10,*) "intermediate_spin_info ", item%label_t1
      !write(10,*) "intermediate_spin_info ", item%idx1
      !write(10,*) "intermediate_spin_info ", item%label_t2
      !write(10,*) "intermediate_spin_info ", item%idx2

      ! Placed in normal ordered {} order not ITF []
      type3 = 0
      t_shift = 1
      ! If the t1 is an intermediate get info about its pairing
      if (item%inter(1)) then
         do i = 1, item%rank1

            if (i>item%rank1/2) then
               j = item%rank1 - i + item%rank1/2 + 1
            else
               j = i
            end if

            type3(t_shift,j) = get_itype(item%idx1(i:i))
         end do
         t_shift = t_shift + 1
      end if

      ! If the t2 is an intermediate get info about its pairing
      if (item%inter(2)) then
         do i = 1, item%rank2

            if (i>item%rank2/2) then
               j = item%rank2 - i + item%rank2/2 + 1
            else
               j = i
            end if

            type3(t_shift,j) = get_itype(item%idx2(i:i))
         end do
      end if

      if (item%permute == 2 .and. item%symm_res) then
         ! Need to transpose by tensors after permutation, to
         ! avoid symmetry problem when using (1 + Pabij)
         ! When we transpose tensors, we get a sign change, however
         ! here, we are just searching for intermediates and don't care
         ! about the final sign - that will come when the intermediate
         ! is printed out
         if (item%rank1==6) then
            !item%idx1 = c_index(item%idx1,1,.true.)
            item%idx1 = c_index(item%idx1,1)
         else
            item%idx1 = f_index(item%idx1,item%rank1/2)
         end if
         item%idx2 = f_index(item%idx2,item%rank2/2)

         if (item%inter(3)) then
            if (item%rank3==6) then
               !item%idx3 = c_index(item%idx3,1,.true.)
               item%idx3 = c_index(item%idx3,1)
            else
               item%idx3 = f_index(item%idx3,item%rank3/2)
            end if
         end if

         item%label_res = trim(item%label_res)//'P'
      end if


      ! Allocate space to store information about intermediates and
      ! their spin cases. Only allocate 2 objects as there can only be
      ! at most two intermediates on a line
      allocate(item%inter_spins(2))


      ! Do not want to print out the lines while gathering info about
      ! intermediates
      item%print_line = .false.


      ! Need to catch lines which don't need to be spin summed
      call assign_simple_spin(item, summed, item%symm_res)

      ! Need to set spin result for intermediate that depends on
      ! another intermediate
      if (.not. summed) then

         if (item%inter(3)) then
            do i = 1, MAXINT
               if (item%label_res == spin_inters(i)%name) then
                  do k = 1, spin_inters(i)%ncase

                     do j = 1, item%rank3/2
                        item%i_spin%spin(1,j) =
     &                                       spin_inters(i)%cases(j,k)
                        item%i_spin%spin(2,j) =
     &                             spin_inters(i)%cases(j+INDEX_LEN/2,k)
                     end do

                     ! Need to clear spins from after the first loop
                     item%t_spin(1)%spin = 0
                     item%t_spin(2)%spin = 0

                     call assign_spin(item)
                  end do
               exit

               end if
            end do
         else
            ! Result is not an intermediate, ie. its a residual, so lets find
            ! out what spin intermediates it needs.
            call assign_spin(item)
         end if
      end if

      item%print_line = .true.

      !write(item%logfile,*) "summed: ", summed, " ninter: ", item%ninter

      if (item%ninter == 0) call line_error("Couldn't find "//
     &                                      "intermediate", item)

      ! Copy information back to array in print_itf()
      do i = 1, item%ninter
         spin_inters(i+n_inter)=item%inter_spins(i)
         spin_inters(i+n_inter)%symm_res=item%symm_res
!         write(item%logfile,*) "hello ", item%label_res, " ",
!     &                        (spin_inters(i+n_inter)%name)
!         write(item%logfile,*) (spin_inters(i+n_inter)%cases)

         ! Set index type information, used to figure out pairing
         do j = 1, INDEX_LEN
            spin_inters(i+n_inter)%itype(j)=type3(i,j)
         end do
         !write(item%logfile,*) "type3 ", (type3(i,j),j=1,INDEX_LEN)
      end do


      ! Overall number of intermediates used to index spin_inters
      n_inter = n_inter + item%ninter

      deallocate(item%inter_spins)

      ! Deallocate memroy used when construcitng item
      call itf_deinit(item)

      return
      end


*----------------------------------------------------------------------*
      subroutine find_spin_intermediate(contr_info,itflog,command,
     &                                  spin_inters,n_inter,symm_res)
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
      logical, intent(in) ::
     &   symm_res

      integer ::
     &    perm_case,   ! Info of permutation factors
     &    i                ! Loop index
      logical ::
     &    inter,            ! True if result is an intermediate
     &    symmetric           ! True is a symmetric tensor
      character(len=MAXLEN_BC_LABEL) ::
     &    old_name


      ! Initialise permutation factors to 0 == no permutation
      perm_case = 0
      do i = 1, ngastp
         if(contr_info%perm(i)) perm_case = perm_case + 1
      end do

      !call check_symmetric(contr_info, command, symmetric)

      if (perm_case == 0) then
         ! No permutations
         call intermediate_spin_info(contr_info,itflog,command,
     &                            spin_inters,n_inter,1)
      else
         if (symm_res) then
            do i=1, perm_case
               call intermediate_spin_info(contr_info,itflog,command,
     &                                  spin_inters,n_inter,i)
            end do
         else
            do i=1, perm_case+1
               call intermediate_spin_info(contr_info,itflog,command,
     &                                  spin_inters,n_inter,i)
            end do
         end if
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine assign_simple_spin(item, summed, symmetric)
*----------------------------------------------------------------------*
!     Assign spin to simple cases using logic conditions
!     This avoid the spin sum routine
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &    item
      logical, intent(inout) ::
     &    summed
      logical, intent(in) ::
     &    symmetric

      integer ::
     &    i,j,k

      summed = .false.

      if (item%rank3 + item%rank1 + item%rank2 == 6) then
         j = 1
         do i = 1, 2
            if (item%inter(i)) then
               ! Only need the aa case
               item%inter_spins(j)%cases(1,1) = 1
               item%inter_spins(j)%cases(1+INDEX_LEN/2,1) = 1
               item%inter_spins(j)%ncase = 1
               if (i == 1) item%inter_spins(j)%name = item%label_t1
               if (i == 2) item%inter_spins(j)%name = item%label_t2
               item%ninter = item%ninter + 1
               j = j + 1
            end if
         end do

         !write(item%logfile,*) "SIMPLE SPIN: ", item%inter_spins
         summed = .true.

      else if (item%rank1 /= 0 .and. item%rank2 == 0) then

         j = 1
         do i = 1, 2
            if (item%inter(i)) then

               if (i == 1) then
                  select case (item%rank1)
                     case(2)
                        ! aa
                        item%inter_spins(j)%cases(1,1) = 1
                        item%inter_spins(j)%cases(1+INDEX_LEN/2,1) = 1
                     case(4)
                        ! abab
                        item%inter_spins(j)%cases(1,1) = 1
                        item%inter_spins(j)%cases(2,1) = 2
                        item%inter_spins(j)%cases(1+INDEX_LEN/2,1) = 1
                        item%inter_spins(j)%cases(2+INDEX_LEN/2,1) = 2
                     case(6)
                        do k = 1, 3
                           ! aaaaaa
                           item%inter_spins(j)%cases(i,1) = 1
                           item%inter_spins(j)%cases(i+INDEX_LEN/2,1)= 1
                        end do
                     case default
                        call line_error("Could not determine tensor "//
     &                                  "rank",item)
                  end select

                  item%inter_spins(j)%name = item%label_t1
                  item%inter_spins(j)%ncase = 1

                  if (item%permute == 2 .and. item%symm_res) then
                     item%inter_spins(j)%name =
     &                             trim(item%inter_spins(j)%name)//'P'
                     item%idx1 = f_index(item%idx1,item%rank1/2,.true.)
                     item%inter_spins(j)%cases(1,1) = 2
                     item%inter_spins(j)%cases(2,1) = 1
                     item%inter_spins(j)%cases(1+INDEX_LEN/2,1) = 2
                     item%inter_spins(j)%cases(2+INDEX_LEN/2,1) = 1
                  end if
               end if

               if (i == 2) then
                  item%inter_spins(j)%name = item%label_t2
                  item%inter_spins(j)%ncase = 0
               end if

               item%ninter = item%ninter + 1
               j = j + 1
            end if
         end do


         !write(item%logfile,*) "SIMPLE SPIN test: ", item%inter_spins
         summed = .true.

      else if (item%rank3==4.and.item%rank1==2.and.item%rank2==2) then
         ! Tensor product
         j = 1
         do i = 1, 2
            if (item%inter(i)) then
               ! Only need the aa case
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
     &                               label,spin_case,itype,ninter,
     &                               symm_res)
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
     &     spin_case(INDEX_LEN),
!     &     itype(2,INDEX_LEN),
     &     itype(INDEX_LEN),
     &   ninter ! Total number of intermediates
      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &     label
      logical, intent(in) ::
     &   symm_res

      type(itf_contr) ::
     &     item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &    i, j, k, l, m
      character(len=INDEX_LEN) ::
     &    spin_name,
     &    tspin_name
      logical ::
     &   found,
     &   upper

      ! Mark begging of spin summed block
      write(itflog,'(a5)') 'BEGIN'

      ! Set index type, which tells us the info about how the
      ! intermediates are paired
      item%itype = itype
      item%symm_res = symm_res

      ! TODO: subroutine to print out all info in itf_contr
      if (scan('P', label)) item%permutation = .true.
      call itf_contr_init(contr_info,item,0,.false.,command,itflog)

      ! Set overall spin case of result
      !write(item%logfile,*) "spin_case ", spin_case
      do i = 1, item%rank3/2
         item%i_spin%spin(1,i) = spin_case(i)
         item%i_spin%spin(2,i) = spin_case(i+INDEX_LEN/2)
      end do

      !call print_spin(item%i_spin%spin, item%rank3, "TEST", 11)

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
      if (scan('P', label) .and. item%symm_res) then

         !item%permutation = .true.
         found = .false.
!         do j = 1, ngastp
!            if (item%nops3(j) > 2) then
!
!               ! Maybe this will lead to problems..
!               ! Rationale - sometimes due to the permuation symmetry of
!               ! tensor (ie. with 3 or 4 indicies of the same type), we
!               ! can skip explicitly writting out all the intermediates
!               ! which are needed. This is like when skipping the result
!               ! (permutation) line above.
!               do k = 1, ngastp
!                  if (item%inter(1) .and. item%nops1(k)>2 .and.
!     &                item%rank2==2) then
!                     return
!                  else if (item%inter(2) .and. item%nops2(k)>2 .and.
!     &                item%rank1==2) then
!                     return
!                  end if
!               end do
!
!
!               do k = 1, len(spin_name)
!                  if (spin_name(k:k)=='a') then
!                     spin_name(k:k)='b'
!                  else if (spin_name(k:k)=='b') then
!                     spin_name(k:k)='a'
!                  end if
!               end do
!
!               !! Find the itype of the 3 indicies
!               if (j == 1) then
!                  l = 3
!               else if (j == 2) then
!                  l = 1
!               else if (j == 3) then
!                  l = 2
!               else if (j == 4) then
!                  l = 4
!               end if
!
!               do k = 1, INDEX_LEN
!                  if (item%itype(1,k)/=l) then
!                     if (k>=item%rank3/2) then
!                        upper = .false.
!                     else
!                        upper = .true.
!                     end if
!                     exit
!                  end if
!               end do
!
!               !item%idx3=f_index(item%idx3,item%rank3/2,.true.)
!               item%idx3=f_index(item%idx3,item%rank3/2,upper)
!
!               ! Also flip i_spin just in case
!               do l = 1, 2
!                  do k = 1, item%rank3/2
!                     if (item%i_spin%spin(l,k)==1) then
!                        item%i_spin%spin(l,k)=2
!                     else if (item%i_spin%spin(l,k)==2) then
!                        item%i_spin%spin(l,k)=1
!                     end if
!                  end do
!               end do
!
!
!               found = .true.
!               exit
!            end if
!         end do


         if (.not.found) then
            if (item%rank3==6) then
               item%idx3 = c_index(item%idx3,1)
               !item%idx3 = c_index(item%idx3,1,.true.)
            else
               item%idx3=f_index(item%idx3,item%rank3/2)
            end if
         end if

         !item%idx1 = f_index(item%idx1,item%rank1/2,.true.)
         if (item%rank1==6) then
            item%idx1 = c_index(item%idx1,1,.true.)
         else
            item%idx1 = f_index(item%idx1,item%rank1/2)
         end if
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

      write(itflog,'(a)') "END"

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
      real(8) ::
     &   c_fact               ! Copy of orginal factor


      ! Reorder integrals into a fixed index order, only need to do this
      ! once for each spin block as the spin summation doesn't depend on
      ! the index string at this point
      if (item%spin_cases==0) then
         call reorder_integral(item%int(1),item%rank1,item%idx1,s1,
     &                         item%j_int,
     &                         item%label_t1,item%nops1)
         call reorder_integral(item%int(2),item%rank2,item%idx2,s2,
     &                         item%j_int,
     &                         item%label_t2,item%nops2)
         if (.not. item%int(1) .and. .not. item%inter(1)) then
            call reorder_amp(item%rank1,item%idx1)
         end if
         if (.not. item%int(2) .and. .not. item%inter(2)) then
            call reorder_amp(item%rank2,item%idx2)
         end if
         if (.not. item%inter(3)) then
            call reorder_amp(item%rank3,item%idx3)
         end if
      end if

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
     &                1,item%binary,item%int(1),item%nops1,item%j_int,
     &                item%logfile)
      call spatial_string(st2,item%idx2,nt2,s2,item%inter(2),item%rank2,
     &                2,item%binary,item%int(2),item%nops2,item%j_int,
     &                item%logfile)


      ! Add factor to sclar result cases (going to skip half the spin
      ! cases as these are the same, so add a factor of two to the
      ! remaining ones)
      c_fact = item%fact
      if (item%rank3 == 0 .and. item%rank1/=0) then
         c_fact = c_fact * 2.0d0
      end if

      ! Convert factor to string, ignore if 1.0 or -1.0
      sfact=''
      sfact_star=''
      if (abs(abs(c_fact) - 1.0d+0) > 1.0d-15) then
            write(sfact,*) c_fact

            if (c_fact < 0.0d+0) then
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
      if (c_fact < 0.0d+0) then
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

      ! Increment number of printed spn cases
      item%spin_cases = item%spin_cases + 1

      return
      end


*----------------------------------------------------------------------*
      subroutine print_spin_cases(item)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
     &     item

      integer ::
     &   i, j,
     &   actual_spin_cases
      logical ::
     &   contains1, contains2,
     &   s1, s2


      ! Cover case where only 1 spin case per line
      if (item%nspin_cases==1) then
         actual_spin_cases = 2
      else
         actual_spin_cases = item%nspin_cases
      end if


      do j = 1, actual_spin_cases -1

            ! Decide if tensor is mixed spin
            contains1 = .false.
            contains2 = .false.
            s1 = .false.
            if (item%rank1>2) then
               do i = 1, size(item%all_spins(j)%t_spin(1)%spin,2)
                  if (item%all_spins(j)%t_spin(1)%spin(1,i)==1) then
                     contains1 = .true.
                  end if
                  if (item%all_spins(j)%t_spin(1)%spin(1,i)==2) then
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
               do i = 1, size(item%all_spins(j)%t_spin(2)%spin,2)
                  if (item%all_spins(j)%t_spin(2)%spin(1,i)==1) then
                     contains1 = .true.
                  end if
                  if (item%all_spins(j)%t_spin(2)%spin(1,i)==2) then
                     contains2 = .true.
                  end if
               end do
               if (contains1 .neqv. contains2) then
                  s2 = .true.
               end if
            end if

            call print_itf_line2(item,s1,s2,item%all_spins(j)%t_spin)

      end do

      return
      end

*----------------------------------------------------------------------*
      subroutine print_itf_line2(item,s1,s2,t_spin)
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

      type(itf_contr2), intent(inout) ::
     &     item
      logical, intent(in) ::
     &     s1,s2
      type(spin_info2), intent(in) ::
     &      t_spin(3)

      character(len=MAXLEN_BC_LABEL) ::
     &     nres, nt1, nt2          ! Name of tensors involved in the contraction
      character(len=INDEX_LEN) ::
     &     slabel1, slabel2, slabel3,
     &     new_idx1, new_idx2
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
      real(8) ::
     &   c_fact               ! Copy of orginal factor
      logical ::
     &   new_j


      new_idx1 = item%idx1
      new_idx2 = item%idx2
      c_fact = item%fact
      new_j = item%j_int

      ! Reorder tensor index into abab blocks
      ! May get factor change here
      call convert_to_abab_block(item, t_spin, new_idx1, new_idx2,
     &                           c_fact)

      ! Change names of specific tensors
      nres=rename_tensor(item%label_res, item%rank3)
      nt1=rename_tensor(item%label_t1, item%rank1)
      nt2=rename_tensor(item%label_t2, item%rank2)


      ! Add intermediate spin strings to names
      if (item%inter(1)) then
         call inter_spin_name2(t_spin(1)%spin,item%rank1/2,slabel1)
         nt1 = trim(nt1)//trim(slabel1)
      end if
      if (item%inter(2)) then
         call inter_spin_name2(t_spin(2)%spin,item%rank2/2,slabel2)
         nt2 = trim(nt2)//trim(slabel2)
      end if
      if (item%inter(3)) then
         call inter_spin_name2(t_spin(3)%spin,item%rank3/2,slabel3)
         nres = trim(nres)//trim(slabel3)
      end if


      ! Reorder integrals into fixed slot order
      ! TODO: I think this is ok to do after converting to abab block,
      !       but need to check...
      call reorder_integral2(item%int(1),item%rank1,new_idx1,s1,
     &                      new_j,
     &                      nt1,item%nops1)
      call reorder_integral2(item%int(2),item%rank2,new_idx2,s2,
     &                      new_j,
     &                      nt2,item%nops2)


      ! Change tensor to spatial orbital quantity, unless it is an
      ! intermediate
      call spatial_string(st1,new_idx1,nt1,s1,item%inter(1),item%rank1,
     &                1,item%binary,item%int(1),item%nops1,new_j,
     &                item%logfile)
      call spatial_string(st2,new_idx2,nt2,s2,item%inter(2),item%rank2,
     &                2,item%binary,item%int(2),item%nops2,new_j,
     &                item%logfile)


      ! Add factor to scalar result cases (going to skip half the spin
      ! cases as these are the same, so add a factor of two to the
      ! remaining ones)
      !c_fact = item%fact
      if (item%rank3 == 0 .and. item%rank1/=0) then
         c_fact = c_fact * 2.0d0
      end if

      ! Convert factor to string, ignore if 1.0 or -1.0
      sfact=''
      sfact_star=''
      if (abs(abs(c_fact) - 1.0d+0) > 1.0d-15) then
            write(sfact,*) c_fact

            if (c_fact < 0.0d+0) then
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
      if (c_fact < 0.0d+0) then
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

      ! Increment number of printed spn cases
      item%spin_cases = item%spin_cases + 1

      return
      end

*----------------------------------------------------------------------*
      subroutine reorder_amp(rank,idx)
*----------------------------------------------------------------------*
!     Print line of ITF code
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=INDEX_LEN), intent(inout) ::
     &   idx
      integer, intent(in) ::
     &   rank

      integer ::
     &   itype(rank),
     &   i,
     &   itmp
      character(len=1) ::
     &   tmp
      logical ::
     &   permute


      if (rank == 2) return

      ! TODO: Factorise this / only do it once...
      do i = 1, rank
         itype(i) = get_itype(idx(i:i))
      end do

      do i = 1, rank/2
         if (itype(i)>itype(i+rank/2)) then
            ! Swap creation and annhilation
            tmp = idx(i:i)
            idx(i:i) = idx(i+rank/2:i+rank/2)
            idx(i+rank/2:i+rank/2) = tmp

            ! Update itype
            itmp = itype(i)
            itype(i) = itype(i+rank/2)
            itype(i+rank/2) = itmp
         end if
      end do


      permute = .false.
      do i = 1, rank/2-1
         if (itype(i)>itype(i+1)) then
            tmp = idx(i:i)
            idx(i:i) = idx(i+1:i+1)
            idx(i+1:i+1) = tmp

            tmp = idx(i+rank/2:i+rank/2)
            idx(i+rank/2:i+rank/2) = idx(i+1+rank/2:i+1+rank/2)
            idx(i+1+rank/2:i+1+rank/2) = tmp

            permute = .true.
         else if (itype(i)<itype(i+1)) then
            permute = .true.
         end if
      end do

      ! Check the annhilations
      if (.not. permute) then
         do i = rank/2+1, rank-1
            if (itype(i)>itype(i+1)) then
               tmp = idx(i:i)
               idx(i:i) = idx(i+1:i+1)
               idx(i+1:i+1) = tmp

               tmp = idx(i-rank/2:i-rank/2)
               idx(i-rank/2:i-rank/2) = idx(i+1-rank/2:i+1-rank/2)
               idx(i+1-rank/2:i+1-rank/2) = tmp
            end if
         end do
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine reorder_integral(integral,rank,idx,s1,j_int,label,nops)
*----------------------------------------------------------------------*
!     Print line of ITF code
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=INDEX_LEN), intent(inout) ::
     &   idx
      integer, intent(in) ::
     &   rank,
     &   nops(ngastp)
      logical, intent(in) ::
     &   integral,
     &   s1,
     &   j_int
      character(len=MAXLEN_BC_LABEL), intent(inout) ::
     &   label

      integer ::
     &   itype(rank),
     &   i, j,
     &   itmp
      character(len=1) ::
     &   tmp
      character(len=INDEX_LEN) ::
     &   tstr
      logical ::
     &   permute,
     &   symmetric


      if (.not. integral) return
      if (rank == 2) return

      do i = 1, rank
         itype(i) = get_itype(idx(i:i))
      end do

      symmetric = .true.
      do i = 1, ngastp
         if (mod(nops(i),2) /= 0 .or. nops(i)>2) then
            symmetric = .false.
            exit
         end if
      end do

      !write(10,*) "idx ", idx
      do i = 1, rank/2
         if (itype(i)>itype(i+rank/2)) then
            ! Swap creation and annhilation
            tmp = idx(i:i)
            idx(i:i) = idx(i+rank/2:i+rank/2)
            idx(i+rank/2:i+rank/2) = tmp

            ! Update itype
            itmp = itype(i)
            itype(i) = itype(i+rank/2)
            itype(i+rank/2) = itmp
         else if (j_int) then
            if (itype(i)==itype(i+rank/2)) then
               ! Catch J-integrals: ecec -> eecc
               tmp = idx(i+1:i+1)
               idx(i+1:i+1) = idx(i+rank/2:i+rank/2)
               idx(i+rank/2:i+rank/2) = tmp
               label = 'J'

               itmp = itype(i+1)
               itype(i+1) = itype(i+rank/2)
               itype(i+rank/2) = itmp

               ! Check index in correct order
               do j = 1, rank/2
                  if (itype(j)>itype(j+rank/2)) then
                     tmp = idx(j:j)
                     idx(j:j) = idx(j+rank/2:j+rank/2)
                     idx(j+rank/2:j+rank/2) = tmp

                     itmp = itype(j)
                     itype(j) = itype(j+rank/2)
                     itype(j+rank/2) = itmp
                  end if
               end do

            end if
         end if
      end do
      !write(10,*) "idx ", idx

      ! Permute pairs
      if (nops(2)<3) then
         permute = .false.
         do i = 1, rank/2-1
            if (itype(i)>itype(i+1)) then
               tmp = idx(i:i)
               idx(i:i) = idx(i+1:i+1)
               idx(i+1:i+1) = tmp

               tmp = idx(i+rank/2:i+rank/2)
               idx(i+rank/2:i+rank/2) = idx(i+1+rank/2:i+1+rank/2)
               idx(i+1+rank/2:i+1+rank/2) = tmp

               permute = .true.
            else if (itype(i)<itype(i+1)) then
               permute = .true.
            end if
         end do

         ! Check the annhilations
         if (.not. permute) then
            do i = rank/2+1, rank-1
               if (itype(i)>itype(i+1)) then
                  tmp = idx(i:i)
                  idx(i:i) = idx(i+1:i+1)
                  idx(i+1:i+1) = tmp

                  tmp = idx(i-rank/2:i-rank/2)
                  idx(i-rank/2:i-rank/2) = idx(i+1-rank/2:i+1-rank/2)
                  idx(i+1-rank/2:i+1-rank/2) = tmp
               end if
            end do
         end if
      !write(10,*) "idx ", idx
      else if (nops(2)==3) then
         ! Need to reorder indices for 3-external integrals
         ! Slot positions 1 and 2 are now paired

         tstr=''
         do i = 1, rank/2
            if (itype(i)==1 .and. itype(i+rank/2)==1) then
               tstr(1:1) = idx(i:i)
               tstr(2:2) = idx(i+rank/2:i+rank/2)
            else
               tstr(3:3) = idx(i:i)
               tstr(4:4) = idx(i+rank/2:i+rank/2)
            end if
         end do

         idx = trim(tstr)
      end if


      if (nops(1)==3 .and. nops(3)==1) then
         ! Need to have special case of K:ccca not K:accc
         tstr = ''
         tstr(1:1) = idx(2:2)
         tstr(2:2) = idx(3:3)
         tstr(3:3) = idx(4:4)
         tstr(4:4) = idx(1:1)
         idx = tstr
      end if

      if (label=='J'.and.nops(1)==1.and.nops(2)==2.and.nops(3)==1) then
         ! Need to have special case of J:eeca not J:eeac
         tstr = ''
         tstr(1:1) = idx(2:2)
         tstr(2:2) = idx(1:1)
         tstr(3:3) = idx(4:4)
         tstr(4:4) = idx(3:3)
         idx = tstr
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine reorder_integral2(integral,rank,idx,s1,j_int,label,
     &                             nops)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=INDEX_LEN), intent(inout) ::
     &   idx
      integer, intent(in) ::
     &   rank,
     &   nops(ngastp)
      logical, intent(in) ::
     &   integral,
     &   s1,
     &   j_int
      character(len=MAXLEN_BC_LABEL), intent(inout) ::
     &   label

      integer ::
     &   itype(rank),
     &   i, j,
     &   itmp
      character(len=1) ::
     &   tmp
      character(len=INDEX_LEN) ::
     &   tstr
      logical ::
     &   permute,
     &   symmetric


      if (.not. integral) return
      if (rank == 2) return

      do i = 1, rank
         itype(i) = get_itype(idx(i:i))
      end do

      symmetric = .true.
      do i = 1, ngastp
         if (mod(nops(i),2) /= 0 .or. nops(i)>2) then
            symmetric = .false.
            exit
         end if
      end do

      !write(10,*) "idx ", idx
      do i = 1, rank/2
         if (itype(i)>itype(i+rank/2)) then
            ! Swap creation and annhilation
            tmp = idx(i:i)
            idx(i:i) = idx(i+rank/2:i+rank/2)
            idx(i+rank/2:i+rank/2) = tmp

            ! Update itype
            itmp = itype(i)
            itype(i) = itype(i+rank/2)
            itype(i+rank/2) = itmp
         else if (j_int) then
            if (itype(i)==itype(i+rank/2)) then
               ! Catch J-integrals: ecec -> eecc
               tmp = idx(i+1:i+1)
               idx(i+1:i+1) = idx(i+rank/2:i+rank/2)
               idx(i+rank/2:i+rank/2) = tmp
               label = 'J'

               itmp = itype(i+1)
               itype(i+1) = itype(i+rank/2)
               itype(i+rank/2) = itmp

               ! Check index in correct order
               do j = 1, rank/2
                  if (itype(j)>itype(j+rank/2)) then
                     tmp = idx(j:j)
                     idx(j:j) = idx(j+rank/2:j+rank/2)
                     idx(j+rank/2:j+rank/2) = tmp

                     itmp = itype(j)
                     itype(j) = itype(j+rank/2)
                     itype(j+rank/2) = itmp
                  end if
               end do

            end if
         end if
      end do
      !write(10,*) "idx ", idx

      ! Permute pairs
      if (nops(2)<3) then
         permute = .false.
         do i = 1, rank/2-1
            if (itype(i)>itype(i+1)) then
               tmp = idx(i:i)
               idx(i:i) = idx(i+1:i+1)
               idx(i+1:i+1) = tmp

               tmp = idx(i+rank/2:i+rank/2)
               idx(i+rank/2:i+rank/2) = idx(i+1+rank/2:i+1+rank/2)
               idx(i+1+rank/2:i+1+rank/2) = tmp

               permute = .true.
            else if (itype(i)<itype(i+1)) then
               permute = .true.
            end if
         end do

         ! Check the annhilations
         if (.not. permute) then
            do i = rank/2+1, rank-1
               if (itype(i)>itype(i+1)) then
                  tmp = idx(i:i)
                  idx(i:i) = idx(i+1:i+1)
                  idx(i+1:i+1) = tmp

                  tmp = idx(i-rank/2:i-rank/2)
                  idx(i-rank/2:i-rank/2) = idx(i+1-rank/2:i+1-rank/2)
                  idx(i+1-rank/2:i+1-rank/2) = tmp
               end if
            end do
         end if
      !write(10,*) "idx ", idx
      else if (nops(2)==3) then
         ! Need to reorder indices for 3-external integrals
         ! Slot positions 1 and 2 are now paired

         tstr=''
         do i = 1, rank/2
            if (itype(i)==1 .and. itype(i+rank/2)==1) then
               tstr(1:1) = idx(i:i)
               tstr(2:2) = idx(i+rank/2:i+rank/2)
            else
               tstr(3:3) = idx(i:i)
               tstr(4:4) = idx(i+rank/2:i+rank/2)
            end if
         end do

         idx = trim(tstr)
      end if


      if (nops(1)==3 .and. nops(3)==1) then
         ! Need to have special case of K:ccca not K:accc
         tstr = ''
         tstr(1:1) = idx(2:2)
         tstr(2:2) = idx(3:3)
         tstr(3:3) = idx(4:4)
         tstr(4:4) = idx(1:1)
         idx = tstr
      end if

      if (label=='J'.and.nops(1)==1.and.nops(2)==2.and.nops(3)==1) then
         ! Need to have special case of J:eeca not J:eeac
         tstr = ''
         tstr(1:1) = idx(2:2)
         tstr(2:2) = idx(1:1)
         tstr(3:3) = idx(4:4)
         tstr(4:4) = idx(3:3)
         idx = tstr
      end if

      return
      end



*----------------------------------------------------------------------*
      subroutine spatial_string(st,idx,nt,spin,inter,rank,tensor,binary,
     &                          integral,nops,j_int,lulog)
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
     &   integral,
     &   j_int
      integer, intent(in) ::
     &   rank,       ! Rank of tensor
     &   tensor,     ! T1 or T2
     &   nops(ngastp),
     &   lulog       ! Logfile

      integer ::
     &   hrank       ! Half rank

      hrank = rank / 2

      if (spin .and. .not.inter) then
         ! Pure spin
         select case (rank)
            case (4)

               if (integral) then
                  if ((nops(1)==3 .or. nops(2)==3 .or. nops(3)==3)) then
                   ! Three-somthing integral (ie. K:eccc, K:accc, ...)
                   st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &             'K'//'['//f_index(idx,hrank,.false.,.true.)//']'//')'

                  else if (j_int) then
                     if (trim(nt)=='K') then
                        ! Need (K:eecc - J:eecc)
                        st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &                     'J'//'['//f_index(idx,hrank,.true.)//']'//')'
                     else if (trim(nt)=='J') then
                        ! Need (J:eecc - K:eecc)
                        st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &                     'K'//'['//f_index(idx,hrank,.true.)//']'//')'
                     end if

                  else if ((nops(1)==1.and.nops(2)==2.and.nops(3)==1))
     &              then
                    ! K:eeac - K:eeac
                    st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &              trimal(nt)//'['//f_index(idx,hrank)//']'//')'

                  else if ((nops(1)==1.and.nops(2)==1.and.nops(3)==2))
     &              then
                     if (get_itype(idx(2:2))==2 .and.
     &                   get_itype(idx(3:3))==2) then
                        ! K:eaac - K:ecaa
                        st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &                      trimal(nt)//'['//
     &                 f_index(idx,hrank,.false.,.false.,.false.,.true.)
     &                 //']'//')'
                     else
                        st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &              trimal(nt)//'['//f_index(idx,hrank,.true.)//']'//')'
                     end if

                  else
                    st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &              trimal(nt)//'['//f_index(idx,hrank,.true.)//']'//')'
                  end if

               else
                  ! Amplitude or density
                  if ((nops(1)==2.and.nops(2)==1.and.nops(3)==1)) then
                    ! T:eacc
                    st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &              trimal(nt)//'['//f_index(idx,hrank,.true.)//']'//')'
                  else
                     st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &                  trimal(nt)//'['//f_index(idx,hrank)//']'//')'
                  end if
               end if

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
               write(lulog,*) "ERROR: Couldn't determine rank: ",
     &                        rank
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
               if (inter) then
                  st=trimal(nt)//'['//trim(idx)//']'
               else
                  st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &                 trimal(nt)//'['//f_index(idx,hrank)//']'//')'
               end if
            case default
               write(lulog,*) "ERROR: Couldn't determine rank: ",
     &                        rank
         end select
         end if
      end if

      return
      end


!*----------------------------------------------------------------------*
!      subroutine assign_add_index(contr_info,item)
!*----------------------------------------------------------------------*
!!     Simple ITF index assignment for lines that only contain a result
!!     and a tensor, ie. COPY and ADD lines. Don't need to pair indices.
!*----------------------------------------------------------------------*
!
!      use itf_utils
!      implicit none
!      include 'opdim.h'
!      include 'def_contraction.h'
!      include 'def_itf_contr.h'
!
!      type(binary_contr), intent(in) ::
!     &     contr_info   ! Information about binary contraction
!      type(itf_contr), intent(inout) ::
!     &     item         ! ITF binary contraction
!
!      integer ::
!     &     e1(ngastp,2),     ! Operator numbers of the first tensor (T1)
!     &     i            ! Loop index
!      character, dimension(4) ::
!     &     hol=(/ 'i','j','k','l' /),
!     &     par=(/ 'a','b','c','d' /)
!      character, dimension(8) ::
!     &     val=(/ 'p','q','r','s','t','u','v','w' /)
!      character(len=INDEX_LEN) ::
!     &     c1, c2, c3,
!     &     a1, a2, a3
!      character(len=INDEX_LEN), dimension(8) ::
!     &     e1_array     ! Index of operator 1
!
!      e1=item%e1
!
!      c1='        '
!      c2='        '
!      c3='        '
!      a1='        '
!      a2='        '
!      a3='        '
!
!      ! Assign e1 (external indices of t1)
!      do i=1, e1(2,1)
!          c1(i:)=par(i)
!      end do
!      e1_array(1)=c1
!      do i=1, e1(3,1)
!          c2(i:)=val(i)
!      end do
!      e1_array(2)=c2
!      do i=1, e1(1,1)
!          c3(i:)=hol(i)
!      end do
!      e1_array(3)=c3
!
!      ! Need to to be shifted to not match assignment of creations above
!      do i=1, e1(2,2)
!          a1(i:)=par(i+e1(2,1))
!      end do
!      e1_array(5)=a1
!      do i=1, e1(3,2)
!          a2(i:)=val(i+e1(3,1))
!      end do
!      e1_array(6)=a2
!      do i=1, e1(1,2)
!          a3(i:)=hol(i+e1(1,1))
!      end do
!      e1_array(7)=a3
!
!      item%idx1=trimal(e1_array(1))//trimal(e1_array(2))//
!     &          trimal(e1_array(3))//trimal(e1_array(5))//
!     &          trimal(e1_array(6))//trimal(e1_array(7))
!
!      item%idx3=item%idx1
!
!      return
!      end


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
      subroutine init_index_str(str, rank, n_cnt)
*----------------------------------------------------------------------*
!     Initalise index_str
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(inout) ::
     &   str

      integer, intent(in) ::
     &   rank,
     &   n_cnt

      allocate(str%str(rank))
      allocate(str%itype(rank))
      allocate(str%cnt_poss(n_cnt))

      return
      end

*----------------------------------------------------------------------*
      subroutine deinit_index_str(str)
*----------------------------------------------------------------------*
!     Deinitalise index_str
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(inout) ::
     &   str

      deallocate(str%str)
      deallocate(str%itype)
      deallocate(str%cnt_poss)

      return
      end


*----------------------------------------------------------------------*
      subroutine create_index_str(idx, cnt, ex, c_shift, e_shift, rank,
     &                             second)
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
      logical, intent(in) ::
     &   second

      integer ::
     &   shift_c,
     &   shift_a,
     &   cnt_shift(ngastp),
     &   i, j, k, l,
     &   shift_p,
     &   op1, op2
      character, dimension(22) ::
     &   ind=(/ 'i','j','k','l','m','n','o','a','b','c','d','e','f',
     &          'g','p','q','r','s','t','u','v','w' /)  ! Letters for index string


      if (second) then
         op1 = 2
         op2 = 1
         shift_c = rank
         shift_a = 1
      else
         op1 = 1
         op2 = 2
         shift_c = 1
         shift_a = rank
      end if

      !shift_c = 1
      !shift_a = rank
      cnt_shift = c_shift
      l = 1

      do i = 1, ngastp
         !if (cnt(i,1)>0) then
         if (cnt(i,op1)>0) then
            !do j = 1, cnt(i,1)
            do j = 1, cnt(i,op1)
               k = 1+(7*(i-1)) + cnt_shift(i)
               idx%str(shift_c) = ind(k)
               idx%cnt_poss(l) = shift_c
               idx%itype(shift_c) = i
               l = l + 1
               cnt_shift(i) = cnt_shift(i) + 1
               if (second) then
                  shift_c = shift_c - 1
               else
                  shift_c = shift_c + 1
               end if
            end do
         end if
         !if (ex(i,1)>0) then
         if (ex(i,op1)>0) then
            !do j = 1, ex(i,1)
            do j = 1, ex(i,op1)
               k = 1+(7*(i-1)) + e_shift(i)
               idx%str(shift_c) = ind(k)
               idx%itype(shift_c) = i
               e_shift(i) = e_shift(i) + 1
               !shift_c = shift_c + 1
               if (second) then
                  shift_c = shift_c - 1
               else
                  shift_c = shift_c + 1
               end if
            end do
         end if

         !if (cnt(i,2)>0) then
         if (cnt(i,op2)>0) then
            !do j = 1, cnt(i,2)
            do j = 1, cnt(i,op2)
               k = 1+(7*(i-1)) + cnt_shift(i)
               idx%str(shift_a) = ind(k)
               idx%cnt_poss(l) = shift_a
               idx%itype(shift_a) = i
               l = l + 1
               cnt_shift(i) = cnt_shift(i) + 1
               if (second) then
                  shift_a = shift_a + 1
               else
                  shift_a = shift_a - 1
               end if
            end do
         end if
         !if (ex(i,2)>0) then
         if (ex(i,op2)>0) then
            !do j = 1, ex(i,2)
            do j = 1, ex(i,op2)
               k = 1+(7*(i-1)) + e_shift(i)
               idx%str(shift_a) = ind(k)
               idx%itype(shift_a) = i
               e_shift(i) = e_shift(i) + 1
               !shift_a = shift_a - 1
               if (second) then
                  shift_a = shift_a + 1
               else
                  shift_a = shift_a - 1
               end if
            end do
         end if
      end do

      ! Change i_type (index type) values to canonical values
      do i = 1, rank
         select case (idx%itype(i))
            case (1)
               idx%itype(i) = 3
            case (2)
               idx%itype(i) = 1
            case (3)
               idx%itype(i) = 2
            case (4)
               idx%itype(i) = 4
         end select
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine nicer_pairing(str1, str2, i1, i2, rank1, rank2, factor,
     &                         n_cnt)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(inout) ::
     &   str1,
     &   str2
      character (len=1), intent(in) ::
     &   i1,         ! Creation index
     &   i2          ! Annhilation index
      integer, intent(in) ::
     &   rank1,
     &   rank2,
     &   n_cnt
      real(8), intent(inout) ::
     &   factor

      character (len=1) ::
     &   tmp
      integer ::
     &   pp1, pp2, pp3,        ! Pair position
     &   distance,
     &   itmp,
     &   i, j, k, l, m, n

      do i = 1, rank1
         if (i1==str1%str(i)) then
            pp1 = rank1-i+1
            do k = 1, rank2
               if (i2==str2%str(k)) then
                  pp2= rank2-k+1
                  if (str1%str(pp1)/=str2%str(pp2)) then
                     do l = 1, rank2
                        if(str1%str(pp1)==str2%str(l))
     &                      then

                  !write(10,*) "str1%str(i) ", str1%str(i)
                  !write(10,*) "str1%str(pp) ", str1%str(pp1)
                  !write(10,*) "str2%str(k) ", str2%str(l)
                  !write(10,*) "str2%str(pp) ", str2%str(pp2)

                           pp3 = rank2-l+1
                           tmp = str2%str(k)
                           str2%str(k) = str2%str(pp3)
                           str2%str(pp3) = tmp

                           !tmp = str2%str(l)
                           !str2%str(l) = str2%str(pp2)
                           !str2%str(pp2) = tmp

                           !write(10,*) "New ", str2%str

                           ! Update factor
!                           distance = abs(l-pp2)
!                           if (mod(distance,2)/=0) then
!!                              factor = factor * -1.0d0
!!!                              write(11,*)"Update factor3.2: ",
!!!     &                                              factor
!                              write(11,*) "not updating factor"
!                           end if


                           ! If swapping a contraction index
                           do m = 1, n_cnt
                              if (pp3==str2%cnt_poss(m)) then
!                                 write(11,*) "Old cnt_poss ",
!     &                                                str2%cnt_poss
                                 str2%cnt_poss(m) = k
!                                 write(11,*) "New cnt_poss ",
!     &                                                str2%cnt_poss
                                 exit
                              end if
                           end do

!                           ! If swapping a contraction index, need to
!                           ! update cnt_poss
!!                           write(11,*) "Old cnt_poss ", str2%cnt_poss
!                           do m = 1, n_cnt
!                              if (l==str2%cnt_poss(m)) then
!                                itmp=str2%cnt_poss(m)
!                                str2%cnt_poss(m)=pp2
!                                n = m
!                              end if
!                           end do
!!                           write(11,*) "Old cnt_poss2 ", str2%cnt_poss
!                           do m = 1, n_cnt
!                              if (rank2-k+1==
!     &                           str2%cnt_poss(m) .and. m/=n) then
!                                 str2%cnt_poss(m)=itmp
!                              end if
!                           end do
!!                           write(11,*) "New cnt_poss ", str2%cnt_poss

                           exit
                        end if
                     end do
                     exit
                  end if
                  exit
               end if
            end do
            exit
         end if
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine nicer_pairing_one_t(str, i1, i2, rank, factor, n_cnt)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(inout) ::
     &   str
      character (len=1), intent(in) ::
     &   i1,
     &   i2
      integer, intent(in) ::
     &   rank,
     &   n_cnt
      real(8), intent(inout) ::
     &   factor

      character (len=1) ::
     &   tmp
      integer ::
     &   pp,
     &   distance,
     &   i, j, k, l, m, n

      do i = 1, rank
         if (i1==str%str(i)) then
            pp = rank-i+1
            do k = 1, rank
               if (i2==str%str(k)) then
                  if (pp/=k) then
                     ! Swap index so it is in paired position (pp)
                     !write(item%logfile,*) "Old str ", str%str
                     tmp = str%str(pp)
                     str%str(pp) = str%str(k)
                     str%str(k) = tmp
                     !write(item%logfile,*) "New str ", str%str

!                     ! Update factor
!                     distance = abs(k-pp)
!                     if (mod(distance,2)/=0) then
!                        !factor = factor * -1.0d0
!                        !write(11,*)"Update factor3.5: ", factor
!                        write(11,*) "Not updating factor"
!                     end if

                     ! If swapping a contraction index, need to
                     ! update cnt_poss
                     do l = 1, n_cnt
                        if (pp==str%cnt_poss(l)) then
!                           write(item%logfile,*) "Old cnt_poss ",
!     &                                          str1%cnt_poss
                           str%cnt_poss(l) = k
!                           write(item%logfile,*) "New cnt_poss ",
!     &                                          str1%cnt_poss
                           exit
                        end if
                     end do
                  end if
                  exit
               end if
            end do
            exit
         end if
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine assign_new_index2(item)
*----------------------------------------------------------------------*
!     Assign an ITF index string to each tensor in a line
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
     &   item           ! ITF binary contraction

      integer ::
     &   c(ngastp,2),        ! Operator numbers of contraction index
     &   ci(ngastp,2),       ! Operator numbers of contraction index (inverse)
     &   e1(ngastp,2),       ! Operator numbers of external index 1
     &   e2(ngastp,2),       ! Operator numbers of external index 2
     &   shift,         ! List shift
     &   shift_a,       ! List shift
     &   shift_c,       ! List shift
     &   n_cnt,         ! Number of contraction operators
     &   rank2, rank1,      ! Rank of tensor - either can be T1 or T2
     &   tensor2, tensor1,  ! Labels a tesor - either 1 or 2
     &   e1ops, e2ops,  ! Number of external ops on T1 and T2
     &   place,         ! Marks which tensor an index was found
     &   distance,      ! Distance from where an index should be
     &   pp,            ! Paired position - position of paired index
     &   ntest = 000,     ! >100 toggles some debug
     &   i, j, k, l, m, z,   ! Loop index
     &   itmp, n,
     &   ipair,
     &   sh1, sh2, sh3,
     &   pp1, pp2, pl1, pl2,
     &   itype(INDEX_LEN),
     &   itype1, itype2, ptype
      character(len=INDEX_LEN) ::
     &   s1, s2, s3,  ! Tmp ITF index strings
     &   tstr
      character(len=1) ::
     &   tmp, tmp2      ! Scratch space to store index letter
      real(8) ::
     &   factor,        ! Factor from equivalent lines
     &   p_factor       ! Overall factor from contraction/rearrangment
      type(pair_list) ::
     &   p_list,       ! Complete list of pairs in binary contraction
     &   p_list2        ! Complete list of pairs in binary contraction
      integer, dimension(4) ::
     &   t_shift,       ! Index shift for external indices
     &   e_shift,       ! Index shift for external indices
     &   c_shift        ! Index shift for contraction indices
      type(index_str) ::
     &   str1,          ! Stores 'normal ordered' index string for T1
     &   str2,          ! Stores 'normal ordered' index string for T2
     &   str3,          ! Stores 'normal ordered' index string for Res
     &   t_str1, t_str2 ! Scratch space

      logical ::        ! These are used when finding pairs of external ops
     &   is_cnt,        ! True if the operator is a contraction op
     &   found_end,     ! True if found external op to make a pair
     &   found_match,   ! True if a matching contraction op has been found on the opposite tensor
     &   found_cnt,     ! True if a contraction operator has been found
     &   p1, p2,
     &   found,
     &   swapped        ! True is exchanged external lines in permutation line
      integer, pointer ::
     &   t1cnt_poss(:) => null(),
     &   t2cnt_poss(:) => null()

      ! Set operator numbers
      c=item%c
      e1=item%e1
      e2=item%e2

      ! Factor due to permuation of annhilation and creation indices
      p_factor = 1.0d0

      ! Set number of contraction indicies
      n_cnt = item%contri

      ! Allocate index_str objects
      call init_index_str(str1, item%rank1, n_cnt)
      call init_index_str(str2, item%rank2, n_cnt)
      call init_index_str(str3, item%rank3, n_cnt)

      allocate(p_list%plist(item%rank3/2))

      ! Set letter shift values for contraction indices
      do i = 1, 4
         c_shift(i) = e1(i,1) + e1(i,2) + e2(i,1) + e2(i,2)
      end do

      ! Make 'tranpose' of c array, this corresponds to the operators on
      ! the second tensor
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

      ! Create an index string for T1 and T2
      ! e_shift updates after each call - this indexs the external
      ! operators; c_shift indexes the contraction operators. These are
      ! needed to ensure the letters are different
      e_shift = 0
      call create_index_str(str1,c,e1, c_shift, e_shift, item%rank1,
     &                      .false.)
      call create_index_str(str2,ci,e2,c_shift, e_shift, item%rank2,
     &                      .true.)

      ! If t1 or t2 is an intermediate, index needs to be consisent with
      ! the item%itype info (which came from the previous line)
      if (item%inter(1)) call match_idx_with_itype(item, str1, n_cnt)


      if (ntest>100) then
         write(item%logfile,*) "STR1: {", str1%str, "}"
         write(item%logfile,*) "STR2: {", str2%str, "}"
         write(item%logfile,*) "STR3: {", str1%str, "}{", str2%str, "}"
         write(item%logfile,*) "CNT POSS1: ", str1%cnt_poss
         write(item%logfile,*) "CNT POSS2: ", str2%cnt_poss
      end if


      ! Find external indicies and correctly pair them
      ! Make a copy of itype which will be modified when a pair is found
      itype = item%itype
      e1ops = sum(sum(e1,dim=2))
      e2ops = sum(sum(e2,dim=2))
      if (e1ops >= e2ops) then
        call find_pairs_wrap2(str1,str2,item%rank1,item%rank2,1,2,n_cnt,
     &                       item,p_list,itype)
      else
        call find_pairs_wrap2(str2,str1,item%rank2,item%rank1,2,1,n_cnt,
     &                       item,p_list,itype)
      end if
      !call print_plist(p_list, item%rank3/2, "PAIRS", item%logfile)


      ! If there is a pair in one string, permute so they are paired
      ! 'imediately'. Also update the position of the contraction index
      do j = 1, item%rank3/2
         p1 = .false.
         p2 = .false.

         if (p_list%plist(j)%ops(1)==1) p1 = .true.
         if (p_list%plist(j)%ops(2)==1) p2 = .true.

         ! If a pair is on one tensor...
         if (p1 .and. p2 .and. item%rank1>2) then
            call nicer_pairing_one_t(str1, p_list%plist(j)%pindex(1),
     &                               p_list%plist(j)%pindex(2),
     &                               item%rank1, p_factor, n_cnt)
         else if (.not. p1 .and. .not. p2 .and. item%rank2>2) then
            call nicer_pairing_one_t(str2, p_list%plist(j)%pindex(1),
     &                               p_list%plist(j)%pindex(2),
     &                               item%rank2, p_factor, n_cnt)

         ! If a pair is split over two tensors...
         else if (p1 .and. .not. p2 .and. item%rank1>2) then
            call nicer_pairing(str1, str2,
     &                         p_list%plist(j)%pindex(1),
     &                         p_list%plist(j)%pindex(2), item%rank1,
     &                         item%rank2, p_factor, n_cnt)
         else if (p2 .and. .not. p1 .and. item%rank2>2) then
            call nicer_pairing(str2, str1,
     &                         p_list%plist(j)%pindex(1),
     &                         p_list%plist(j)%pindex(2), item%rank2,
     &                         item%rank1, p_factor, n_cnt)
         end if
      end do
      !call print_plist(p_list,item%rank3/2,"NICER PAIRS",item%logfile)


      ! Work out the factor due to permuation of contraction indicies
      if (.not. item%den(1) .and. .not. item%den(2)) then
         do i = 1, n_cnt
          if (mod(item%rank1-str1%cnt_poss(i)+str2%cnt_poss(i)-1,2)
     &        /=0)then
             p_factor = p_factor * -1.0d0
             if (ntest>100) then
               write(item%logfile,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
             end if
          end if
         end do

      else if (item%den(1) .and. item%rank1/=0) then
!         write(item%logfile,*) "str poss1 ", str1%cnt_poss, str1%str
!         write(item%logfile,*) "str poss2 ", str2%cnt_poss, str2%str

         do i = 1, n_cnt
            if (str1%cnt_poss(i)<=item%rank1/2) then
               ! Creation operator of density is contracted
               if (mod(item%rank1-str1%cnt_poss(i)+str2%cnt_poss(i)-1
     &                 -item%rank1/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                write(item%logfile,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
                  end if
               end if

            else if (str1%cnt_poss(i)>item%rank1/2) then
               ! Annhilation operator of density is contracted
               ! NOTE: not removing previously conracted indicies from
               ! the list - so this may be wrong in the long term. Seems
               ! to produce the correct answer for now...
               if (mod(item%rank2-str2%cnt_poss(i)+str1%cnt_poss(i)-1
     &                 -item%rank1/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                write(item%logfile,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
                  end if
               end if

            end if
         end do

      else if (item%den(2) .and. item%rank2/=0) then
         !write(item%logfile,*) "str poss1 ", str1%cnt_poss, str1%str
         !write(item%logfile,*) "str poss2 ", str2%cnt_poss, str2%str

         do i = 1, n_cnt
            if (str2%cnt_poss(i)<=item%rank2/2) then
               ! Creation operator of density is contracted
               if (mod(item%rank2-str2%cnt_poss(i)+str1%cnt_poss(i)-1
     &                 -item%rank2/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                write(item%logfile,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
                  end if
               end if

            else if (str2%cnt_poss(i)>item%rank2/2) then
               ! Annhilation operator of density is contracted
               if (mod(item%rank1-str1%cnt_poss(i)+str2%cnt_poss(i)-1
     &                 -item%rank2/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                write(item%logfile,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
                  end if
               end if

            end if
         end do

      end if


      ! Create result index string from only external operators. Order
      ! is not final...This is basically splicing str1 and str2 together
      shift = 1
      do i = 1, item%rank1
         is_cnt = .false.
         do j = 1, n_cnt
            if (i == str1%cnt_poss(j)) then
               is_cnt = .true.
               exit
            end if
         end do
         if (.not. is_cnt) then
            str3%str(shift) = str1%str(i)
            shift = shift + 1
         end if
      end do

      do i = 1, item%rank2
         is_cnt = .false.
         do j = 1, n_cnt
            if (i == str2%cnt_poss(j)) then
               is_cnt = .true.
               exit
            end if
         end do
         if (.not. is_cnt) then
            str3%str(shift) = str2%str(i)
            shift = shift + 1
         end if
      end do

      !write(item%logfile,*) "T1 string: {", str1%str, "}", str1%cnt_poss
      !write(item%logfile,*) "T2 string: {", str2%str, "}", str2%cnt_poss
      !write(item%logfile,*) "Result string: {", str3%str, "}"

      ! Rearrange the result string so it is in normal order (all
      ! creation operators to the left of the annhilation). This can
      ! also introduce a sign change.
      do j = 1, item%rank3/2
         shift = 1
         !do i = 0, item%rank3/2-1
         do i = 0, item%rank3/2
          if (p_list%plist(j)%pindex(1) == str3%str(item%rank3-i)) then

               tstr(shift:shift) = str3%str(item%rank3-i)
               shift = shift + 1

               do k = 1, item%rank3
                  if (p_list%plist(j)%pindex(1) /=
     &                                      str3%str(k)) then
                     tstr(shift:shift) = str3%str(k)
                     shift = shift + 1
                  end if
               end do
               do k = 1, item%rank3
                  str3%str(k) = tstr(k:k)
               end do

               !write(item%logfile,*) "New result string {",str3%str, "}"

               ! Update factor. If index is in even position, requires
               ! odd number of permuations; so get a negative
               if (mod(item%rank3-i,2)==0) then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                     write(item%logfile,*) "Update factor 2 (rearrange",
     &               " the result string to normal order): ", p_factor
                  end if
               end if

               exit
          end if
         end do
      end do

      ! Move all annhilations to the right
      do j = 1, item%rank3/2
         shift = item%rank3
         do i = 1, item%rank3/2
          if (p_list%plist(j)%pindex(2) == str3%str(i)) then

               tstr(shift:shift) = str3%str(i)
               shift = shift - 1

               do k = 0, item%rank3-1
                  if (p_list%plist(j)%pindex(2) /=
     &                                      str3%str(item%rank3-k)) then
                     tstr(shift:shift) = str3%str(item%rank3-k)
                     shift = shift - 1
                  end if
               end do
               do k = 1, item%rank3
                  str3%str(k) = tstr(k:k)
               end do

               !write(item%logfile,*) "New result string {",str3%str, "}"

               ! Update factor. If index is in odd position, requires
               ! odd number of permuations; so get a negative
               if (mod(i,2)/=0) then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                     write(item%logfile,*) "Update factor 3 (move
     &               annhilations to the right in the result string): ",
     &               p_factor
                  end if
               end if

               exit
          end if
         end do
      end do

      ! Due to how the R:ea residual is defined, {p a^+} instead of
      ! {a^+ p}, we need an extra minus to flip the normal ordered
      ! string.
      !if (.not. item%inter(3)) then
         if (item%nops3(1)==0 .and. item%nops3(2)==1 .and.
     &       item%nops3(3)==1 .and. item%nops3(4)==0) then
         if (.not. item%inter(3)) then

            p_factor = p_factor * -1.0d0
            if (ntest>100) then
               write(item%logfile,*) "Update factor (R:ea)", p_factor
            end if
         else
            !write(item%logfile,*) "kidder "
         end if
         end if
      !end if

      ! Rearrange result operators to get correct pairing
      do j = 1, item%rank3/2
         do i = 0, item%rank3/2-1

          ! Search for annhilation ops
          if (p_list%plist(j)%pindex(2) == str3%str(item%rank3-i)) then
!            write(item%logfile,*) "Found ", str3%str(item%rank3-i)
            do k = 1, item%rank3/2

               ! Search for creation ops
               if (p_list%plist(j)%pindex(1) ==
     &                                      str3%str(k)) then
                  ! Check if position is in the pair position
!                  write(item%logfile,*) "Found ", str3%str(k), " in ", k
!                  write(item%logfile,*) "Pair ", str3%str(item%rank3-i),
!     &                                  " in ",item%rank3-i
!                  write(item%logfile,*) "PP ",
!     &                                  item%rank3-(item%rank3-i)+1
                  pp = item%rank3-(item%rank3-i)+1

                  if (item%rank3-k+1 /= item%rank3-i) then
                     ! Get distance from paired position (pp) and factor
                     distance = abs(k-pp)
                     if (mod(distance,2)/=0) then
                        p_factor = p_factor * -1.0d0
                        if (ntest>100) then
                           write(item%logfile,*)"Update factor 4 ",
     &                      "(rearrange reuslt string to correct ",
     &                      "pairing): ", p_factor
                        end if
                     end if
                     ! Swapped letters between current k and pp
                     tmp = str3%str(pp)
                     str3%str(pp) = str3%str(k)
                     str3%str(k) = tmp

!                     write(11,*) "swapping ", str3%str(pp),
!     &                         " and ", str3%str(k)
!
!                     write(11,*) "Swap ops string {", str3%str, "}"
                  end if
                  exit
               end if
            end do
          end if
         end do
      end do


      ! Permute strings into nicer order
      ! Get canocial values for result string
      do i = 1, item%rank3
         str3%itype(i) = get_itype(str3%str(i))
      end do

      ! Permute string: {baij} => {abji}, {ia} => {ai} etc.
      call permute_index(str1, item%rank1)
      call permute_index(str2, item%rank2)
      if (item%inter(3)) then
         ! For intermeidiates, the index pairing and order must match how
         ! it is used (it must match the itype string). Therefore the pair
         ! ordering needs to be checked. This doesn't introduce a sign
         ! change as pairs of indices are swapped.
         !write(item%logfile,*) "Result string{", str3%str, "}"
         !write(item%logfile,*) "Result itype{", str3%itype, "}"
         !write(item%logfile,*) "Real itype{", item%itype, "}"
         !call permute_slot_order(str3, item%rank3, item%itype)
         call permute_index(str3, item%rank3)
      else
         call permute_index(str3, item%rank3)
      end if


      ! Check intermediates have correct itypes, if not, permute the
      ! pairs
      if (item%inter(1)) call match_idx_with_itype2(item, str1, n_cnt)


      s1 = ""
      s2 = ""
      s3 = ""

      do i = 1, item%rank1/2
         s1(i:i) = str1%str(i)
         s1(item%rank1/2+i:item%rank1/2+i)=str1%str(item%rank1-(i-1))
      end do
      do i = 1, item%rank2/2
         s2(i:i) = str2%str(i)
         s2(item%rank2/2+i:item%rank2/2+i)=str2%str(item%rank2-(i-1))
      end do
      do i = 1, item%rank3/2
         s3(i:i) = str3%str(i)
         s3(item%rank3/2+i:item%rank3/2+i)=str3%str(item%rank3-(i-1))
      end do

      if (ntest>10) then
         if (p_factor>0.0d0) then
            tmp = '+'
         else
            tmp = '-'
         end if

         write(item%logfile,*)
         write(item%logfile,*) "---------------------------"
         write(item%logfile,*) "{",str3%str,"} ",tmp,"= {",str1%str,
     &                         "}{",str2%str,"}"

         write(item%logfile,*) "[",trim(s3),"] ",tmp,"= [",trim(s1),
     &                         "][",trim(s2),"]"
         write(item%logfile,*) "---------------------------"
      end if

      item%idx1=trim(s1)
      item%idx2=trim(s2)
      item%idx3=trim(s3)
      item%fact = item%fact * p_factor

      deallocate(p_list%plist)

      call deinit_index_str(str1)
      call deinit_index_str(str2)
      call deinit_index_str(str3)

      return
      end


*----------------------------------------------------------------------*
      subroutine match_idx_with_itype(item, idx, n_cnt)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
     &   item           ! ITF binary contraction
      type(index_str), intent(inout) ::
     &   idx
      integer, intent(in) ::
     &   n_cnt

      character(len=1) ::
     &   tmp
      integer ::
     &   itmp, i, cnt_poss(n_cnt)


      if (item%rank1<=2) return
      ! TODO: this is shit
      !write(item%logfile,*) "item: ", item%itype
      !write(item%logfile,*) "idx: ", idx%itype
      !write(item%logfile,*) "cnt_poss: ", idx%cnt_poss
      !write(item%logfile,*) "str: ", idx%str

      cnt_poss = idx%cnt_poss

      if (item%itype(1)/=idx%itype(1)) then
         tmp = idx%str(2)
         idx%str(2) = idx%str(1)
         idx%str(1) = tmp
         itmp = idx%itype(2)
         idx%itype(2) = idx%itype(1)
         idx%itype(1) = itmp

         do i = 1, n_cnt
            if (idx%cnt_poss(i)==1) then
               cnt_poss(i)=2
               exit
            end if
         end do
         do i = 1, n_cnt
            if (idx%cnt_poss(i)==2) then
               cnt_poss(i)=1
               exit
            end if
         end do

      end if

      if (item%rank1>2) then
         if (item%itype(3)/=idx%itype(3)) then
            tmp = idx%str(4)
            idx%str(4) = idx%str(3)
            idx%str(3) = tmp
            itmp = idx%itype(4)
            idx%itype(4) = idx%itype(3)
            idx%itype(3) = itmp
         end if

         do i = 1, n_cnt
            if (idx%cnt_poss(i)==3) then
               cnt_poss(i)=4
               exit
            end if
         end do
         do i = 1, n_cnt
            if (idx%cnt_poss(i)==4) then
               cnt_poss(i)=3
               exit
            end if
         end do
      end if

      idx%cnt_poss = cnt_poss

      !write(item%logfile,*) "item: ", item%itype
      !write(item%logfile,*) "idx: ", idx%itype
      !write(item%logfile,*) "cnt_poss: ", idx%cnt_poss
      !write(item%logfile,*) "str: ", idx%str

      return
      end


*----------------------------------------------------------------------*
      subroutine match_idx_with_itype2(item, idx, n_cnt)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
     &   item           ! ITF binary contraction
      type(index_str), intent(inout) ::
     &   idx
      integer, intent(in) ::
     &   n_cnt

      character(len=1) ::
     &   tmp1, tmp2
      integer ::
     &   itmp1, itmp2, i, cnt_poss(n_cnt)


      if (item%rank1<=2) return

      ! Update itype which may have changed
      do i = 1, item%rank1
         idx%itype(i) = get_itype(idx%str(i))
      end do

      !write(item%logfile,*) "itype: ", item%itype
      !write(item%logfile,*) "idx itype: ", idx%itype
      !write(item%logfile,*) "idx: ", idx%str
      do i = 1, 4
         if (item%itype(i)/=idx%itype(i)) then

            !write(item%logfile,*) "changing ", idx%str
            !write(item%logfile,*) "changing ", idx%itype

            tmp1 = idx%str(2)
            tmp2 = idx%str(4)
            idx%str(2) = idx%str(1)
            idx%str(4) = idx%str(3)
            idx%str(1) = tmp1
            idx%str(3) = tmp2

            itmp1 = idx%itype(2)
            itmp2 = idx%itype(4)
            idx%itype(2) = idx%itype(1)
            idx%itype(4) = idx%itype(3)
            idx%itype(1) = itmp1
            idx%itype(3) = itmp2

            !write(item%logfile,*) "to ", idx%str
            !write(item%logfile,*) "to ", idx%itype

         end if
      end do

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
     &   c(ngastp,2),        ! Operator numbers of contraction index
     &   ci(ngastp,2),       ! Operator numbers of contraction index (inverse)
     &   e1(ngastp,2),       ! Operator numbers of external index 1
     &   e2(ngastp,2),       ! Operator numbers of external index 2
     &   shift,         ! List shift
     &   shift_a,       ! List shift
     &   shift_c,       ! List shift
     &   n_cnt,         ! Number of contraction operators
     &   rank2, rank1,      ! Rank of tensor - either can be T1 or T2
     &   tensor2, tensor1,  ! Labels a tesor - either 1 or 2
     &   e1ops, e2ops,  ! Number of external ops on T1 and T2
     &   place,         ! Marks which tensor an index was found
     &   distance,      ! Distance from where an index should be
     &   pp,            ! Paired position - position of paired index
     &   ntest = 000,     ! >100 toggles some debug
     &   i, j, k, l, m, z,   ! Loop index
     &   itmp, n,
     &   ipair,
     &   sh1, sh2, sh3,
     &   pp1, pp2, pl1, pl2,
     &   itype(INDEX_LEN),
     &   itype1, itype2, ptype
      character(len=INDEX_LEN) ::
     &   s1, s2, s3,  ! Tmp ITF index strings
     &   tstr
      character(len=1) ::
     &   tmp, tmp2      ! Scratch space to store index letter
      real(8) ::
     &   factor,        ! Factor from equivalent lines
     &   p_factor       ! Overall factor from contraction/rearrangment
      type(pair_list) ::
     &   p_list,       ! Complete list of pairs in binary contraction
     &   p_list2        ! Complete list of pairs in binary contraction
      integer, dimension(4) ::
     &   t_shift,       ! Index shift for external indices
     &   e_shift,       ! Index shift for external indices
     &   c_shift        ! Index shift for contraction indices
      type(index_str) ::
     &   str1,          ! Stores 'normal ordered' index string for T1
     &   str2,          ! Stores 'normal ordered' index string for T2
     &   str3,          ! Stores 'normal ordered' index string for Res
     &   t_str1, t_str2 ! Scratch space

      logical ::        ! These are used when finding pairs of external ops
     &   is_cnt,        ! True if the operator is a contraction op
     &   found_end,     ! True if found external op to make a pair
     &   found_match,   ! True if a matching contraction op has been found on the opposite tensor
     &   found_cnt,     ! True if a contraction operator has been found
     &   p1, p2,
     &   found,
     &   swapped        ! True is exchanged external lines in permutation line
      integer, pointer ::
     &   t1cnt_poss(:) => null(),
     &   t2cnt_poss(:) => null()

      ! Set operator numbers
      c=item%c
      e1=item%e1
      e2=item%e2

      ! Factor due to permuation of annhilation and creation indices
      p_factor = 1.0d0

      ! Set number of contraction indicies
      n_cnt = item%contri

      ! Allocate index_str objects
      call init_index_str(str1, item%rank1, n_cnt)
      call init_index_str(str2, item%rank2, n_cnt)
      call init_index_str(str3, item%rank3, n_cnt)

      allocate(p_list%plist(item%rank3/2))

      ! Set letter shift values for contraction indices
      do i = 1, 4
         c_shift(i) = e1(i,1) + e1(i,2) + e2(i,1) + e2(i,2)
      end do

      ! Make 'tranpose' of c array, this corresponds to the operators on
      ! the second tensor
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

      ! Create an index string for T1 and T2
      ! e_shift updates after each call - this indexs the external
      ! operators; c_shift indexes the contraction operators. These are
      ! needed to ensure the letters are different
      e_shift = 0
      call create_index_str(str1,c,e1, c_shift, e_shift, item%rank1,
     &                      .false.)
      call create_index_str(str2,ci,e2,c_shift, e_shift, item%rank2,
     &                      .true.)

      if (ntest>100) then
         write(item%logfile,*) "STR1: {", str1%str, "}"
         write(item%logfile,*) "STR2: {", str2%str, "}"
         write(item%logfile,*) "STR3: {", str1%str, "}{", str2%str, "}"
         write(item%logfile,*) "CNT POSS1: ", str1%cnt_poss
         write(item%logfile,*) "CNT POSS2: ", str2%cnt_poss
      end if


      ! Find external indicies and correctly pair them
      ! Make a copy of itype which will be modified when a pair is found
      itype = item%itype
      e1ops = sum(sum(e1,dim=2))
      e2ops = sum(sum(e2,dim=2))
      if (e1ops >= e2ops) then
        call find_pairs_wrap(str1,str2,item%rank1,item%rank2,1,2,n_cnt,
     &                       item,p_list,itype)
      else
        call find_pairs_wrap(str2,str1,item%rank2,item%rank1,2,1,n_cnt,
     &                       item,p_list,itype)
      end if
      !call print_plist(p_list, item%rank3/2, "PAIRS", item%logfile)


      ! If there is a pair in one string, permute so they are paired
      ! 'imediately'. Also update the position of the contraction index
      do j = 1, item%rank3/2
         p1 = .false.
         p2 = .false.

         if (p_list%plist(j)%ops(1)==1) p1 = .true.
         if (p_list%plist(j)%ops(2)==1) p2 = .true.

         ! If a pair is on one tensor...
         if (p1 .and. p2 .and. item%rank1>2) then
            call nicer_pairing_one_t(str1, p_list%plist(j)%pindex(1),
     &                               p_list%plist(j)%pindex(2),
     &                               item%rank1, p_factor, n_cnt)
         else if (.not. p1 .and. .not. p2 .and. item%rank2>2) then
            call nicer_pairing_one_t(str2, p_list%plist(j)%pindex(1),
     &                               p_list%plist(j)%pindex(2),
     &                               item%rank2, p_factor, n_cnt)

         ! If a pair is split over two tensors...
         else if (p1 .and. .not. p2 .and. item%rank1>2) then
            call nicer_pairing(str1, str2,
     &                         p_list%plist(j)%pindex(1),
     &                         p_list%plist(j)%pindex(2), item%rank1,
     &                         item%rank2, p_factor, n_cnt)
         else if (p2 .and. .not. p1 .and. item%rank2>2) then
            call nicer_pairing(str2, str1,
     &                         p_list%plist(j)%pindex(1),
     &                         p_list%plist(j)%pindex(2), item%rank2,
     &                         item%rank1, p_factor, n_cnt)
         end if
      end do
      !call print_plist(p_list,item%rank3/2,"NICER PAIRS",item%logfile)


      ! Work out the factor due to permuation of contraction indicies
      if (.not. item%den(1) .and. .not. item%den(2)) then
         do i = 1, n_cnt
          if (mod(item%rank1-str1%cnt_poss(i)+str2%cnt_poss(i)-1,2)
     &        /=0)then
             p_factor = p_factor * -1.0d0
             if (ntest>100) then
               write(item%logfile,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
             end if
          end if
         end do

      else if (item%den(1) .and. item%rank1/=0) then
!         write(item%logfile,*) "str poss1 ", str1%cnt_poss, str1%str
!         write(item%logfile,*) "str poss2 ", str2%cnt_poss, str2%str

         do i = 1, n_cnt
            if (str1%cnt_poss(i)<=item%rank1/2) then
               ! Creation operator of density is contracted
               if (mod(item%rank1-str1%cnt_poss(i)+str2%cnt_poss(i)-1
     &                 -item%rank1/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                write(item%logfile,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
                  end if
               end if

            else if (str1%cnt_poss(i)>item%rank1/2) then
               ! Annhilation operator of density is contracted
               ! NOTE: not removing previously conracted indicies from
               ! the list - so this may be wrong in the long term. Seems
               ! to produce the correct answer for now...
               if (mod(item%rank2-str2%cnt_poss(i)+str1%cnt_poss(i)-1
     &                 -item%rank1/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                write(item%logfile,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
                  end if
               end if

            end if
         end do

      else if (item%den(2) .and. item%rank2/=0) then
         !write(item%logfile,*) "str poss1 ", str1%cnt_poss, str1%str
         !write(item%logfile,*) "str poss2 ", str2%cnt_poss, str2%str

         do i = 1, n_cnt
            if (str2%cnt_poss(i)<=item%rank2/2) then
               ! Creation operator of density is contracted
               if (mod(item%rank2-str2%cnt_poss(i)+str1%cnt_poss(i)-1
     &                 -item%rank2/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                write(item%logfile,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
                  end if
               end if

            else if (str2%cnt_poss(i)>item%rank2/2) then
               ! Annhilation operator of density is contracted
               if (mod(item%rank1-str1%cnt_poss(i)+str2%cnt_poss(i)-1
     &                 -item%rank2/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                write(item%logfile,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
                  end if
               end if

            end if
         end do

      end if



      ! Create result index string from only external operators. Order
      ! is not final...This is basically splicing str1 and str2 together
      shift = 1
      do i = 1, item%rank1
         is_cnt = .false.
         do j = 1, n_cnt
            if (i == str1%cnt_poss(j)) then
               is_cnt = .true.
               exit
            end if
         end do
         if (.not. is_cnt) then
            str3%str(shift) = str1%str(i)
            shift = shift + 1
         end if
      end do

      do i = 1, item%rank2
         is_cnt = .false.
         do j = 1, n_cnt
            if (i == str2%cnt_poss(j)) then
               is_cnt = .true.
               exit
            end if
         end do
         if (.not. is_cnt) then
            str3%str(shift) = str2%str(i)
            shift = shift + 1
         end if
      end do

      !write(item%logfile,*) "T1 string: {", str1%str, "}", str1%cnt_poss
      !write(item%logfile,*) "T2 string: {", str2%str, "}", str2%cnt_poss
      !write(item%logfile,*) "Result string: {", str3%str, "}"

      ! Rearrange the result string so it is in normal order (all
      ! creation operators to the left of the annhilation). This can
      ! also introduce a sign change.
      do j = 1, item%rank3/2
         shift = 1
         !do i = 0, item%rank3/2-1
         do i = 0, item%rank3/2
          if (p_list%plist(j)%pindex(1) == str3%str(item%rank3-i)) then

               tstr(shift:shift) = str3%str(item%rank3-i)
               shift = shift + 1

               do k = 1, item%rank3
                  if (p_list%plist(j)%pindex(1) /=
     &                                      str3%str(k)) then
                     tstr(shift:shift) = str3%str(k)
                     shift = shift + 1
                  end if
               end do
               do k = 1, item%rank3
                  str3%str(k) = tstr(k:k)
               end do

               !write(item%logfile,*) "New result string {",str3%str, "}"

               ! Update factor. If index is in even position, requires
               ! odd number of permuations; so get a negative
               if (mod(item%rank3-i,2)==0) then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                     write(item%logfile,*) "Update factor 2 (rearrange",
     &               " the result string to normal order): ", p_factor
                  end if
               end if

               exit
          end if
         end do
      end do

      ! Move all annhilations to the right
      do j = 1, item%rank3/2
         shift = item%rank3
         do i = 1, item%rank3/2
          if (p_list%plist(j)%pindex(2) == str3%str(i)) then

               tstr(shift:shift) = str3%str(i)
               shift = shift - 1

               do k = 0, item%rank3-1
                  if (p_list%plist(j)%pindex(2) /=
     &                                      str3%str(item%rank3-k)) then
                     tstr(shift:shift) = str3%str(item%rank3-k)
                     shift = shift - 1
                  end if
               end do
               do k = 1, item%rank3
                  str3%str(k) = tstr(k:k)
               end do

               !write(item%logfile,*) "New result string {",str3%str, "}"

               ! Update factor. If index is in odd position, requires
               ! odd number of permuations; so get a negative
               if (mod(i,2)/=0) then
                  p_factor = p_factor * -1.0d0
                  if (ntest>100) then
                     write(item%logfile,*) "Update factor 3 (move
     &               annhilations to the right in the result string): ",
     &               p_factor
                  end if
               end if

               exit
          end if
         end do
      end do

      ! Due to how the R:ea residual is defined, {p a^+} instead of
      ! {a^+ p}, we need an extra minus to flip the normal ordered
      ! string.
      !if (.not. item%inter(3)) then
         if (item%nops3(1)==0 .and. item%nops3(2)==1 .and.
     &       item%nops3(3)==1 .and. item%nops3(4)==0) then
         if (.not. item%inter(3)) then

            p_factor = p_factor * -1.0d0
            if (ntest>100) then
               write(item%logfile,*) "Update factor (R:ea)", p_factor
            end if
         else
            !write(item%logfile,*) "kidder "
         end if
         end if
      !end if

      ! Rearrange result operators to get correct pairing
      do j = 1, item%rank3/2
         do i = 0, item%rank3/2-1

          ! Search for annhilation ops
          if (p_list%plist(j)%pindex(2) == str3%str(item%rank3-i)) then
!            write(item%logfile,*) "Found ", str3%str(item%rank3-i)
            do k = 1, item%rank3/2

               ! Search for creation ops
               if (p_list%plist(j)%pindex(1) ==
     &                                      str3%str(k)) then
                  ! Check if position is in the pair position
!                  write(item%logfile,*) "Found ", str3%str(k), " in ", k
!                  write(item%logfile,*) "Pair ", str3%str(item%rank3-i),
!     &                                  " in ",item%rank3-i
!                  write(item%logfile,*) "PP ",
!     &                                  item%rank3-(item%rank3-i)+1
                  pp = item%rank3-(item%rank3-i)+1

                  if (item%rank3-k+1 /= item%rank3-i) then
                     ! Get distance from paired position (pp) and factor
                     distance = abs(k-pp)
                     if (mod(distance,2)/=0) then
                        p_factor = p_factor * -1.0d0
                        if (ntest>100) then
                           write(item%logfile,*)"Update factor 4 ",
     &                      "(rearrange reuslt string to correct ",
     &                      "pairing): ", p_factor
                        end if
                     end if
                     ! Swapped letters between current k and pp
                     tmp = str3%str(pp)
                     str3%str(pp) = str3%str(k)
                     str3%str(k) = tmp

!                     write(11,*) "swapping ", str3%str(pp),
!     &                         " and ", str3%str(k)
!
!                     write(11,*) "Swap ops string {", str3%str, "}"
                  end if
                  exit
               end if
            end do
          end if
         end do
      end do


      ! Permute strings into nicer order
      ! Get canocial values for result string
      do i = 1, item%rank3
         str3%itype(i) = get_itype(str3%str(i))
      end do

      ! Permute string: {baij} => {abji}, {ia} => {ai} etc.
      call permute_index(str1, item%rank1)
      call permute_index(str2, item%rank2)
      if (item%inter(3)) then
         ! For intermeidiates, the index pairing and order must match how
         ! it is used (it must match the itype string). Therefore the pair
         ! ordering needs to be checked. This doesn't introduce a sign
         ! change as pairs of indices are swapped.
         !write(item%logfile,*) "Result string{", str3%str, "}"
         !write(item%logfile,*) "Result itype{", str3%itype, "}"
         !write(item%logfile,*) "Real itype{", item%itype, "}"
         call permute_slot_order(str3, item%rank3, item%itype)
      else
         call permute_index(str3, item%rank3)
      end if


      ! Need to swap annhilation ops (1-P_ij)
      ! Don't swap index when there is a tensor product
      ! TODO: this will not work if the result is greater than rank 4
      if (item%permute==2 .and. .not. item%product
     &    .and. .not. item%inter(3) .and. item%symm_res) then
!     &    .and. .not. item%inter(3)) then
         !write(item%logfile,*) "hello ", item%label_res
         !write(item%logfile,*) "hello ", item%label_t1
         !write(item%logfile,*) "hello ", item%label_t2
         !write(item%logfile,*) "str1 ", str1%str
         !write(item%logfile,*) "str2 ", str2%str

         ! TODO: This information is not being handed to the lines
         !if (item%rank1/=2) then
            tstr=''
            do i = 1, item%rank1
               if (str1%str(i)==p_list%plist(1)%pindex(2)) then
                  tstr(i:i) = p_list%plist(2)%pindex(2)
               else if (str1%str(i)==p_list%plist(2)%pindex(2)) then
                  tstr(i:i) = p_list%plist(1)%pindex(2)
               else
                  tstr(i:i) = str1%str(i)
               end if
            end do
            do i = 1, item%rank1
               str1%str(i) = tstr(i:i)
            end do
         !end if

         !if (item%rank2/=2) then
            tstr=''
            do i = 1, item%rank2
               if (str2%str(i)==p_list%plist(1)%pindex(2)) then
                  tstr(i:i) = p_list%plist(2)%pindex(2)
               else if (str2%str(i)==p_list%plist(2)%pindex(2)) then
                  tstr(i:i) = p_list%plist(1)%pindex(2)
               else
                  tstr(i:i) = str2%str(i)
               end if
            end do

            do i = 1, item%rank2
               str2%str(i) = tstr(i:i)
            end do
         !end if
         !write(item%logfile,*) "new str1 ", str1%str
         !write(item%logfile,*) "new str2 ", str2%str
      end if


      ! Need permuatation line for non-symmetric residuals.
      ! But keep same spin
      if (item%permutation .and. .not. item%symm_res) then

         do i = 1, ngastp
            if (item%perm_case(i)) then
               ptype = i
               exit
            end if
         end do

         do i = 1, item%rank1
            ! Loop through str1 until not a contraction
            if (any(str1%cnt_poss==i)) cycle

            itype1 = get_itype(str1%str(i), .true.)

            ! Check itype matches perm_case
            if (ptype==itype1) then
               do j = 1, item%rank2
                  ! Loop through str2 until not a contraction
                  if (any(str2%cnt_poss==j)) cycle

                  itype2 = get_itype(str2%str(j), .true.)

                  ! Check external lines of same itype, and swap
                  if (itype1 == itype2) then
                     tmp = str1%str(i)
                     str1%str(i) = str2%str(j)
                     str2%str(j) = tmp
                     swapped = .true.
                  end if

                  if (swapped) exit
               end do

            end if

            ! Update factor and transpose index
            if (swapped) then
               ! Factor from (1-P_xy)
               p_factor = p_factor * -1.0d0

               ! Need to transpose indicies to maintin correct pairing
               ! TODO: have an f_string???
               if (item%rank1>2) then
                  tmp = str1%str(1)
                  str1%str(1) = str1%str(2)
                  str1%str(2) = tmp
               end if

               if (item%rank2>2) then
                  tmp = str2%str(1)
                  str2%str(1) = str2%str(2)
                  str2%str(2) = tmp
               end if

               ! Update factor due to transposition
               if (item%rank1>2) p_factor = p_factor * -1.0d0
               if (item%rank2>2) p_factor = p_factor * -1.0d0

               exit
            end if
         end do

      end if


      s1 = ""
      s2 = ""
      s3 = ""

      do i = 1, item%rank1/2
         s1(i:i) = str1%str(i)
         s1(item%rank1/2+i:item%rank1/2+i)=str1%str(item%rank1-(i-1))
      end do
      do i = 1, item%rank2/2
         s2(i:i) = str2%str(i)
         s2(item%rank2/2+i:item%rank2/2+i)=str2%str(item%rank2-(i-1))
      end do
      do i = 1, item%rank3/2
         s3(i:i) = str3%str(i)
         s3(item%rank3/2+i:item%rank3/2+i)=str3%str(item%rank3-(i-1))
      end do

      if (ntest>10) then
         if (p_factor>0.0d0) then
            tmp = '+'
         else
            tmp = '-'
         end if

         write(item%logfile,*)
         write(item%logfile,*) "---------------------------"
         write(item%logfile,*) "{",str3%str,"} ",tmp,"= {",str1%str,
     &                         "}{",str2%str,"}"

         write(item%logfile,*) "[",trim(s3),"] ",tmp,"= [",trim(s1),
     &                         "][",trim(s2),"]"
         write(item%logfile,*) "---------------------------"
      end if

      item%idx1=trim(s1)
      item%idx2=trim(s2)
      item%idx3=trim(s3)
      item%fact = item%fact * p_factor

      deallocate(p_list%plist)

      call deinit_index_str(str1)
      call deinit_index_str(str2)
      call deinit_index_str(str3)

      return
      end


*----------------------------------------------------------------------*
      subroutine find_pairs_wrap(str1, str2, rank1, rank2, t1,
     &                           t2, n_cnt, item, p_list, itype)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(in) ::
     &   str1,
     &   str2
      type(itf_contr), intent(in) ::
     &   item
      integer, intent(in) ::
     &   rank1,
     &   rank2,
     &   t1,
     &   t2,
     &   n_cnt
      integer, intent(inout) ::
     &   itype(INDEX_LEN)
      type(pair_list), intent(inout) ::
     &   p_list

      integer ::
     &   shift,
     &   i

      shift = 1

      call find_pairs_new(str1, str2, rank1, rank2, t1, t2, n_cnt,
     &                    item, shift, p_list, itype)

      ! If all the external pairs havn't been found, then there are two
      ! seperate pairs on each tensor. Go back and find them...
      if (shift-1/=item%rank3/2) then
      call find_pairs_new(str2, str1, rank2, rank1, t2, t1, n_cnt,
     &                    item, shift, p_list, itype)
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine find_pairs_new(str1, str2, rank1, rank2, t1, t2, n_cnt,
     &                          item, shift, p_list, itype)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(in) ::
     &   str1,
     &   str2
      type(itf_contr), intent(in) ::
     &   item
      integer, intent(in) ::
     &   rank1,
     &   rank2,
     &   t1,
     &   t2,
     &   n_cnt
      integer, intent(inout) ::
     &   itype(INDEX_LEN)
      integer, intent(inout) ::
     &   shift
      type(pair_list), intent(inout) ::
     &   p_list

      logical ::
     &   is_cnt,
     &   found_ex,
     &   already_found,
     &   correct_pair
      integer ::
     &   i,j,k,l,
     &   ntest = 000

      ! Search only the creations of the first string for ex ops
      do i = 1, rank1/2

         ! Check if creation has already been paired
         already_found = .false.
         do j = 1, shift
            if (str1%str(i)==p_list%plist(j)%pindex(1)) then
               already_found = .true.
               exit
            end if
         end do
         if (already_found) cycle

         if (already_found) cycle
         is_cnt = .false.
         do j = 1, n_cnt
            if (i==str1%cnt_poss(j)) then
               is_cnt = .true.
            end if
         end do

         if (.not. is_cnt) then

            ! Search the annhilations of the first string
            do j = rank1, rank1/2+1, -1

               if (ntest>0) then
               write(item%logfile,*) "searching with ", str1%str(i), t1
               write(item%logfile,*) "matching with ", str1%str(j), t1
               end if

               call suitable_pair(found_ex, str1, str1, rank1, rank1,
     &                            i, j, 2, shift, n_cnt,
     &                            p_list, item, itype)

               if (found_ex) then
                  !write(10,*)"found pair 1 ",str1%str(i)," ",str1%str(j)
                  p_list%plist(shift)%pindex(1)=str1%str(i)
                  p_list%plist(shift)%pindex(2)=str1%str(j)
                  p_list%plist(shift)%ops(1)=t1
                  p_list%plist(shift)%ops(2)=t1
                  shift = shift + 1
                  exit
               end if
            end do

            ! If it didn't find an operator in the first annhilations,
            ! look on the second annhilation ops
            if (.not. found_ex) then
               do j = rank2, rank2/2+1, -1

                  if (ntest>0) then
              write(item%logfile,*) "searching 2 with ", str1%str(i), t1
              write(item%logfile,*) "matching 2 with ", str2%str(j), t2
                  end if
                  call suitable_pair(found_ex, str1, str2, rank1, rank2,
     &                            i, j, 2, shift, n_cnt,
     &                            p_list, item, itype)

                  if (found_ex) then
                  !write(10,*)"found pair 2 ",str1%str(i)," ",str2%str(j)
                     p_list%plist(shift)%pindex(1)=str1%str(i)
                     p_list%plist(shift)%pindex(2)=str2%str(j)
                     p_list%plist(shift)%ops(1)=t1
                     p_list%plist(shift)%ops(2)=t2
         !call print_plist(p_list, shift, "p_list", item%logfile)
                     shift = shift + 1
                     exit
                  end if
               end do
            end if


         end if

         if (.not. is_cnt .and. .not. found_ex) then
            write(item%logfile,*) "Failed to find creation/annhilation",
     &                            " pair 1"
            exit
         end if

      end do


      ! Havn't found all external pairs yet
      if (shift-1/=item%rank3/2) then
         !write(11,*) "Continuing the search"

      ! Search the annhilation operators of the first string
      ! TODO: factorise this
      do i = rank1, rank1/2+1, -1

         ! Check if annhilation has already been paired
         already_found = .false.
         do j = 1, shift
            if (str1%str(i)==p_list%plist(j)%pindex(2)) then
               !write(11,*) "already found!"
               already_found = .true.
               exit
            end if
         end do
         if (already_found) cycle

         is_cnt = .false.
         do j = 1, n_cnt
            if (i==str1%cnt_poss(j)) then
               is_cnt = .true.
            end if
         end do

         if (.not. is_cnt) then

            ! Search the creations of the first string
            do j = 1, rank1/2

               if (ntest>0) then
               write(item%logfile,*) "searching 3 with ",str1%str(i), t1
               write(item%logfile,*) "matching 3 with ", str1%str(j), t1
               end if

               call suitable_pair(found_ex, str1, str1, rank1, rank1,
     &                            i, j, 1, shift, n_cnt,
     &                            p_list, item, itype)

               if (found_ex) then
                  !write(10,*)"found pair 3 ",str1%str(i)," ",str1%str(j)
                  p_list%plist(shift)%pindex(2)=str1%str(i)
                  p_list%plist(shift)%pindex(1)=str1%str(j)
                  p_list%plist(shift)%ops(2)=t1
                  p_list%plist(shift)%ops(1)=t1
                  shift = shift + 1
                  exit
               end if
            end do

            ! If it didn't find an operator in the first creations,
            ! look on the second creations ops
            if (.not. found_ex) then

               do j = 1, rank2/2

               if (ntest>0) then
              write(item%logfile,*) "searching 4 with ", str1%str(i), t1
              write(item%logfile,*) "matching 4 with ", str2%str(j), t2
               end if

                  call suitable_pair(found_ex, str1, str2, rank1, rank2,
     &                               i, j, 1, shift, n_cnt,
     &                               p_list, item, itype)

                  if (found_ex) then
                  !write(10,*)"found pair 4 ",str1%str(i)," ",str2%str(j)
                     p_list%plist(shift)%pindex(2)=str1%str(i)
                     p_list%plist(shift)%pindex(1)=str2%str(j)
                     p_list%plist(shift)%ops(2)=t1
                     p_list%plist(shift)%ops(1)=t2
                     shift = shift + 1
                     exit
                  end if
               end do
            end if


         end if

         if (.not. is_cnt .and. .not. found_ex) then
            write(item%logfile,*) "Failed to find creation/annhilation",
     &                            " pair 2"
            exit
         end if

      end do
      end if

      return
      end



*----------------------------------------------------------------------*
      subroutine suitable_pair(found_ex, str1, str2, rank1, rank2,
     &                         place1, place2, ann_cre, shift, n_cnt,
     &                         p_list, item, itype)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      logical, intent(inout) ::
     &   found_ex
      type(itf_contr), intent(in) ::
     &   item
      type(index_str), intent(in) ::
     &   str1, str2
      type(pair_list), intent(in) ::
     &   p_list
      integer, intent(in) ::
     &   place1, place2,
     &   rank1, rank2,
     &   ann_cre,
     &   shift,
     &   n_cnt
      integer, intent(inout) ::
     &   itype(INDEX_LEN)

      logical ::
     &   correct_pair
      integer ::
     &   k


      ! False if the matching index is a contraction index
      found_ex = .true.
      do k = 1, n_cnt
         if (place2==str2%cnt_poss(k)) then
            found_ex = .false.
            !write(item%logfile,*) "false here 1"
            return
         end if
      end do

      ! Check if the annhilation operator has already been paried
      do k = 1, shift
         if (str2%str(place2)==p_list%plist(k)%pindex(ann_cre)) then
            found_ex = .false.
            !write(item%logfile,*) "false here 2"
            return
         end if
      end do

      ! False if the matching index is a correct pair
      if (item%inter(3)) then
         correct_pair = .false.
         call check_pairing(correct_pair,str1,str2,rank1,
     &                      rank2,place1,place2,item, itype)
         if (.not. correct_pair) then
            found_ex = .false.
            !write(item%logfile,*) "false here 3"
            return
         end if
      end if


      return
      end


*----------------------------------------------------------------------*
      subroutine find_pairs_wrap2(str1, str2, rank1, rank2, t1,
     &                           t2, n_cnt, item, p_list, itype)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(in) ::
     &   str1,
     &   str2
      type(itf_contr), intent(in) ::
     &   item
      integer, intent(in) ::
     &   rank1,
     &   rank2,
     &   t1,
     &   t2,
     &   n_cnt
      integer, intent(inout) ::
     &   itype(INDEX_LEN)
      type(pair_list), intent(inout) ::
     &   p_list

      integer ::
     &   shift,
     &   i

      shift = 1

      call find_pairs_new2(str1, str2, rank1, rank2, t1, t2, n_cnt,
     &                    item, shift, p_list, itype)

      ! If all the external pairs havn't been found, then there are two
      ! seperate pairs on each tensor. Go back and find them...
      if (shift-1/=item%rank3/2) then
      call find_pairs_new2(str2, str1, rank2, rank1, t2, t1, n_cnt,
     &                    item, shift, p_list, itype)
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine find_pairs_new2(str1,str2, rank1, rank2, t1, t2, n_cnt,
     &                          item, shift, p_list, itype)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(in) ::
     &   str1,
     &   str2
      type(itf_contr2), intent(in) ::
     &   item
      integer, intent(in) ::
     &   rank1,
     &   rank2,
     &   t1,
     &   t2,
     &   n_cnt
      integer, intent(inout) ::
     &   itype(INDEX_LEN)
      integer, intent(inout) ::
     &   shift
      type(pair_list), intent(inout) ::
     &   p_list

      logical ::
     &   is_cnt,
     &   found_ex,
     &   already_found,
     &   correct_pair
      integer ::
     &   i,j,k,l,
     &   ntest = 000

      ! Search only the creations of the first string for ex ops
      do i = 1, rank1/2

         ! Check if creation has already been paired
         already_found = .false.
         do j = 1, shift
            if (str1%str(i)==p_list%plist(j)%pindex(1)) then
               already_found = .true.
               exit
            end if
         end do
         if (already_found) cycle

         if (already_found) cycle
         is_cnt = .false.
         do j = 1, n_cnt
            if (i==str1%cnt_poss(j)) then
               is_cnt = .true.
            end if
         end do

         if (.not. is_cnt) then

            ! Search the annhilations of the first string
            do j = rank1, rank1/2+1, -1

               if (ntest>0) then
               write(item%logfile,*) "searching with ", str1%str(i), t1
               write(item%logfile,*) "matching with ", str1%str(j), t1
               end if

               call suitable_pair3(found_ex, str1, str1, rank1, rank1,
     &                            i, j, 2, shift, n_cnt,
     &                            p_list, item, itype, t1, t2)

               if (found_ex) then
        !write(item%logfile,*)"found pair 1 ",str1%str(i)," ",str1%str(j)
                  p_list%plist(shift)%pindex(1)=str1%str(i)
                  p_list%plist(shift)%pindex(2)=str1%str(j)
                  p_list%plist(shift)%ops(1)=t1
                  p_list%plist(shift)%ops(2)=t1
                  shift = shift + 1
                  exit
               end if
            end do

            ! If it didn't find an operator in the first annhilations,
            ! look on the second annhilation ops
            if (.not. found_ex) then
               do j = rank2, rank2/2+1, -1

                  if (ntest>0) then
              write(item%logfile,*) "searching 2 with ", str1%str(i), t1
              write(item%logfile,*) "matching 2 with ", str2%str(j), t2
                  end if
                  call suitable_pair2(found_ex,str1, str2, rank1, rank2,
     &                            i, j, 2, shift, n_cnt,
     &                            p_list, item, itype, t1, t2)

                  if (found_ex) then
        !write(item%logfile,*)"found pair 2 ",str1%str(i)," ",str2%str(j)
                     p_list%plist(shift)%pindex(1)=str1%str(i)
                     p_list%plist(shift)%pindex(2)=str2%str(j)
                     p_list%plist(shift)%ops(1)=t1
                     p_list%plist(shift)%ops(2)=t2
         !call print_plist(p_list, shift, "p_list", item%logfile)
                     shift = shift + 1
                     exit
                  end if
               end do
            end if


         end if

         if (.not. is_cnt .and. .not. found_ex) then
            write(item%logfile,*) "Failed to find creation/annhilation",
     &                            " pair 1"
            exit
         end if

      end do


      ! Havn't found all external pairs yet
      if (shift-1/=item%rank3/2) then
         !write(11,*) "Continuing the search"

      ! Search the annhilation operators of the first string
      ! TODO: factorise this
      do i = rank1, rank1/2+1, -1

         ! Check if annhilation has already been paired
         already_found = .false.
         do j = 1, shift
            if (str1%str(i)==p_list%plist(j)%pindex(2)) then
               !write(11,*) "already found!"
               already_found = .true.
               exit
            end if
         end do
         if (already_found) cycle

         is_cnt = .false.
         do j = 1, n_cnt
            if (i==str1%cnt_poss(j)) then
               is_cnt = .true.
            end if
         end do

         if (.not. is_cnt) then

            ! Search the creations of the first string
            do j = 1, rank1/2

               if (ntest>0) then
               write(item%logfile,*) "searching 3 with ",str1%str(i), t1
               write(item%logfile,*) "matching 3 with ", str1%str(j), t1
               end if

               call suitable_pair2(found_ex, str1, str1, rank1, rank1,
     &                            i, j, 1, shift, n_cnt,
     &                            p_list, item, itype, t1, t2)

               if (found_ex) then
                  !write(10,*)"found pair 3 ",str1%str(i)," ",str1%str(j)
                  p_list%plist(shift)%pindex(2)=str1%str(i)
                  p_list%plist(shift)%pindex(1)=str1%str(j)
                  p_list%plist(shift)%ops(2)=t1
                  p_list%plist(shift)%ops(1)=t1
                  shift = shift + 1
                  exit
               end if
            end do

            ! If it didn't find an operator in the first creations,
            ! look on the second creations ops
            if (.not. found_ex) then

               do j = 1, rank2/2

               if (ntest>0) then
              write(item%logfile,*) "searching 4 with ", str1%str(i), t1
              write(item%logfile,*) "matching 4 with ", str2%str(j), t2
               end if

                  call suitable_pair2(found_ex,str1, str2, rank1, rank2,
     &                               i, j, 1, shift, n_cnt,
     &                               p_list, item, itype, t1, t2)

                  if (found_ex) then
                  !write(10,*)"found pair 4 ",str1%str(i)," ",str2%str(j)
                     p_list%plist(shift)%pindex(2)=str1%str(i)
                     p_list%plist(shift)%pindex(1)=str2%str(j)
                     p_list%plist(shift)%ops(2)=t1
                     p_list%plist(shift)%ops(1)=t2
                     shift = shift + 1
                     exit
                  end if
               end do
            end if


         end if

         if (.not. is_cnt .and. .not. found_ex) then
            write(item%logfile,*) "Failed to find creation/annhilation",
     &                            " pair 2"
            exit
         end if

      end do
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine suitable_pair2(found_ex, str1, str2, rank1, rank2,
     &                         place1, place2, ann_cre, shift, n_cnt,
     &                         p_list, item, itype, t1, t2)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      logical, intent(inout) ::
     &   found_ex
      type(itf_contr2), intent(in) ::
     &   item
      type(index_str), intent(in) ::
     &   str1, str2
      type(pair_list), intent(in) ::
     &   p_list
      integer, intent(in) ::
     &   place1, place2,
     &   rank1, rank2,
     &   ann_cre,
     &   shift,
     &   n_cnt,
     &   t1, t2
      integer, intent(inout) ::
     &   itype(INDEX_LEN)

      logical ::
     &   correct_pair
      integer ::
     &   k


      ! False if the matching index is a contraction index
      found_ex = .true.
      do k = 1, n_cnt
         if (place2==str2%cnt_poss(k)) then
            found_ex = .false.
            !write(item%logfile,*) "false here 1"
            return
         end if
      end do

      ! Check if the annhilation operator has already been paried
      do k = 1, shift
         if (str2%str(place2)==p_list%plist(k)%pindex(ann_cre)) then
            found_ex = .false.
            !write(item%logfile,*) "false here 2"
            return
         end if
      end do

!      ! False if the matching index is a correct pair
!      if (item%inter(1) .and. t1==1 .and. t2==1) then
!         write(item%logfile,*) "t1, t2: ", t1, t2
!         correct_pair = .false.
!         call check_pairing2(correct_pair,str1,str2,rank1,
!     &                      rank2,place1,place2,item, itype)
!         if (.not. correct_pair) then
!            found_ex = .false.
!            !write(item%logfile,*) "false here 3"
!            return
!         end if
!      end if


      return
      end


*----------------------------------------------------------------------*
      subroutine suitable_pair3(found_ex, str1, str2, rank1, rank2,
     &                         place1, place2, ann_cre, shift, n_cnt,
     &                         p_list, item, itype, t1, t2)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      logical, intent(inout) ::
     &   found_ex
      type(itf_contr2), intent(in) ::
     &   item
      type(index_str), intent(in) ::
     &   str1, str2
      type(pair_list), intent(in) ::
     &   p_list
      integer, intent(in) ::
     &   place1, place2,
     &   rank1, rank2,
     &   ann_cre,
     &   shift,
     &   n_cnt,
     &   t1, t2
      integer, intent(inout) ::
     &   itype(INDEX_LEN)

      logical ::
     &   correct_pair
      integer ::
     &   k


      found_ex = .true.

      ! False if the matching index is a contraction index
      found_ex = .true.
      do k = 1, n_cnt
         if (place2==str2%cnt_poss(k)) then
            found_ex = .false.
            !write(item%logfile,*) "false here 1"
            return
         end if
      end do

      ! Check if the annhilation operator has already been paried
      do k = 1, shift
         if (str2%str(place2)==p_list%plist(k)%pindex(ann_cre)) then
            found_ex = .false.
            !write(item%logfile,*) "false here 2"
            return
         end if
      end do

      ! False if the matching index is a correct pair
      if (item%inter(1) .and. item%rank1>2 .and. t1/=2) then
         ! TODO: t2 are useless
         correct_pair = .false.
         call check_pairing2(correct_pair,str1,str2,rank1,
     &                      rank2,place1,place2,item, itype)
         if (.not. correct_pair) then
            found_ex = .false.
            !write(item%logfile,*) "false here 3"
            return
         end if
      end if


      return
      end


*----------------------------------------------------------------------*
      subroutine check_pairing2(correct_pair, str1, str2, rank1, rank2,
     &                         place1, place2, item, itype)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      logical, intent(inout) ::
     &   correct_pair
      type(itf_contr2), intent(in) ::
     &   item
      type(index_str), intent(in) ::
     &   str1, str2
      integer, intent(in) ::
     &   place1, place2,
     &   rank1, rank2
      integer, intent(inout) ::
     &   itype(INDEX_LEN)

      integer ::
     &   i, j, extent, k,
     &   pp2

      ! Remeber - the itype info refers to the first tensor positions

      !write(item%logfile,*) "is it a correct pairing?"
      !write(item%logfile,*) "place1 ", place1
      !write(item%logfile,*) "place2 ", place2

      !write(item%logfile, *) "itype: ", itype

      if (place1>rank1/2) then
         j = item%rank3/2+1
         extent = item%rank3
      else
         j = 1
         extent = item%rank3/2
      end if

      do i = j, extent
!         write(item%logfile,*) "what ", str1%itype(place1), itype(i)
         if (str1%itype(place1)==itype(i)) then
            pp2 = item%rank3 - i + 1
!          write(item%logfile,*) "pp2 ", pp2
!          write(item%logfile,*) "what2 ", str2%itype(place2), itype(pp2)
            if (str2%itype(place2)==
     &                             itype(pp2)) then
!               write(item%logfile,*)"correct pair ",str1%str(place1),
!     &                              " ", str2%str(place2)
               correct_pair = .true.
               ! Remove pair from itype copy
               itype(i) = 0
               itype(pp2) = 0
               !write(10,*) "itype after: ", (itype(k),k=1,INDEX_LEN)
               exit
            end if
         end if
      end do


      return
      end

*----------------------------------------------------------------------*
      subroutine check_pairing(correct_pair, str1, str2, rank1, rank2,
     &                         place1, place2, item, itype)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      logical, intent(inout) ::
     &   correct_pair
      type(itf_contr), intent(in) ::
     &   item
      type(index_str), intent(in) ::
     &   str1, str2
      integer, intent(in) ::
     &   place1, place2,
     &   rank1, rank2
      integer, intent(inout) ::
     &   itype(INDEX_LEN)

      integer ::
     &   i, j, extent, k,
     &   pp2

      ! Remeber - the itype info refers to the result positions

      !write(item%logfile,*) "is it a correct pairing?"
      !write(item%logfile,*) "place1 ", place1
      !write(item%logfile,*) "place2 ", place2


      if (item%rank3==2) then
         ! For rank 2, sometimes creation and annhilation ops are
         ! switched, so need to check both pairing combinations
         !write(10,*) "str1%itype ", str1%itype(place1)
         !write(10,*) "str2%itype ", str2%itype(place2)
         !write(10,*) "str1%itype: ", (str1%itype(k),k=1,rank1)
         !write(10,*) "str2%itype: ", (str2%itype(k),k=1,rank2)
         !write(10,*) "itype ", (itype(j),j=1,INDEX_LEN)
         if (str1%itype(place1)==item%itype(1)) then
            if (str2%itype(place2)==
     &                             item%itype(2)) then
               !write(11,*) "correct pair"
               correct_pair = .true.
            end if
         else if (str1%itype(place1)==item%itype(2)) then
            if (str2%itype(place2)==
     &                             item%itype(1)) then
               !write(11,*) "correct pair"
               correct_pair = .true.
            end if
         end if
      else
         if (place1>rank1/2) then
            j = item%rank3/2+1
            extent = item%rank3
         else
            j = 1
            extent = item%rank3/2
         end if

         !write(10,*) "str1 ", str1%str
         !write(10,*) "str2 ", str2%str
         !write(10,*) "str1%itype ", str1%itype(place1)
         !write(10,*) "str2%itype ", str2%itype(place2)
         !write(10,*) "itype: ", (itype(k),k=1,INDEX_LEN)
         !write(10,*) "str1%itype: ", (str1%itype(k),k=1,rank1)
         !write(10,*) "str2%itype: ", (str2%itype(k),k=1,rank2)
         do i = j, extent
            !write(10,*) "what ", str1%itype(place1), itype(i)
            if (str1%itype(place1)==itype(i)) then
               pp2 = item%rank3 - i + 1
               !write(10,*) "pp2 ", pp2
               !write(10,*) "what2 ", str2%itype(place2), itype(pp2)
               if (str2%itype(place2)==
     &                                itype(pp2)) then
!                  write(item%logfile,*)"correct pair ",str1%str(place1),
!     &                                 " ", str2%str(place2)
                  correct_pair = .true.
                  ! Remove pair from itype copy
                  itype(i) = 0
                  itype(pp2) = 0
                  !write(10,*) "itype after: ", (itype(k),k=1,INDEX_LEN)
                  exit
               end if
            end if
         end do
      end if


      return
      end


*----------------------------------------------------------------------*
      subroutine permute_index(idx, rank)
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
     &   rank

      integer ::
     &   i, j,
     &   sum1, sum2
      character (len=1) ::
     &   tmp

      if (rank == 6) then
         return
      end if

      if (rank == 2) then
         if (idx%itype(1)>idx%itype(2)) then
            tmp = idx%str(1)
            idx%str(1) = idx%str(2)
            idx%str(2) = tmp
         end if
         return
      end if

      sum1 = 0
      sum2 = 0
      do i = 1, rank/2
         sum1 = sum1 + idx%itype(i)
         sum2 = sum2 + idx%itype(i+rank/2)
      end do

      if (sum1<sum2) then
         do i = 1, rank/2-1
            if (idx%itype(i)>=idx%itype(i+1)) then
               if (idx%str(i)>idx%str(i+1)) then
                  tmp = idx%str(i)
                  idx%str(i) = idx%str(i+1)
                  idx%str(i+1) = tmp

                  tmp = idx%str(rank/2+i)
                  idx%str(rank/2+i) = idx%str(rank/2+i+1)
                  idx%str(rank/2+i+1) = tmp

                  ! Update possition of contraction indicies
                  do j = 1, size(idx%cnt_poss)
                     if (idx%cnt_poss(j)==i) then
                        idx%cnt_poss(j)=i+1
                     else if (idx%cnt_poss(j)==i+1) then
                        idx%cnt_poss(j)=i
                     end if
                  end do

                  !itmp = idx%itype(i)
                  !idx%itype(i) = idx%itype(i+1)
                  !idx%itype(i+1) = itmp

                  !itmp = idx%itype(rank/2+i)
                  !idx%itype(rank/2+i) = idx%itype(rank/2+i+1)
                  !idx%itype(rank/2+i+1) = itmp
               end if
            end if
         end do
      else
         do i = 1, rank/2-1
            if (idx%itype(i+rank/2)<=idx%itype(i+rank/2+1)) then
               if (idx%str(i+rank/2)<idx%str(i+rank/2+1)) then
                  tmp = idx%str(i+rank/2)
                  idx%str(i+rank/2) = idx%str(i+rank/2+1)
                  idx%str(i+rank/2+1) = tmp

                  tmp = idx%str(i)
                  idx%str(i) = idx%str(i+1)
                  idx%str(i+1) = tmp

                  do j = 1, size(idx%cnt_poss)
                     if (idx%cnt_poss(j)==i+rank/2) then
                        idx%cnt_poss(j)=i+rank/2+1
                     else if (idx%cnt_poss(j)==i+1) then
                        idx%cnt_poss(j)=i+rank/2
                     end if
                  end do

                  !itmp = idx%itype(i+rank/2)
                  !idx%itype(i+rank/2) = idx%itype(i+rank/2+1)
                  !idx%itype(i+rank/2+1) = itmp

                  !itmp = idx%itype(i)
                  !idx%itype(i) = idx%itype(i+1)
                  !idx%itype(i+1) = itmp
               end if
            end if
         end do
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine permute_slot_order(idx, rank, itype)
*----------------------------------------------------------------------*
!     Swap around pairs of indices, so as to match a given itype
!     sequence. This changes which indicies will ultimiatley end up in
!     each ITF tensor slot.
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(inout) ::
     &   idx         ! Index and itype of tensor
      integer, intent(in) ::
     &   rank,       ! Rank of tensor
     &   itype(rank)    ! Index type to compare to

      integer ::
     &   i, j, k,
     &   pp1, pp2    ! Paired position
      character (len=1), pointer ::
     &   new_str(:) => null()    ! New index string with correct slot positions
      character (len=1) ::
     &   tmp
      logical ::
     &   incorrect,  ! True if the indicies are not in the correct slots
     &   used        ! True if the pair has already been used


      ! Check if the indicies are already in the correct slots
      incorrect = .false.
      do i = 1, rank
         if (idx%itype(i) /= itype(i)) incorrect = .true.
      end do
      if (.not. incorrect) return

      ! Always have itypes of lower value on the left
      if (rank == 2) then
         if (idx%itype(1)>idx%itype(2)) then
            tmp = idx%str(1)
            idx%str(1) = idx%str(2)
            idx%str(2) = tmp
         end if
         return
      end if

      ! For higher rank tensors, search for index pair which matches the
      ! itype of a certain position (enumerated by i) and place that
      ! pair it that position in new_str
      allocate(new_str(rank))
      new_str = ''

      do i = 1, rank/2
         do j = 1, rank/2
            pp1 = rank - i + 1
            pp2 = rank - j + 1

            used = .false.
            do k = 1, rank/2
               if (new_str(k) == idx%str(j)) then
                  used = .true.
                  exit
               end if
            end do

            if (used) then
               cycle
            end if

            if (idx%itype(j)==itype(i) .and.
     &                                  idx%itype(pp2)==itype(pp1)) then
               !write(10,*) "idx%itype ", idx%itype(j)
               !write(10,*) "idx%itype ", idx%itype(pp2)
               !write(10,*) "itype ", itype(i)
               !write(10,*) "itype ", itype(pp1)
               !write(10,*) "idx%str(j) ", idx%str(j)
               !write(10,*) "idx%str(pp2) ", idx%str(pp2)
               new_str(i) = idx%str(j)
               new_str(pp1) = idx%str(pp2)
               !write(10,*) "new_str ", new_str
               exit

            else
               continue
            end if

         end do
      end do

      ! Update old result string and itype
      idx%str = new_str
      idx%itype = itype
      deallocate(new_str)

      return
      end


!*----------------------------------------------------------------------*
!      subroutine assign_index(contr_info,item)
!*----------------------------------------------------------------------*
!!     Assign an ITF index string to each tensor in a line
!*----------------------------------------------------------------------*
!
!      use itf_utils
!      implicit none
!      include 'opdim.h'
!      include 'def_contraction.h'
!      include 'def_itf_contr.h'
!
!      type(binary_contr), intent(in) ::
!     &   contr_info     ! Information about binary contraction
!      type(itf_contr), intent(inout) ::
!     &   item           ! ITF binary contraction
!
!      integer ::
!     &   c(4,2),        ! Operator numbers of contraction index
!     &   e1(4,2),       ! Operator numbers of external index 1
!     &   e2(4,2),       ! Operator numbers of external index 2
!     &   e3(4,2),       ! Operator numbers of result index
!     &   e4(4,2),       ! Operator numbers of result index
!     &   e5(4,2),       ! Operator numbers of result index
!     &   te1(4,2),      ! Test Operator numbers of external index 1
!     &   te2(4,2),      ! Test Operator numbers of external index 2
!     &   tc(4,2),       ! Test Operator numbers of external index 2
!     &   i, j, k, l,    ! Loop index
!     &   ii,            ! Letter index
!     &   nloop,         ! Number of contraction loops
!     &   i1, i2,        ! Search creation or annihilation first
!     &   shift,         ! List shift
!     &   sp             ! Pair list shift
!!      character, dimension(21) ::
!!     &   ind=(/ 'a','b','c','d','e','f','g','i','j','k','l','m','n',
!!     &          'o','p','q','r','s','t','u','v' /)   ! Letters for index string
!      character, dimension(21) ::
!     &   ind=(/ 'a','b','c','d','e','f','g','p','q','r','s','t','u',
!     &          'v','i','j','k','l','m','n','o' /)   ! Letters for index string
!      character(len=INDEX_LEN) ::
!     &     s1, s2, s3   ! Tmp ITF index strings
!      character(len=1) ::
!     &   tmp
!      real(8) ::
!     &   factor         ! Factor from equivalent lines
!      logical ::
!     &   found,         ! True if found pairing index
!     &   sort           ! Used in bubble sort
!      type(pair_list) ::
!     &   p_list,        ! Complete list of pairs in binary contraction
!     &   t1_list,       ! List of pairs in first tensor (T1)
!     &   t2_list,       ! List of pairs in second tensor (T2)
!     &   r_list,        ! List of pairs in result tensor
!     &   d_list         ! debug list
!      integer, dimension(4) ::
!     &   t_shift,       ! Index shift for external indices
!     &   c_shift        ! Index shift for contraction indices
!      integer, dimension(2) ::
!     &   cops,          ! Number of creation/annihilation contraction operators
!     &   e1ops,         ! Number of C/A external T1 operators
!     &   e2ops,         ! Number of C/A external T2 operators
!     &   e3ops          ! Total number of external operators
!      integer, dimension(4) ::
!     &   tops           ! Scratch space
!
!      c=0
!      e1=0
!      e2=0
!      e3=0
!
!      e4=0
!      e5=0
!
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
!      do i = 1, contr_info%nj_res
!        call count_index(i,
!     &     contr_info%occ_res(1:,1:,i),
!     &     contr_info%rst_res(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,e3)
!      end do
!
!      do i = 1, contr_info%nj_op1
!        call count_index(i,
!     &     contr_info%occ_op1(1:,1:,i),
!     &     contr_info%rst_op1(1:,1:,1:,1:,1:,i),
!     &     contr_info%ngas,contr_info%nspin,e4)
!      end do
!      e5 = e3
!
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
!      ! Find out number of creation/annihilation operators per operator
!      cops = sum(c, dim=1)
!      e1ops = sum(e1, dim=1)
!      e2ops = sum(e2, dim=1)
!      e3ops = sum(e3, dim=1)
!
!      ! Set ranks of tensors
!      call itf_rank(e1, c, item%rank1, .false.)
!      call itf_rank(e2, c, item%rank2, .false.)
!      call itf_rank(e3, c, item%rank3, .true.)
!
!      ! Set number of indcies
!      item%nops1 = sum(e1, dim=2) + sum(c, dim=2)
!      item%nops2 = sum(e2, dim=2) + sum(c, dim=2)
!      item%nops3 = sum(e3, dim=2)
!
!      ! Set if it has three operators of the same type
!      ! TODO: probably not needed anymore
!      tops = sum(e1, dim=2)
!      do i = 1, ngastp
!         if (tops(i)==3) then
!            item%three(1) = .true.
!         end if
!      end do
!
!      tops = sum(e1, dim=2)
!      do i = 1, ngastp
!         if (tops(i)==3) then
!            item%three(2) = .true.
!         end if
!      end do
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
!!      ! New code ======================================================
!!      allocate(d_list%plist(item%rank2/2))
!!      sp=1
!!      do i = 1, 4
!!         c_shift(i) = e1(i,1) + e1(i,2) + e2(i,1) + e2(i,2)
!!      end do
!!      t_shift=0
!!      te1 = e1
!!      te2 = e2
!!      tc = c
!!      e1ops = sum(te1, dim=1)
!!      e2ops = sum(te2, dim=1)
!!
!!      !write(item%logfile,*) "e1 ", te1
!!      !write(item%logfile,*) "e2 ", te2
!!      !write(item%logfile,*) "c ", tc
!!      ! Assume exciation operator always comes second
!!      if (trim(item%label_t2)=='T2g' .and. item%rank2>2) then
!!
!!         do while (sum(e2ops) /= 0)
!!            ! Search on second tensor
!!            ! Note the number of creation/annihilation ops has been
!!            ! switched
!!            call find_T_pairs(d_list,sp,2, te2, te1,tc,t_shift, c_shift,
!!     &                      e2ops, e1ops, item)
!!
!!            ! Check for more operators
!!            e2ops = sum(te2, dim=1)
!!         end do
!!
!!         !call print_plist(d_list, item%rank2/2, "d_list", item%logfile)
!!      end if
!!
!!
!!      e1ops = sum(e1, dim=1)
!!      e2ops = sum(e2, dim=1)
!!      deallocate(d_list%plist)
!!      ! End new code ==================================================
!
!      ! To start, we pair off contraction loops
!      ! Pair list index (shift pair)
!      sp = 1
!
!      ! Keep searching until all paired loops are found (ie. the number
!      ! of creation or annihilation operators becomes 0)
!
!      ! Set letter shift values for contraction indices
!      do i = 1, 4
!         c_shift(i) = e1(i,1) + e1(i,2) + e2(i,1) + e2(i,2)
!      end do
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
!
!      return
!      end


!*----------------------------------------------------------------------*
!      subroutine find_pairs(list, sp, tensor, e1, e2, c, t_shift,
!     &                      c_shift, e1ops, e2ops, item)
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
!!      character, dimension(21) ::
!!     &   ind=(/ 'a','b','c','d','e','f','g','i','j','k','l','m','n',
!!     &          'o','p','q','r','s','t','u','v' /)   ! Letters for index string
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
!     &   i,j,k,l,m,n,      ! Loop indices
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
!            ii = 1+(8*(i-1)) + t_shift(i)
!            list%plist(sp)%pindex(i1) = ind(ii)
!            call assign_nval(i, i1, sp, list)
!            list%plist(sp)%ops(i1) = tensor
!
!            t_shift(i) = t_shift(i) + 1
!            e1(i,i1) = e1(i,i1) - 1
!
!            ! Look for matching operator
!            if (n1 > 0) then
!               ! Look on the same tensor
!               do k = 1, 4
!                  do l = 1, e1(k,i2)
!                     ii = 1+(8*(k-1)) + t_shift(k)
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
!            else if (n2 > 0) then
!               ! Look on the opposite tensor
!               do k = 1, 4
!                  do l = 1, e2(k,i2)
!                     ii = 1+(8*(k-1))+t_shift(k)
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
!                     do m=1, 4
!                        do n=1, c(m,i2)
!                            ii = 1+(8*(m-1)) + c_shift(m)
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
!            else
!               ! No ops of opposite type to match...
!               call line_error("Particle number not conserving",item)
!            end if
!            ! A pair loop has been found so exit
!            if (found) exit
!         end do
!         if (found) exit
!      end do
!
!      return
!      end


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
      !if (item%rank1 + item%rank2 + item%rank3 == 4) then
      if (item%rank1 == 2 .and. item%rank2 == 0 .or.
     &    item%rank1 == 0 .and. item%rank2 == 2) then
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

         call print_itf_line(item,.false.,.false.)
         return
      else if (item%rank1 + item%rank2 + item%rank3 == 6) then
         ! This covers all cases where we have three rank-2 tensors

         if (item%inter(1)) then
            item%label_t1 = trim(item%label_t1)//'aa'
         end if
         if (item%inter(2)) then
            item%label_t2 = trim(item%label_t2)//'aa'
         end if

         call print_itf_line(item,.false.,.false.)
         return
      else if (item%rank1 + item%rank2 + item%rank3 == 0) then
         ! Scalar contributions
         call print_itf_line(item,.false.,.false.)
         return
      else if (item%rank1 == 0 .or. item%rank2 ==0) then
         ! Tensor multiplied by a scalar (not involving an intermediate)

         if (item%inter(1)) then
            call simple_spin_name(item%inter1,item%rank1,item%permute)
         else if (item%inter(2)) then
            call simple_spin_name(item%inter2,item%rank2,item%permute)
         end if

         call print_itf_line(item,.false.,.false.)
         return
      else if (item%rank3==4 .and. item%rank1==2
     &         .and. item%rank2==2) then

         ! Tensor product
         if (item%inter(1)) then
            item%label_t1 = trim(item%label_t1)//'aa'
         else if (item%inter(2)) then
            item%label_t2 = trim(item%label_t2)//'aa'
         end if

         call print_itf_line(item,.false.,.false.)
         return
      end if

      !write(item%logfile,*) "inside assign_spin"

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
               do i=1, 2
                  ! aaaaaa
                  item%t_spin(3)%spin(1,i) = 1
                  item%t_spin(3)%spin(2,i) = 1
               end do
               item%t_spin(3)%spin(1,3) = 1
               item%t_spin(3)%spin(2,3) = 1
            case default
               call line_error("Could not determine tensor rank",item)
         end select
      end if

!      call print_spin(item%t_spin(3)%spin,item%rank3,"Result",
!     &                item%logfile)

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

!      call print_spin(item%t_spin(1)%spin, item%rank1, "T1 before",
!     &                item%logfile)
!      call print_spin(item%t_spin(2)%spin, item%rank2, "T2 before",
!     &                item%logfile)

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
      character(len=1) ::
     &   p1, p2
      logical ::
     &   eloop,
     &   error

      allocate(poss(2,item%contri))

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

      !call print_spin(item%t_spin(1)%spin,item%rank1,"T1",item%logfile)
      !call print_spin(item%t_spin(2)%spin,item%rank2,"T2",item%logfile)

      ! Get position information of contraction indicies, ie. look for
      ! zeros
      ! TODO: know this information already when assigning index?
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
                  ! For scalar results, only need half of the spin
                  ! cases, the rest are the same
                  if (item%rank3 == 0 .and. i == 2) exit
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
                              if (item%rank3 == 0 .and. i == 2) exit
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
                                  if (item%rank3 == 0 .and. i == 2) exit
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
         ! Check other reasons why the spin case wasn't printed
         error = .false.
         ! Check all indicies are assigned
         do i = 1, item%rank3
            if(item%idx3(i:i)==' ') then
               error = .true.
               exit
            end if
         end do

         ! Check correctly paired
         ! Warning: this doesn't check pairs across the tensors!
         if (.not. error) then
            do i = 1, item%rank3/2
               p1 = item%idx3(i:i)
               p2 = item%idx3(i+item%rank3/2:i+item%rank3/2)

               do j = 1, item%rank1/2
                  if (p1 == item%idx1(j:j)) then
                     if (p2 /=
     &                   item%idx1(j+item%rank1/2:j+item%rank1/2)) then
                        do k = 1, item%rank2/2
                           if (p2 ==
     &                    item%idx2(k+item%rank2/2:k+item%rank2/2)) then
                              if
     &                        (item%idx1(j+item%rank1/2:j+item%rank1/2)
     &                        /=
     &                        item%idx2(k:k))
     &                        then
                                error = .true.
                                exit
                              end if

                              exit
                           end if

                        end do

                     end if
                  end if
               end do

               do j = 1, item%rank2/2
                  if (p1 == item%idx2(j:j)) then
                     if (p2 /=
     &                   item%idx2(j+item%rank2/2:j+item%rank2/2)) then
                        do k = 1, item%rank1/2
                           if (p2 ==
     &                    item%idx1(k+item%rank1/2:k+item%rank1/2)) then
                              if
     &                        (item%idx2(j+item%rank2/2:j+item%rank2/2)
     &                        /=
     &                        item%idx1(k:k))
     &                        then
                                error = .true.
                                exit
                              end if

                              exit
                           end if

                        end do

                     end if
                  end if
               end do

            end do
         end if

         if (error) then
           call line_error("Didn't print out spin case", item)
         else
           call line_error("This spin case possibly doesn't exist",item)
         end if
      end if

      deallocate(poss)

      return
      end


*----------------------------------------------------------------------*
      subroutine assign_spin2(item)
*----------------------------------------------------------------------*
!     Assign spin to tensors, then sum remaining contraction indices
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
     &   item

      integer ::
     &   i,j,k,l,m,n,             ! Loop index
     &   result_spin(4),   ! Spin case of result
     &   hr1, hr2, hr3     ! Half tensor ranks

      ! Calculate half ranks for use in indexing letter index strings
      hr1 = item%rank1/2
      hr2 = item%rank2/2
      hr3 = item%rank3/2

      ! Assign spin to indicies on the result tensor
      if (.not. item%inter(3)) then
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
               do i=1, 2
                  ! aaaaaa
                  item%t_spin(3)%spin(1,i) = 1
                  item%t_spin(3)%spin(2,i) = 1
               end do
               item%t_spin(3)%spin(1,3) = 1
               item%t_spin(3)%spin(2,3) = 1
            case default
               call line_error("Could not determine tensor rank",item)
         end select

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
         ! Sum over the remaining contraction indicies and print out the
         ! line
!      call print_spin(item%t_spin(3)%spin,item%rank3,item%label_res,
!     &                item%logfile)
!      call print_spin(item%t_spin(1)%spin,item%rank3,"T1",
!     &                item%logfile)
!      call print_spin(item%t_spin(2)%spin,item%rank3,"T2",
!     &                item%logfile)
         call spin_index2(item)

      else

      if (item%rank3>0) then
         do n = 1, 2
            item%t_spin(3)%spin(1,1) = n
            do m = 1, 2
               item%t_spin(3)%spin(2,1) = m
               if (item%rank3==2) then
      ! Assign spin of external indicies to T1 and T2
      item%t_spin(1)%spin = 0
      item%t_spin(2)%spin = 0
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

!      call print_spin(item%t_spin(3)%spin,item%rank3,"Result",
!     &                item%logfile)
      call spin_index2(item)

               else
                  do k = 1, 2
                     item%t_spin(3)%spin(1,2) = k
                     do l = 1, 2
                        item%t_spin(3)%spin(2,2) = l
      ! Assign spin of external indicies to T1 and T2
      item%t_spin(1)%spin = 0
      item%t_spin(2)%spin = 0
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

!      call print_spin(item%t_spin(3)%spin,item%rank3,"Result2",
!     &                item%logfile)
      call spin_index2(item)
                     end do
                  end do
               end if
            end do
         end do

      end if

      end if

!      call print_spin(item%t_spin(3)%spin,item%rank3,"Result",
!     &                item%logfile)

!      ! Assign spin of external indicies to T1 and T2
!      do j=1, hr3
!         do i=1, hr1
!            ! Assign spin of first tensor
!            if (item%idx3(j:j)==item%idx1(i:i)) then
!               item%t_spin(1)%spin(1,i) = item%t_spin(3)%spin(1,j)
!            else if (item%idx3(j:j)==item%idx1(i+hr1:i+hr1)) then
!               item%t_spin(1)%spin(2,i) = item%t_spin(3)%spin(1,j)
!            end if
!
!            if (item%idx3(j+hr3:j+hr3)==item%idx1(i:i)) then
!               item%t_spin(1)%spin(1,i) = item%t_spin(3)%spin(2,j)
!            else if (item%idx3(j+hr3:j+hr3)==item%idx1(i+hr1:i+hr1))then
!               item%t_spin(1)%spin(2,i) = item%t_spin(3)%spin(2,j)
!            end if
!         end do
!
!         do i=1, hr2
!            ! Assign spin of second tensor
!            if (item%idx3(j:j)==item%idx2(i:i)) then
!               item%t_spin(2)%spin(1,i) = item%t_spin(3)%spin(1,j)
!            else if (item%idx3(j:j)==item%idx2(i+hr2:i+hr2)) then
!               item%t_spin(2)%spin(2,i) = item%t_spin(3)%spin(1,j)
!            end if
!
!            if (item%idx3(j+hr3:j+hr3)==item%idx2(i:i)) then
!               item%t_spin(2)%spin(1,i) = item%t_spin(3)%spin(2,j)
!            else if (item%idx3(j+hr3:j+hr3)==item%idx2(i+hr2:i+hr2))then
!               item%t_spin(2)%spin(2,i) = item%t_spin(3)%spin(2,j)
!            end if
!         end do
!      end do

!      call print_spin(item%t_spin(1)%spin, item%rank1, "T1 before",
!     &                item%logfile)
!      call print_spin(item%t_spin(2)%spin, item%rank2, "T2 before",
!     &                item%logfile)

      ! Sum over the remaining contraction indicies and print out the
      ! line
!      call spin_index2(item)

      return
      end


*----------------------------------------------------------------------*
      subroutine spin_index2(item)
*----------------------------------------------------------------------*
!     Find contraction index used in spin summation
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
     &     item

      type(twodarray), pointer ::
     &   poss(:,:) => null()
      integer ::
     &   i, j, k, l, m, n, shift,
     &   i1, i2, i3, i4,
     &   z1, z2, r1, r2
      character(len=INDEX_LEN) ::
     &   str1, str2
      character(len=1) ::
     &   p1, p2
      logical ::
     &   eloop,
     &   error

      if (item%contri>0)  then
         allocate(poss(2,item%contri))

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

      !call print_spin(item%t_spin(1)%spin,item%rank1,"T1",item%logfile)
      !call print_spin(item%t_spin(2)%spin,item%rank2,"T2",item%logfile)

      ! Get position information of contraction indicies, ie. look for
      ! zeros
      ! TODO: know this information already when assigning index?
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

      if (shift == 1 .and. item%contri>0) then
         call line_error("Didn't find contraction index", item)
      end if

!      if (item%contri>0) then
!      write(item%logfile,*) "spin z1: ", size(item%t_spin(z1)%spin,2)
!      write(item%logfile,*) "spin z2: ", size(item%t_spin(z2)%spin,2)
!      write(item%logfile,*) "CONTRI: ", item%contri
!      write(item%logfile,*) "SHIFT: ", shift
!      write(item%logfile,*) "========================"
!      do i = 1, 2
!         do j = 1, item%contri
!            write(item%logfile,*) "DEBUG1 poss: ", poss(i,j)%elements
!         end do
!      end do
!      write(item%logfile,*) "========================"
!      end if

      end if

      ! Main spin summation loop
      shift = shift - 1
      eloop = .false.

      if (item%contri/=0) then
      do i = 1, 2
         i1 = poss(1,1)%elements(1)
         i2 = poss(1,1)%elements(2)
         i3 = poss(2,1)%elements(1)
         i4 = poss(2,1)%elements(2)
         item%t_spin(z1)%spin(i1, i2) = i
         item%t_spin(z2)%spin(i3, i4) = i
         if (shift <= 1) then
            call print_spin_case2(item,eloop)
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
                  ! For scalar results, only need half of the spin
                  ! cases, the rest are the same
                  if (item%rank3 == 0 .and. i == 2) exit
                  call print_spin_case2(item,eloop)
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
                        call print_spin_case2(item,eloop)
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
                              if (item%rank3 == 0 .and. i == 2) exit
                              call print_spin_case2(item,eloop)
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
                                    call print_spin_case2(item,eloop)
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
                                  if (item%rank3 == 0 .and. i == 2) exit
                                       call print_spin_case2(item,eloop)
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

      else
         ! No contraction indicies - tensors have same spin as result
         if (r1>0) item%t_spin(z1)%spin = item%t_spin(3)%spin
         if (r2>0) item%t_spin(z2)%spin = item%t_spin(3)%spin

         call print_spin_case2(item,eloop)
      end if

!      if (.not. eloop) then
!         ! Check other reasons why the spin case wasn't printed
!         error = .false.
!         ! Check all indicies are assigned
!         do i = 1, item%rank3
!            if(item%idx3(i:i)==' ') then
!               error = .true.
!               exit
!            end if
!         end do
!
!         ! Check correctly paired
!         ! Warning: this doesn't check pairs across the tensors!
!         if (.not. error) then
!            do i = 1, item%rank3/2
!               p1 = item%idx3(i:i)
!               p2 = item%idx3(i+item%rank3/2:i+item%rank3/2)
!
!               do j = 1, item%rank1/2
!                  if (p1 == item%idx1(j:j)) then
!                     if (p2 /=
!     &                   item%idx1(j+item%rank1/2:j+item%rank1/2)) then
!                        do k = 1, item%rank2/2
!                           if (p2 ==
!     &                    item%idx2(k+item%rank2/2:k+item%rank2/2)) then
!                              if
!     &                        (item%idx1(j+item%rank1/2:j+item%rank1/2)
!     &                        /=
!     &                        item%idx2(k:k))
!     &                        then
!                                error = .true.
!                                exit
!                              end if
!
!                              exit
!                           end if
!
!                        end do
!
!                     end if
!                  end if
!               end do
!
!               do j = 1, item%rank2/2
!                  if (p1 == item%idx2(j:j)) then
!                     if (p2 /=
!     &                   item%idx2(j+item%rank2/2:j+item%rank2/2)) then
!                        do k = 1, item%rank1/2
!                           if (p2 ==
!     &                    item%idx1(k+item%rank1/2:k+item%rank1/2)) then
!                              if
!     &                        (item%idx2(j+item%rank2/2:j+item%rank2/2)
!     &                        /=
!     &                        item%idx1(k:k))
!     &                        then
!                                error = .true.
!                                exit
!                              end if
!
!                              exit
!                           end if
!
!                        end do
!
!                     end if
!                  end if
!               end do
!
!            end do
!         end if
!
!         if (error) then
!           call line_error("Didn't print out spin case", item)
!         else
!           call line_error("This spin case possibly doesn't exist",item)
!         end if
!      end if

      if (item%contri>0) deallocate(poss)

      return
      end


*----------------------------------------------------------------------*
      subroutine print_spin_case2(item,eloop)
*----------------------------------------------------------------------*
!     Print spin case
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
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

!            call print_spin(item%t_spin(3)%spin, item%rank3, "Spin 3",
!     &                      item%logfile)
!            call print_spin(item%t_spin(1)%spin, item%rank1, "Spin 1",
!     &                      item%logfile)
!            call print_spin(item%t_spin(2)%spin, item%rank2, "Spin 2",
!     &                      item%logfile)

            item%all_spins(item%nspin_cases)%t_spin(1)%spin=
     &                                          item%t_spin(1)%spin
            item%all_spins(item%nspin_cases)%t_spin(2)%spin=
     &                                          item%t_spin(2)%spin
            item%all_spins(item%nspin_cases)%t_spin(3)%spin=
     &                                          item%t_spin(3)%spin
            item%nspin_cases = item%nspin_cases + 1



            ! TODO: need eloop??
            eloop=.true.
         end if
      end if


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

!            call print_spin(item%t_spin(1)%spin, item%rank1, "Spin 1",
!     &                      item%logfile)
!            call print_spin(item%t_spin(2)%spin, item%rank2, "Spin 2",
!     &                      item%logfile)

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

               item%ninter = ishift
            end if

            ! Print the spin summed line
            if (item%print_line) then
               ! Mark the start of the spin summed block, if the current
               ! line is a permutation of the previous line (ie. permute >1)
               ! then we should not begin a new block
               call print_itf_line(item,s1,s2)
            end if

            eloop=.true.
         end if
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine inter_spin_name2(spin,hrank,label)
*----------------------------------------------------------------------*
!     Add spin name to intermediate (ie. STIN001 -> STIN001abab)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   hrank,            ! Rank/2
     &   spin(2,INDEX_LEN/2)    ! Spin info
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
      subroutine simple_spin_name(spin_name,rank,permute)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=INDEX_LEN), intent(inout) ::
     &   spin_name             ! Spin name of intermediate
      integer, intent(in) ::
     &   rank,
     &   permute

      select case (rank)
         case (0)
            return
         case (2)
            spin_name = 'aa'
         case (4)
            if (permute==2) then
               spin_name = 'baba'
            else
               spin_name = 'abab'
            end if
         case (6)
            spin_name = 'aaaaaa'
      end select

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
     &   i, j,
     &   shift

      logical ::
     &   co_v,
     &   contr_v

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

      ! Should check here if repeated spin cases are being added
      ! If the algo requests a already delcared spin case, then just
      ! skip it - we only need it once
      if (item%inter_spins(ishift)%ncase>=1) then
      do i = 1, item%inter_spins(ishift)%ncase
         co_v = .false.
         contr_v = .false.
         do j = 1, hrank
            if (item%inter_spins(ishift)%cases(j,i)/=spin(1,j)) then
               co_v = .true.
               exit
            end if
            if (item%inter_spins(ishift)%cases(j+INDEX_LEN/2,i)
     &          /=spin(2,j)) then
               contr_v = .true.
               exit
            end if
         end do
         if (.not. co_v .and. .not. contr_v) then
            return
         end if
      end do

      end if

      do i=1, hrank
         item%inter_spins(ishift)%cases(i,shift)=spin(1,i)
         item%inter_spins(ishift)%cases(i+INDEX_LEN/2,shift)=spin(2,i)
      end do

      ! Update number of spin cases
      item%inter_spins(ishift)%ncase =
     &                           item%inter_spins(ishift)%ncase + 1

      return
      end


*----------------------------------------------------------------------*
      subroutine itf_contr_init(contr_info,item,perm,itin,comm,lulog)
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
      logical, intent(in) ::
     &     itin

      integer :: i, j, ct(ngastp,2)

      ! Assign output file
      item%logfile=lulog

      ! Assign command type
      item%command=comm

      ! Get number of contraction and external indices on each tensor
      call itf_ops(contr_info, item%c, item%e1, item%e2, item%command)

      ! Set ranks of tensors using matricies from itf_ops
      call itf_rank(item%e1, item%c, item%rank1, item%nops1, .false.)
      call itf_rank(item%e2, item%c, item%rank2, item%nops2, .false.)
      call itf_rank(item%e1, item%e2, item%rank3, item%nops3, .false.)

      ! Set number of contraction indicies
      item%contri = sum(sum(item%c, dim=1))

      ! Set external lines (used for permuations)
      item%perm_case = contr_info%perm

      ! Determine factor from equivalent lines
      item%fact = 1.0d+0
      call itf_equiv_lines_factor(item%c, item%fact)

      ! Get any remaining factors from GeCCo
      item%fact = item%fact * abs(contr_info%fact_itf)
      !item%fact = item%fact * contr_info%fact_itf


      ! Assign labels
      item%label_t1=contr_info%label_op1

      ! Check if an intermediate
      item%inter(1) = check_inter(item%label_t1)
      if (.not. item%inter(1)) then
         ! Check if a density matrix
         item%den(1) = check_den(item%label_t1)

         if (.not. item%den(1)) then
            ! Check if an integral
            item%int(1) = check_int(item%label_t1)
            if (item%int(1) .and. item%rank1==4) then
               call check_j_integral(item%e1, item%c, item%nops1,
     &                               item%j_int, .false.)
            end if
         end if
      end if


      if (item%command/=command_cp_intm .or.
     &    item%command/=command_add_intm) then

         ! Operator 2 does not exist in [ADD] or [COPY]
         item%label_t2=contr_info%label_op2

         ! Assign permutation number
         item%permute=perm
         if (perm>1) then
            item%permutation = .true.
         end if

         item%inter(2) = check_inter(item%label_t2)

         if (.not. item%inter(2)) then
            item%den(2) = check_den(item%label_t2)

            if (.not. item%den(2)) then
               item%int(2) = check_int(item%label_t2)
               if (item%int(2) .and. item%rank2==4) then
                  call check_j_integral(item%e2, item%c, item%nops2,
     &                                  item%j_int, .true.)
               end if
            end if
         end if
      end if

      item%label_res=contr_info%label_res
      item%inter(3) = check_inter(item%label_res)
      item%int(3) = .false.


      ! If a residual, is it symmetric (R_{ab}^{ij] = R_{ba}^{ji})?
      if (.not. item%inter(3)) then
         call check_symmetric(contr_info, item%command, item%symm_res)

         ! Add factor to account for factor of two when the residual is
         ! symmetrised:
         ! R:eecc[abij] += G:eecc[abij]
         ! R:eecc[abij] += G:eecc[baji]
         if (item%symm_res .and. item%permute==0) then
            if (contr_info%label_res/='INTpp') then
               item%fact = item%fact * 0.5d+0
            end if
         end if
      end if


      ! Remove zeorth-body density from equations
      if (item%label_t2 == 'GAM0' .and. item%rank2 == 0) then
         item%label_t2 = ''
         item%command = command_add_intm
      else if (item%label_t1 == 'GAM0' .and. item%rank1 == 0
     &         .and. item%command /= command_add_intm) then
         item%label_t1 = item%label_t2
         item%rank1 = item%rank2
         item%e1 = item%e2
         item%nops1 = item%nops2
         item%int(1) = item%int(2)
         item%inter(1) = item%inter(2)
         item%label_t2 = ''
         item%e2 = 0
         item%nops2 = 0
         item%int(2) = .false.
         item%inter(2) = .false.
         item%command = command_add_intm
      end if


      ! Assign factor --- use special ITF factor
      ! the ITF factor is closer to the value expected from standard
      ! diagram rules, but still some care has to be taken when translating
      ! to ITF tensors; e.g. the (HP;HP) integrals are stored by GeCCo as
      ! <aj||bi> while the standard sign for ring terms assumes <aj||ib>
      !item%fact=contr_info%fact_itf


      ! Assign index string. Tensor ranks and number
      ! of operators are also set here
      if (item%command==command_cp_intm .or.
     &    item%command==command_add_intm) then
         ! For [ADD] and [COPY]
         ! Not a binary contraction
         item%binary = .false.
         !call assign_add_index(contr_info,item)
         call assign_new_index(contr_info,item)
      else
         ! For other contractions
         !call assign_index(contr_info,item)
         call assign_new_index(contr_info,item)
      end if


      ! Set up arrays to store information about intermediates
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

      item%spin_cases = 0



      ! Check if a tensor product
      if (item%rank3==4 .and. item%rank1==2 .and. item%rank2==2) then
         item%product=.true.
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine itf_contr_init2(contr_info,item,perm,itin,comm,lulog,
     &                           itype)
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
      type(itf_contr2), intent(inout) ::
     &     item     ! Object which holds information necessary to print out an ITF algo line
      integer, intent(in) ::
     &     perm,        ! Permutation information
     &     comm,        ! formula_item command
     &     lulog,        ! Output file
     &   itype(INDEX_LEN)
      logical, intent(in) ::
     &     itin

      integer :: i, j, ct(ngastp,2)

      ! Assign output file
      item%logfile=lulog

      ! Assign command type
      item%command=comm

      ! Get number of contraction and external indices on each tensor
      call itf_ops(contr_info, item%c, item%e1, item%e2, item%command)

      ! Set ranks of tensors using matricies from itf_ops
      call itf_rank(item%e1, item%c, item%rank1, item%nops1, .false.)
      call itf_rank(item%e2, item%c, item%rank2, item%nops2, .false.)
      call itf_rank(item%e1, item%e2, item%rank3, item%nops3, .false.)

      ! Set number of contraction indicies
      item%contri = sum(sum(item%c, dim=1))

      ! Set external lines (used for permuations)
      item%perm_case = contr_info%perm

      ! Determine factor from equivalent lines
      item%fact = 1.0d+0
      call itf_equiv_lines_factor(item%c, item%fact)

      ! Get any remaining factors from GeCCo
      item%fact = item%fact * abs(contr_info%fact_itf)
      !item%fact = item%fact * contr_info%fact_itf


      ! Assign labels
      item%label_t1=contr_info%label_op1

      ! Check if an intermediate
      item%inter(1) = check_inter(item%label_t1)
      if (.not. item%inter(1)) then
         ! Check if a density matrix
         item%den(1) = check_den(item%label_t1)

         if (.not. item%den(1)) then
            ! Check if an integral
            item%int(1) = check_int(item%label_t1)
            if (item%int(1) .and. item%rank1==4) then
               call check_j_integral(item%e1, item%c, item%nops1,
     &                               item%j_int, .false.)
            end if
         end if
      end if


      if (item%command/=command_cp_intm .or.
     &    item%command/=command_add_intm) then

         ! Operator 2 does not exist in [ADD] or [COPY]
         item%label_t2=contr_info%label_op2

         ! Assign permutation number
         item%permute=perm
         if (perm>1) then
            item%permutation = .true.
         end if

         item%inter(2) = check_inter(item%label_t2)

         if (.not. item%inter(2)) then
            item%den(2) = check_den(item%label_t2)

            if (.not. item%den(2)) then
               item%int(2) = check_int(item%label_t2)
               if (item%int(2) .and. item%rank2==4) then
                  call check_j_integral(item%e2, item%c, item%nops2,
     &                                  item%j_int, .true.)
               end if
            end if
         end if
      end if

      item%label_res=contr_info%label_res
      item%inter(3) = check_inter(item%label_res)
      item%int(3) = .false.

      call check_symmetric(contr_info, item%command, item%symmetric)

      ! If a residual, is it symmetric (R_{ab}^{ij] = R_{ba}^{ji})?
      if (.not. item%inter(3)) then
         call check_symmetric(contr_info, item%command, item%symm_res)

         ! Add factor to account for factor of two when the residual is
         ! symmetrised:
         ! R:eecc[abij] += G:eecc[abij]
         ! R:eecc[abij] += G:eecc[baji]
         if (item%symm_res .and. item%permute==0) then
            if (contr_info%label_res/='INTpp') then
               item%fact = item%fact * 0.5d+0
            end if
         end if
      end if


      ! Remove zeorth-body density from equations
      if (item%label_t2 == 'GAM0' .and. item%rank2 == 0) then
         item%label_t2 = ''
         item%command = command_add_intm
      else if (item%label_t1 == 'GAM0' .and. item%rank1 == 0
     &         .and. item%command /= command_add_intm) then
         item%label_t1 = item%label_t2
         item%rank1 = item%rank2
         item%e1 = item%e2
         item%nops1 = item%nops2
         item%int(1) = item%int(2)
         item%inter(1) = item%inter(2)
         item%label_t2 = ''
         item%e2 = 0
         item%nops2 = 0
         item%int(2) = .false.
         item%inter(2) = .false.
         item%command = command_add_intm
      end if


      ! Assign factor --- use special ITF factor
      ! the ITF factor is closer to the value expected from standard
      ! diagram rules, but still some care has to be taken when translating
      ! to ITF tensors; e.g. the (HP;HP) integrals are stored by GeCCo as
      ! <aj||bi> while the standard sign for ring terms assumes <aj||ib>
      !item%fact=contr_info%fact_itf


      ! Assign index string. Tensor ranks and number
      ! of operators are also set here
      if (item%command==command_cp_intm .or.
     &    item%command==command_add_intm) then
         ! For [ADD] and [COPY]
         ! Not a binary contraction
         item%binary = .false.
         !call assign_add_index(contr_info,item)
         !call assign_new_index(contr_info,item)
      else
         ! For other contractions
         !call assign_index(contr_info,item)
         !call assign_new_index(contr_info,item)
      end if


      ! Set up arrays to store information about intermediates
      item%inter1 = ''
      item%inter2 = ''

      ! TODO: bit of a hack to avoid seg faults for rank == 0
      ! Instead have a check when assiging spins if == 0
      if (item%rank1>0) then
         allocate(item%t_spin(1)%spin(2, item%rank1/2))
      else
         allocate(item%t_spin(1)%spin(2, 2))
      end if

      if (item%rank2>0) then
         allocate(item%t_spin(2)%spin(2, item%rank2/2))
      else
         allocate(item%t_spin(2)%spin(2, 2))
      end if

      if (item%rank3>0) then
         allocate(item%t_spin(3)%spin(2, item%rank3/2))
      else
         allocate(item%t_spin(3)%spin(2, 2))
      end if


      item%t_spin(1)%spin = 0
      item%t_spin(2)%spin = 0
      item%t_spin(3)%spin = 0

      if (item%inter(3)) then
         allocate(item%i_spin%spin(2, item%rank3/2))
      end if

      ! Initalise the number of different spin cases
      ! Each line has a minimum of one spin case
      item%spin_cases = 1



      ! Check if a tensor product
      if (item%rank3==4 .and. item%rank1==2 .and. item%rank2==2) then
         item%product=.true.
      end if

      ! Number of spin cases associated with this line
      item%nspin_cases=1


      ! Set itype from previous line (can be 0)
      item%itype = itype

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
      subroutine itf_deinit2(item)
*----------------------------------------------------------------------*
!     Deinitialise ITF contraction object
*----------------------------------------------------------------------*

      use itf_utils
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
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
      subroutine print_itf_contr2(item)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr2), intent(inout) ::
     &     item     ! Object which holds information necessary to print out an ITF algo line

      integer ::
     &   i,j,k,
     &   l1,l2,l3

      write(item%logfile,*)
      write(item%logfile,*) "================================="
      write(item%logfile,*) trim(item%label_res)//"["//
     & trim(item%idx3)//"]"//
     & " = "//trim(item%label_t1)//"["//
     & trim(item%idx1)//"]"//" "//
     & trim(item%label_t2)//"["//
     & trim(item%idx2)//"]"
      write(item%logfile,*) "FACTOR: ", item%fact
      write(item%logfile,*) "SPIN CASES: ", item%spin_cases-1

      l1 = item%rank1/2
      l2 = item%rank2/2
      l3 = item%rank3/2

      do i = 1, item%spin_cases-1
      write(item%logfile,*)
      write(item%logfile,*) "---------------------------------"
      do k = 2, 1, -1
      write(item%logfile,'(2a)',advance='no') "  "
      do j = 1, l3
      write(item%logfile,'(i1)',advance='no')
     &     item%all_spins(i)%t_spin(3)%spin(k,j)
      end do
      write(item%logfile,'(2a)',advance='no') "  "
      do j = 1, l1
      write(item%logfile,'(i1)',advance='no')
     &     item%all_spins(i)%t_spin(1)%spin(k,j)
      end do
      write(item%logfile,'(2a)',advance='no') "  "
      do j = 1, l2
      write(item%logfile,'(i1)',advance='no')
     &     item%all_spins(i)%t_spin(2)%spin(k,j)
      end do
      write(item%logfile,*)
      end do
      end do

      write(item%logfile,*) "================================="
      write(item%logfile,*)

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
      write(item%logfile,'(a)') trim(line)

      if (item%product) then
         ! We don't need to symmetrise from a tensor product
         return
      end if

      tindex = ' '
      tindex(1:1) = item%idx3(2:2)
      tindex(2:2) = item%idx3(1:1)
      tindex(3:3) = item%idx3(4:4)
      tindex(4:4) = item%idx3(3:3)

      line = '.'//trim(new)//'['//trim(item%idx3)//'] += '//
     &       trim(item%label_res)//'['//trimal(tindex)//']'
      write(item%logfile,'(a)') trim(line)


      return
      end


*----------------------------------------------------------------------*
      subroutine print_symmetrise2(old_result, item)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &     old_result
      type(itf_contr), intent(inout) ::
     &     item

      character(len=70) ::
     &     line
      character(len=INDEX_LEN) ::
     &     tindex
      character(len=MAXLEN_BC_LABEL) ::
     &     new

      if (item%inter(3)) return
      if (item%rank3<=2) return

      new = rename_tensor(old_result, item%rank3)

      line = '.'//trim(new)//'['//trim(item%idx3)//'] += '//
     &       trim(item%label_res)//'['//trim(item%idx3)//']'
      write(item%logfile,'(a)') trim(line)

      if (item%product) then
         ! We don't need to symmetrise from a tensor product
         return
      end if

      tindex = ' '
      tindex(1:1) = item%idx3(2:2)
      tindex(2:2) = item%idx3(1:1)
      tindex(3:3) = item%idx3(4:4)
      tindex(4:4) = item%idx3(3:3)

      line = '.'//trim(new)//'['//trim(item%idx3)//'] += '//
     &       trim(item%label_res)//'['//trimal(tindex)//']'
      write(item%logfile,'(a)') trim(line)


      return
      end


*----------------------------------------------------------------------*
      subroutine check_symmetric(contr_info, command, symmetric)
*----------------------------------------------------------------------*
!     Check if a tensor has permutational symmetry
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &   contr_info     ! Information about binary contraction
      integer, intent(in) ::
     &   command

      integer ::
     &   c(ngastp,2),
     &   e1(ngastp,2),
     &   e2(ngastp,2),
     &   nops(ngastp),
     &   i
      logical ::
     &   symmetric

      call itf_ops(contr_info, c, e1, e2, command)

      nops = sum(e1, dim=2) + sum(e2, dim=2)

      if (sum(nops)==0) then
         ! Scalars don't have symmetry
         symmetric = .false.
         return
      end if

      symmetric = .true.
      do i = 1, ngastp
         if (mod(nops(i),2) /= 0) then
            symmetric = .false.
         end if
      end do

      return
      end

*----------------------------------------------------------------------*
      subroutine check_j_integral(e, c, nops, j_int, second)
*----------------------------------------------------------------------*
!     Check if the two electron integrals need to be mapped to the J
!     array in Molpro
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   e(ngastp,2),
     &   c(ngastp,2),
     &   nops(ngastp)
      logical, intent(in) ::
     &   second
      logical, intent(inout) ::
     &   j_int

      integer ::
     &   i, j,
     &   ct(ngastp,2)

      j_int = .true.

      if (nops(1)==1 .and. nops(2)==1 .and. nops(3)==2) then
         ! K:eaca - K:eaac (there is no J:eaac)
         j_int = .false.
         return
      end if

      if (second) then
         do i = 1, ngastp
            ct(i,1) = c(i,2)
            ct(i,2) = c(i,1)
         end do
      else
         ct = c
      end if

      do i = 1, ngastp
         do j = 1, 2
            if (e(i,j) + ct(i,j) > 1) then
               j_int = .false.
            end if
         end do
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine itf_rank(ops1, ops2, rank, nops, flag)
*----------------------------------------------------------------------*
!     Calculate rank of tensor
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   ops1(ngastp,2),
     &   ops2(ngastp,2)
      integer, intent(inout) ::
     &   rank,
     &   nops(ngastp)
      logical, intent(in) ::
     &   flag

      if (flag) then
         ! Only use the first occupation array to calculate the rank
         nops = sum(ops1, dim=2)
         rank = sum(nops, dim=1)
      else
         nops = sum(ops1, dim=2) + sum(ops2, dim=2)
         rank = sum(nops, dim=1)
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine itf_ops(contr_info, c, e1, e2, command)
*----------------------------------------------------------------------*
!     Assign contraction (c), external indicies 1 (e1) and external
!     indicies 2 (e2) to ift_contr item
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'mdef_operator_info.h' ! For def_formular_item.h
      include 'def_contraction.h'
      include 'def_formula_item.h' ! For command parameters
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &   contr_info     ! Information about binary contraction

      integer, intent(inout) ::
     &   c(ngastp,2),        ! Operator numbers of contraction index
     &   e1(ngastp,2),       ! Operator numbers of external index 1
     &   e2(ngastp,2)        ! Operator numbers of external index 2

      integer, intent(in) ::
     &   command

      integer ::
     &   i

      c=0
      e1=0
      e2=0

      if (command==command_cp_intm .or.
     &    command==command_add_intm) then
         do i = 1, contr_info%nj_op1
           call count_index2(contr_info%occ_op1(1:,1:,i), e1)
         end do
      else
         ! Get occupation info
         do i = 1, contr_info%n_cnt
           call count_index2(contr_info%occ_cnt(1:,1:,i), c)
         end do
         do i = 1, contr_info%nj_op1
           call count_index2(contr_info%occ_ex1(1:,1:,i), e1)
         end do
         do i = 1, contr_info%nj_op2
           call count_index2(contr_info%occ_ex2(1:,1:,i), e2)
         end do
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine itf_equiv_lines_factor(ops, factor)
*----------------------------------------------------------------------*
!     Figure out factor from equivalent lines
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   ops(ngastp,2)     ! Matrix of contraction indicies
      real(8), intent(inout) ::
     &   factor

      integer ::
     &   i, j

      do j = 1, 2
         do i = 1, ngastp
            if (ops(i,j) == 0) cycle
            if (mod(ops(i,j),2) == 0) then
               factor = factor * (1.0d+0/real(ops(i,j),8))
            end if
         end do
      end do

      return
      end


!     Old alternative to nicer_pairing
!      write(item%logfile,*) "T1 old: {", str1%str, "}", str1%cnt_poss
!      write(item%logfile,*) "T2 old: {", str2%str, "}", str2%cnt_poss
!      write(item%logfile,*) "Result old: {", str3%str, "}"
!
!
!
!      ! TODO: new idea instead of nicer pairing
!      ! 1. Arrange the pairs in the strings (with or without contraction
!      ! indicies)
!      ! 2. Place in remaining contraction indices, making sure creation
!      ! and annhilation are not swapped
!      allocate(t_str1%str(item%rank1))
!      allocate(t_str2%str(item%rank2))
!      allocate(t_str1%cnt_poss(n_cnt))
!      allocate(t_str2%cnt_poss(n_cnt))
!      allocate(t_str1%itype(item%rank1))
!      allocate(t_str2%itype(item%rank2))
!
!      allocate(t1cnt_poss(n_cnt))
!      allocate(t2cnt_poss(n_cnt))
!
!      t1cnt_poss=str1%cnt_poss
!      t2cnt_poss=str2%cnt_poss
!
!      write(10,*) "t1cnt_poss ", t1cnt_poss
!      write(10,*) "t2cnt_poss ", t2cnt_poss
!
!      ! Loop over pairs
!      sh1 = 1
!      sh2 = 1
!      sh3 = 1
!      t_str1%str=''
!      t_str2%str=''
!      do i = 1, item%rank3/2
!         found = .false.
!
!         ! Place in creation index
!         if (p_list%plist(i)%ops(1)==1) then
!            t_str1%str(sh1) = p_list%plist(i)%pindex(1)
!
!            ! Place in annhilation index
!            if (p_list%plist(i)%ops(2)==1) then
!               pp1 = item%rank1 - sh1 + 1
!               t_str1%str(pp1) = p_list%plist(i)%pindex(2)
!            else
!               pp2 = item%rank2 - sh2 + 1
!               t_str2%str(pp2) = p_list%plist(i)%pindex(2)
!
!
!
!               ! Need to place in a contraction annhilation into t1 and
!               ! a contraction creation into t2
!               ! TODO: Factor this
!               found_cnt = .false.
!               do j = 1, n_cnt
!                  if (t1cnt_poss(j)>item%rank1/2) then
!                     found_cnt = .true.
!                     write(10,*) "cnt ", t1cnt_poss(j)
!                     pl1 = str1%cnt_poss(j)
!                     t_str1%str(item%rank1-sh1+1) = str1%str(pl1)
!                     t_str2%str(sh2) = str1%str(pl1)
!
!                     ! set new cnt_poss to t_str1 and t_str2
!                     t_str1%cnt_poss(sh3) = item%rank1-sh1+1
!                     t_str2%cnt_poss(sh3) = sh2
!                     sh3 = sh3 + 1
!
!                     ! set cnt_poss to 0?, so don't get same op twice (in both
!                     ! tensors)
!                     t1cnt_poss(j) = 0
!                     ! Check same contraction exists on t2
!                     do k = 1, item%rank2/2
!                        if (str1%str(pl1)==str2%str(k)) then
!                           found = .true.
!
!                           do l = 1, n_cnt
!                              if (t2cnt_poss(l)==k) then
!                                 t2cnt_poss(l) = 0
!                                 exit
!                              end if
!                           end do
!
!                           exit
!                        end if
!                     end do
!                     if (.not. found) then
!                        write(item%logfile,*) "ERROR: Couldn't find"//
!     &                                        " matching creation on t2"
!                     end if
!
!                     ! Only place one contraction index at a time, exit
!                     ! the loop
!                     exit
!
!                  end if
!               end do
!
!               if (.not. found_cnt) then
!                  write(item%logfile,*) "ERROR: Couldn't find"//
!     &                                  " contraction annhilation on t1"
!               end if
!
!               ! Placed an index in t_str2, so update shift
!               sh2 = sh2 + 1
!
!            end if
!
!            ! Placed an index in t_str1, so update shift
!            sh1 = sh1 + 1
!
!         else
!            pp2 = item%rank2 - shift + 1
!            t_str2%str(sh2) = p_list%plist(i)%pindex(1)
!
!            if (p_list%plist(i)%ops(2)==1) then
!               pp1 = item%rank1 - sh1 + 1
!               t_str1%str(pp1) = p_list%plist(i)%pindex(2)
!
!               ! Need to place in a contraction annhilation into t2 and
!               ! a contraction creation into t1
!               found_cnt = .false.
!               do j = 1, n_cnt
!                  if (t2cnt_poss(j)>item%rank2/2) then
!                     found_cnt = .true.
!                     pl1 = str2%cnt_poss(j)
!                     t_str1%str(sh1) = str2%str(pl1)
!                     t_str2%str(item%rank2-sh2+1) = str2%str(pl1)
!
!                     ! set new cnt_poss to t_str1 and t_str2
!                     t_str1%cnt_poss(sh3) = sh1
!                     t_str2%cnt_poss(sh3) = item%rank2-sh2+1
!                     sh3 = sh3 + 1
!
!                     ! set cnt_poss to 0?, so don't get same op twice (in both
!                     ! tensors)
!                     t2cnt_poss(j) = 0
!                     ! Check same contraction exists on t1
!                     do k = 1, item%rank1/2
!                        if (str2%str(pl1)==str1%str(k)) then
!                           found = .true.
!
!                           do l = 1, n_cnt
!                              if (t1cnt_poss(l)==k) then
!                                 t1cnt_poss(l) = 0
!                                 exit
!                              end if
!                           end do
!
!                        end if
!                     end do
!                     if (.not. found) then
!                        write(item%logfile,*) "ERROR: Couldn't find"//
!     &                                        " matching creation on t1"
!                     end if
!
!                     ! Only place one contraction index at a time, exit
!                     ! the loop
!                     exit
!
!                  end if
!               end do
!
!               if (.not. found_cnt) then
!                  write(item%logfile,*) "ERROR: Couldn't find"//
!     &                                  " contraction annhilation on t2"
!               end if
!
!               ! Placed an index in t_str1, so update shift
!               sh1 = sh1 + 1
!
!
!            else
!               pp2 = item%rank2 - sh2 + 1
!               t_str2%str(pp2) = p_list%plist(i)%pindex(2)
!            end if
!
!            ! Placed an index in t_str2, so update shift
!            sh2 = sh2 + 1
!
!         end if
!      end do
!
!      ! 2.5 Place in remaing contration loops
!      ! 3. Update old stings and update contraction postitions
!      ! TODO: decide which one is bigger, then only search on that
!      ! tensor
!      if (sh1-1 < item%rank1/2) then
!         do i = 1, n_cnt
!            if (t1cnt_poss(i)>0) then
!               if (t1cnt_poss(i)>item%rank1/2) then
!                  t_str1%str(item%rank1-sh1+1) = str1%str(t1cnt_poss(i))
!
!                  do j = 1, n_cnt
!                     if (t1cnt_poss(j)>0 .and.
!     &                                 t1cnt_poss(j)<=item%rank1/2) then
!
!                        ! Place in creation contraction
!                        t_str1%str(sh1) = str1%str(t1cnt_poss(j))
!                        t_str2%str(item%rank2-sh2+1) =
!     &                                           str1%str(t1cnt_poss(j))
!                        t1cnt_poss(j) = 0
!                     end if
!                  end do
!
!                  t_str2%str(sh2) = str1%str(t1cnt_poss(i))
!                  sh1 = sh1 + 1
!                  sh2 = sh2 + 1
!
!                  ! Update contration positions
!                  t_str1%cnt_poss(sh3) = sh1
!                  t_str2%cnt_poss(sh3) = sh2
!                  sh3 = sh3 + 1
!                  t_str1%cnt_poss(sh3) = item%rank1-sh1+1
!                  t_str2%cnt_poss(sh3) = item%rank2-sh2+1
!
!                  t1cnt_poss(i) = 0
!               end if
!            end if
!         end do
!      end if
!
!
!      write(10,*) "Some str1 {", t_str1%str, "} ", t_str1%cnt_poss
!      write(10,*) "Some str2 {", t_str2%str, "} ", t_str2%cnt_poss
!
!      !str1%str = t_str1%str
!      !str2%str = t_str2%str
!      !str1%cnt_poss = t_str1%cnt_poss
!      !str2%cnt_poss = t_str2%cnt_poss
!
!
!      deallocate(t1cnt_poss)
!      deallocate(t2cnt_poss)
!
!      deallocate(t_str1%str)
!      deallocate(t_str2%str)
!      deallocate(t_str1%cnt_poss)
!      deallocate(t_str2%cnt_poss)
!      deallocate(t_str1%itype)
!      deallocate(t_str2%itype)

