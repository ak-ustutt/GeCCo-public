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

      write(logfile,*)
      write(logfile,*) "SPIN OF ", trim(label)
      write(logfile,*) "---------------------------------"
      write(logfile,*) (spins(2,i), i=1, rank/2)
      write(logfile,*) (spins(1,i), i=1, rank/2)
      write(logfile,*) "---------------------------------"
      write(logfile,*)

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

      write(logfile,*)
      write(logfile,*) "P_list: ", label
      write(logfile,*) "------------------------------"
      do i = 1, size
      write(logfile,*) p_list%plist(i)%pindex(1), p_list%plist(i)%ops(1)
      write(logfile,*) p_list%plist(i)%pindex(2), p_list%plist(i)%ops(2)
      write(logfile,*) p_list%plist(i)%link
      write(logfile,'(3i3)') (p_list%plist(i)%nval(j), j=1, 3)
      write(logfile,*) "---------------------------"
      end do
      write(logfile,*)

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

      write(item%out,*) "ERROR: ", error
      write(item%out,*) "================================="
      write(item%out,*) "Result: ", item%label_res, trim(item%idx3)
      write(item%out,*) "Tensor1: ", item%label_t1, trim(item%idx1)
      write(item%out,*) "Tensor2: ", item%label_t2, trim(item%idx2)
      write(item%out,*)

      return
      end


*----------------------------------------------------------------------*
      subroutine debug_header(label, logfile)
*----------------------------------------------------------------------*
!     Print error message in the itf logfile
*----------------------------------------------------------------------*

      implicit none

      character(len=*), intent(in) ::
     &   label
      integer, intent(in) ::
     &   logfile

      write(logfile,*)
      write(logfile,*) "================================="
      write(logfile,*) "Debug from:"
      write(logfile,*)
      write(logfile,*) trim(label)
      write(logfile,*) "================================="
      write(logfile,*)

      return
      end


*----------------------------------------------------------------------*
      subroutine count_index(iocc, nops)
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
      subroutine set_itype(item, itype, itype2, ntest)
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

      type(itf_contr), intent(inout) ::
     &   item        ! ITF contraction object; holds all info about the ITF algo line
      integer, intent(inout) ::
     &   itype(INDEX_LEN)
      integer, intent(inout) ::
     &   itype2(MAXINT, INDEX_LEN)
      integer, intent(in) ::
     &   ntest

      integer ::
     &   i, j,
     &   ninter

      ! If the current result is not an intermeidte, don't need an itype
      ! for the next line
      if (.not. item%inter(3)) then
         itype = 0
         return
      end if

      if (item%rank3==2) then
         itype(1) = get_itype(item%idx3(1:1))
         itype(2) = get_itype(item%idx3(2:2))
      else if (item%rank3==4) then
         itype(1) = get_itype(item%idx3(1:1))
         itype(2) = get_itype(item%idx3(2:2))
         itype(3) = get_itype(item%idx3(4:4))
         itype(4) = get_itype(item%idx3(3:3))
      else if (item%rank3==6) then
         itype(1) = get_itype(item%idx3(1:1))
         itype(2) = get_itype(item%idx3(2:2))
         itype(3) = get_itype(item%idx3(3:3))
         itype(4) = get_itype(item%idx3(6:6))
         itype(5) = get_itype(item%idx3(5:5))
         itype(6) = get_itype(item%idx3(4:4))
      end if

      if (ntest>=100) then
         call debug_header("set_itype", item%out)
         write(item%out,*) "itype from intermediate: ", item%idx3
         write(item%out,'(a6)',advance='no') "     {"
         do i = 1, item%rank3
            write(item%out,'(i0, 1x)', advance='no') itype(i)
         end do
         write(item%out,'(a1)') "}"
      end if


      ! More advanced itype for handling n intermeidates which are
      ! defined and used anywhere within the complete diagram
      i = len(trim(item%label_res))
      read(item%label_res(i:i),'(1i)') ninter

      ! Clear all itype info from previous diagram
      if (i == 1) then
         itype2 = 0
      end if

      if (item%rank3==2) then
         itype2(ninter,1) = get_itype(item%idx3(1:1))
         itype2(ninter,2) = get_itype(item%idx3(2:2))
      else if (item%rank3==4) then
         itype2(ninter,1) = get_itype(item%idx3(1:1))
         itype2(ninter,2) = get_itype(item%idx3(2:2))
         itype2(ninter,3) = get_itype(item%idx3(4:4))
         itype2(ninter,4) = get_itype(item%idx3(3:3))
      else if (item%rank3==6) then
         itype2(ninter,1) = get_itype(item%idx3(1:1))
         itype2(ninter,2) = get_itype(item%idx3(2:2))
         itype2(ninter,3) = get_itype(item%idx3(3:3))
         itype2(ninter,4) = get_itype(item%idx3(6:6))
         itype2(ninter,5) = get_itype(item%idx3(5:5))
         itype2(ninter,6) = get_itype(item%idx3(4:4))
      end if


      if (ntest>=100) then
         call debug_header("set_itype", item%out)
         do j = 1, MAXINT
            write(item%out,'(a6)',advance='no') "     {"
            do i = 1, item%rank3
               write(item%out,'(i0, 1x)', advance='no') itype2(j, i)
            end do
            write(item%out,'(a1)') "}"
         end do
      end if


      return
      end


*----------------------------------------------------------------------*
      subroutine command_to_itf(contr_info, itin, itflog, command,
     &                           inter_itype, contr_no, nk4e, itype2)
*----------------------------------------------------------------------*
!     Take GeCco binary contraction and produce ITF algo code.
!     Includes antisymmetry of residual equations and spin summation.
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(inout) ::
     &   contr_info          ! Information about binary contraction
      logical, intent(in) ::
     &   itin                ! Print ITIN lines or not
      integer, intent(in) ::
     &   itflog,             ! Output file
     &   command,            ! Type of formula item command, ie. contraction, copy etc.
     &   contr_no            ! Formula number
      integer, intent(inout) ::
     &   nk4e                ! K4E counter
      integer, intent(inout) ::
     &   inter_itype(*)      ! Store itypes of intermediates between lines
      integer, intent(inout) ::
     &   itype2(MAXINT, INDEX_LEN)      ! Store itypes of intermediates between lines

      type(itf_contr) ::
     &   item,               ! ITF contraction object; holds all info about the ITF algo line
     &   pitem               ! Permutation ITF contraction object
      integer ::
     &   i,                  ! Loop index
     &   perm_case,          ! Info of permutation factors
     &   ntest1 =  00,       ! Control debug: init_itf_contr
     &   ntest2 =  00,       ! Control debug: assign_index
     &   ntest3 =  00,       ! Control debug: determine permuation
     &   ntest4 =  00,       ! Control debug: prepare_symmetrise
     &   ntest5 =  00,       ! Control debug: prepare permutation
     &   ntest6 =  00,       ! Control debug: assign_spin
     &   ntest7 =  00,       ! Control debug: print_spin_cases
     &   ntest8 =  00,       ! Control debug: print_symmetrise
     &   ntest9 =  00,       ! Control debug: set_itype
     &   ntest10 = 00        ! Control debug
      logical ::
     &   intpp,              ! Use INTPP kext intermeidate
     &   pline               ! True if including a permuation line

      if (ntest1>=100 .or. ntest2>=100) then
         write(itflog,*) "FORMULA NUMBER: ", contr_no
      end if

      ! Begin a special block which the python processor will pull out
      ! into its own code block
      intpp = .false.
      if (contr_info%label_res=='INTpp') then
         intpp = .true.
         write(itflog,'(a)') "BEGIN_INTPP"
      end if

      ! Mark begining of spin summed block
      write(itflog,'(a5)') 'BEGIN'


      ! 1. Initalise itf_contr
      call itf_contr_init(contr_info,item,1,itin,command,itflog,
     &                    inter_itype,itype2,nk4e,ntest1)


      ! 2. Assign index / Determine sign
      call assign_index(item,ntest2)


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

      if (ntest3>=100) then
         call debug_header("Determine permutation", item%out)
         write(item%out,*) "perm_case: ", perm_case
      end if


      ! 4. Decide whether to symmetrise after every term
      item%old_name = item%label_res
      call prepare_symmetrise(perm_case, item, ntest4)



      ! 5. Prepare permuation line if required
      if (item%symmetric .and. perm_case==2 .or.
     &    .not. item%symmetric .and. perm_case==1 .and.
     &    .not. item%nosym) then

         ! Don't need perm line for intermediates or tensor products
         if (.not. item%inter(3) .and. .not. item%product) then
            pline = .true.

            call itf_contr_init(contr_info,pitem,1,itin,command,itflog,
     &                           inter_itype,itype2,nk4e,ntest5)
            call assign_index(pitem,ntest5)

            pitem%old_name = pitem%label_res
            pitem%label_res = item%label_res

            ! Permute indicies and update factor
            ! -1 factor for permuation line
            call create_permutation(pitem, contr_info%perm,ntest5)
         end if

      end if

      if (item%nosym) then
         ! Need to get both ab/ab and ab/ba spin cases for R[apiq]
         ! For now, this is stored in pline as it is not used in R[apiq]
         pline = .true.

         call itf_contr_init(contr_info,pitem,1,itin,command,itflog,
     &                       inter_itype,itype2,nk4e,ntest5)
         call assign_index(pitem,ntest5)

         pitem%old_name = pitem%label_res
         pitem%label_res = item%label_res
         pitem%abba_line = .true.
      end if


      ! 6. Spin sum
      call assign_spin(item, ntest6)
      if (pline) call assign_spin(pitem, ntest6)


      ! 7. Loop over spin cases and print out each line
      call print_spin_cases(item, ntest7)
      if (pline) call print_spin_cases(pitem, ntest7)


      ! 8. Print symmetrisation term
      if (.not. item%inter(3)) then
         call print_symmetrise(item, ntest8)
      end if


      ! 9. If an intermediate, set inter_itype for use in next line
      ! where the intermediate is created
      call set_itype(item, inter_itype, itype2, ntest9)


      if (ntest10>=100) then
         call print_itf_contr(item)
         if (pline) call print_itf_contr(pitem)
      end if


      ! Mark end of spin block
      write(itflog,'(a)') "END"

      if (intpp) then
         write(itflog,'(a)') "END_INTPP"
      end if


      ! Update K4E counter
      if (item%k4e_line) then
         nk4e = nk4e + 1
      end if


      ! Deallocate memroy used when construcitng item
      call itf_deinit(item)

      return
      end


*----------------------------------------------------------------------*
      subroutine itf_contr_init(contr_info,item,perm,itin,comm,lulog,
     &                          itype,itype2,nk4e,ntest)
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
     &   contr_info   ! Information about binary contraction
      type(itf_contr), intent(inout) ::
     &   item     ! Object which holds information necessary to print out an ITF algo line
      integer, intent(in) ::
     &   perm,        ! Permutation information
     &   comm,        ! formula_item command
     &   lulog,        ! Output file
     &   itype(INDEX_LEN),
     &   itype2(MAXINT,INDEX_LEN),
     &   nk4e,
     &   ntest
      logical, intent(in) ::
     &   itin

      ! Assign output file
      item%out=lulog

      ! Assign command type
      item%command=comm

      ! Initalise index strings
      item%idx1 = ''
      item%idx2 = ''
      item%idx3 = ''

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

      call check_symmetric(contr_info, item%command, item%symmetric,
     &                     item%nosym)

      ! If a residual, is it symmetric (R_{ab}^{ij] = R_{ba}^{ji})?
      if (.not. item%inter(3)) then
         call check_symmetric(contr_info, item%command, item%symm_res,
     &                        item%nosym)

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
      end if


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
      item%itype2 = itype2

      ! Set K4E info
      item%nk4e = nk4e
      if (item%int(1) .or. item%int(2))then
         call check_k4e(contr_info, item%command, item%k4e_line)
      end if

      if (ntest>=100) then
         call debug_header("itf_contr_init", item%out)
         call print_itf_contr(item)
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

      return
      end



*----------------------------------------------------------------------*
      subroutine create_permutation(item, perm, ntest)
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
      integer, intent(in) ::
     &   ntest

      integer ::
     &   ex_itype,
     &   i, j,
     &   shift
      character(len=1) ::
     &   ex_ind(2)
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

      if (shift < 3) write(item%out,*) "ERROR in create_permuation"

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

      ! Update orginal index
      item%idx1 = tmp1
      item%idx2 = tmp2

      ! Negative sign from (1-P_xy)
      item%fact = item%fact * -1.0d+0

      if (ntest>=100) then
         call debug_header("create_permuation", item%out)
         write(item%out,*) "external index: ", ex_ind
         write(item%out,*) "external itype: ", ex_itype
         write(item%out,*) "new factor: ", item%fact
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine convert_to_abab_block(item, t_spin, new_idx1, new_idx2,
     &                                 new_idx3, new_fact)
*----------------------------------------------------------------------*
!     Convert integrals and amplitudes to abab spn blocks, this may also
!     intoduce a sign change
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &   item
      type(spin_info2), intent(in) ::
     &      t_spin(3)
      character(len=INDEX_LEN), intent(inout) ::
     &     new_idx1, new_idx2, new_idx3
      real(8), intent(inout) ::
     &     new_fact

      character(len=1) ::
     &   tmp
      integer ::
     &   ops(4,2)

      ! TODO: Don't bother for eeee or cccc integrals - a hack for now
      ! to allow the summation of K:eeee terms in simplfy.py
      !ops = item%e1 + item%c
      !if (ops(2,2)==2 .and. ops(2,1)==2) then
      !   if (.not. item%inter(1)) then
      !      return
      !   end if
      !end if

      !if (ops(1,2)==2 .and. ops(1,1)==2) then
      !   if (.not. item%inter(1) .and. .not. item%inter(3)) then
      !      return
      !   end if
      !end if

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

      ! For the R[apiq] abba case, we need to rearange
      if (item%rank3>2 .and. item%abba_line) then
         if (t_spin(3)%spin(1,1)>
     &       t_spin(3)%spin(1,2)) then
            tmp = new_idx3(2:2)
            new_idx3(2:2) = new_idx3(1:1)
            new_idx3(1:1) = tmp
            new_fact = new_fact * -1.0d+0
         end if
         if (t_spin(3)%spin(2,1)>
     &       t_spin(3)%spin(2,2)) then
            tmp = new_idx3(3:3)
            new_idx3(3:3) = new_idx3(4:4)
            new_idx3(4:4) = tmp
            new_fact = new_fact * -1.0d+0
         end if
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine prepare_symmetrise(perm_case, item, ntest)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   perm_case,                     ! Info of permutation factors
     &   ntest                          ! Control debug
      type(itf_contr), intent(inout) ::
     &   item                           ! ITF contration info

      ! Return if an intermediate or result rank less than 2
      if (item%inter(3) .or. item%rank3<=2) return

      ! Return if not a symmetric result
      if (.not. item%symmetric) return

      ! Return if INTPP block
      if (item%intpp) return

      if (ntest>=100) then
         call debug_header("prepare_symmetrise", item%out)
         write(item%out,*) "Old factor: ", item%fact
         write(item%out,*) "Old label: ", item%label_res
      end if

      ! Check permuational symmetry, if permutational symmetry, then
      ! multiply by 0.5 (we will symmetrise each result)
      if (perm_case == 0) item%fact = item%fact * 0.5d+0

      ! Introduce ITIN intermeidate to collect terms
      ! .R[abij] += ITIN[abij]
      ! .R[abij] += ITIN[baji]
      item%old_name = item%label_res
      item%label_res = "ITIN"

      if (ntest>=100) then
         call debug_header("prepare_symmetrise", item%out)
         write(item%out,*) "New factor: ", item%fact
         write(item%out,*) "New label: ", item%label_res
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine print_spin_cases(item, ntest)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &   item
      integer, intent(in) ::
     &   ntest

      integer ::
     &   i, j,
     &   actual_spin_cases
      logical ::
     &   contains1, contains2,
     &   s1, s2

      if (ntest>=100) call debug_header("print_spin_cases", item%out)

      ! Cover case where only 1 spin case per line
      !if (item%nspin_cases==1) then
      !   actual_spin_cases = 2
      !   write(item%out,*) "hello"
      !else
      !   actual_spin_cases = item%nspin_cases
      !end if

      actual_spin_cases = item%nspin_cases

      do j = 1, actual_spin_cases -1

            ! Decide if tensor is mixed spin
            contains1 = .false.
            contains2 = .false.
            s1 = .false.
            if (item%rank1>2) then
               do i = 1, item%rank1/2
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

            ! True if pure spin
            contains1 = .false.
            contains2 = .false.
            s2 = .false.
            if (item%rank2>2) then
               do i = 1, item%rank2/2
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

            if (ntest>=100) then
               write(item%out,*) "T1 is pure spin: ", s1
               write(item%out,*) "T2 is pure spin: ", s2
            end if

            call print_itf_line(item,s1,s2,item%all_spins(j)%t_spin,
     &                          ntest)

      end do

      return
      end

*----------------------------------------------------------------------*
      subroutine print_itf_line(item,s1,s2,t_spin,ntest)
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
     &   item
      logical, intent(in) ::
     &   s1,s2
      type(spin_info2), intent(in) ::
     &   t_spin(3)
      integer, intent(in) ::
     &   ntest

      character(len=MAXLEN_BC_LABEL) ::
     &   nres, nt1, nt2          ! Name of tensors involved in the contraction
      character(len=INDEX_LEN) ::
     &   slabel1, slabel2, slabel3,
     &   new_idx1, new_idx2, new_idx3
      character(len=264) ::
     &   itf_line,          ! Line of ITF code
     &   st1, st2           ! Name of spin summed tensors + index
      character(len=2) ::
     &   equal_op           ! ITF contraction operator; ie. +=, -=, :=
      character(len=25) ::
     &   sfact,             ! String representation of factor
     &   sfact_star,        ! String representation of factor formatted for output
     &   k4e_no             ! Counter of K4E tensors
      integer ::
     &   i
      real(8) ::
     &   c_fact               ! Copy of orginal factor
      logical ::
     &   new_j

      new_idx1 = item%idx1
      new_idx2 = item%idx2
      new_idx3 = item%idx3
      c_fact = item%fact
      new_j = item%j_int

      if (ntest>=100) then
         call debug_header("print_itf_line", item%out)
         write(item%out,*) "Idx1: ", new_idx1
         write(item%out,*) "Idx2: ", new_idx2
         write(item%out,*) "Idx3: ", new_idx3
         write(item%out,*) "Factor: ", c_fact
         write(item%out,*) "j_int: ", new_j
      end if

      ! Reorder tensor index into abab blocks
      ! May get factor change here
      call convert_to_abab_block(item, t_spin, new_idx1, new_idx2,
     &                           new_idx3, c_fact)

      if (ntest>=100) then
         write(item%out,*) "---------------------------"
         write(item%out,*) "After convert_to_abab_block"
         write(item%out,*) "Idx1: ", new_idx1
         write(item%out,*) "Idx2: ", new_idx2
         write(item%out,*) "Idx3: ", new_idx3
         write(item%out,*) "Factor: ", c_fact
         write(item%out,*) "---------------------------"
      end if

      ! Change names of specific tensors
      nres=rename_tensor(item%label_res, item%rank3)
      nt1=rename_tensor(item%label_t1, item%rank1)
      nt2=rename_tensor(item%label_t2, item%rank2)

      ! Add intermediate spin strings to names
      if (item%inter(1)) then
         call inter_spin_name(t_spin(1)%spin,item%rank1/2,slabel1)
         nt1 = trim(nt1)//trim(slabel1)
         call reorder_intermediate(item%rank1,new_idx1,item%nops1)
      end if
      if (item%inter(2)) then
         call inter_spin_name(t_spin(2)%spin,item%rank2/2,slabel2)
         nt2 = trim(nt2)//trim(slabel2)
         call reorder_intermediate(item%rank2,new_idx2,item%nops2)
      end if
      if (item%inter(3)) then
         call inter_spin_name(t_spin(3)%spin,item%rank3/2,slabel3)
         nres = trim(nres)//trim(slabel3)
         call reorder_intermediate(item%rank3,new_idx3,item%nops3)
      end if

      ! Reorder integrals into fixed slot order
      call reorder_integral(item%int(1),item%rank1,new_idx1,
     &                      new_j,
     &                      nt1,item%nops1)
      call reorder_integral(item%int(2),item%rank2,new_idx2,
     &                      new_j,
     &                      nt2,item%nops2)

      if (.not. item%inter(1) .and. .not. item%int(1)) then
         call reorder_amp(item%rank1, new_idx1)
      end if
      if (.not. item%inter(2) .and. .not. item%int(2)) then
         call reorder_amp(item%rank2, new_idx2)
      end if
      if (.not. item%inter(3) .and. .not. item%int(3)) then
         call reorder_amp(item%rank3, new_idx3)
      end if

      if (ntest>=100) then
         write(item%out,*) "---------------------------"
         write(item%out,*) "After reorder()"
         write(item%out,*) "Idx1: ", new_idx1
         write(item%out,*) "Idx2: ", new_idx2
         write(item%out,*) "Idx3: ", new_idx3
         write(item%out,*) "Factor: ", c_fact
         write(item%out,*) "---------------------------"
      end if

      if (ntest>=200) then
         call print_spin(t_spin(1)%spin, item%rank1, "T1", item%out)
         call print_spin(t_spin(2)%spin, item%rank2, "T2", item%out)
         write(item%out,*) "s1: ", s1
         write(item%out,*) "s2: ", s2
      end if

      ! Change tensor to spatial orbital quantity, unless it is an
      ! intermediate
      call spatial_string(st1,new_idx1,nt1,s1,item%inter(1),item%rank1,
     &                1,item%binary,item%int(1),item%nops1,new_j,
     &                item%out)
      call spatial_string(st2,new_idx2,nt2,s2,item%inter(2),item%rank2,
     &                2,item%binary,item%int(2),item%nops2,new_j,
     &                item%out)


      ! Add factor to scalar result cases (going to skip half the spin
      ! cases as these are the same, so add a factor of two to the
      ! remaining ones)
      if (item%rank3 == 0 .and. item%rank1/=0) then
         c_fact = c_fact * 2.0d0
         if (ntest>=100) write(item%out,*) "Changing factor: ", c_fact
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
!     &    '['//trim(item%idx3)//'] '//equal_op//' '//
     &    '['//trim(new_idx3)//'] '//equal_op//' '//
     &    trim(sfact_star)//trimal(st1)//' '//trimal(st2)


      ! Handle K4E
      ! Print orginal line as a comment, then line with K4E
      if (item%k4e_line) then
         !itf_line = "// "//trim(itf_line)
         !write(item%out,'(a)') trim(itf_line)

         write(k4e_no,*) item%nk4e

         itf_line='.'//trimal(nres)//
     &      '['//trim(new_idx3)//'] '//equal_op//' '//
     &      trim(sfact_star)//'K4E'//trim(adjustl(k4e_no))//
     &      '['//trim(new_idx3)//']'
      end if

      ! Print it to bcontr.tmp
      write(item%out,'(a)') trim(itf_line)

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
      subroutine reorder_integral(integral,rank,idx,j_int,label,
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

      ! TODO: is this the same as reorder amplitude?

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
      subroutine reorder_intermediate(rank,idx,nops)
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

      integer ::
     &   itype(rank),
     &   i,
     &   itmp
      character(len=1) ::
     &   tmp
      logical ::
     &   permute

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
                  else if ((nops(1)==1.and.nops(2)==1.and.
     &                      nops(3)==2)) then
                    ! T:eaca or T:eaac
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
      subroutine assign_index(item,ntest)
*----------------------------------------------------------------------*
!     Assign an ITF index string to each tensor in a line
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &   item           ! ITF binary contraction
      integer, intent(in) ::
     &   ntest

      integer ::
     &   c(ngastp,2),        ! Operator numbers of contraction index
     &   ci(ngastp,2),       ! Operator numbers of contraction index (inverse)
     &   e1(ngastp,2),       ! Operator numbers of external index 1
     &   e2(ngastp,2),       ! Operator numbers of external index 2
     &   shift,         ! List shift
     &   n_cnt,         ! Number of contraction operators
     &   e1ops, e2ops,  ! Number of external ops on T1 and T2
     &   distance,      ! Distance from where an index should be
     &   pp, pp2,            ! Paired position - position of paired index
     &   i, j, k, l,
     &   itype(INDEX_LEN)
      character(len=INDEX_LEN) ::
     &   s1, s2, s3,  ! Tmp ITF index strings
     &   tstr
      character(len=1) ::
     &   tmp       ! Scratch space to store index letter
      real(8) ::
     &   p_factor       ! Overall factor from contraction/rearrangment
      type(pair_list) ::
     &   p_list       ! Complete list of pairs in binary contraction
      integer, dimension(4) ::
     &   e_shift,       ! Index shift for external indices
     &   c_shift        ! Index shift for contraction indices
      type(index_str) ::
     &   str1,          ! Stores 'normal ordered' index string for T1
     &   str2,          ! Stores 'normal ordered' index string for T2
     &   str3           ! Stores 'normal ordered' index string for Res

      logical ::        ! These are used when finding pairs of external ops
     &   is_cnt,        ! True if the operator is a contraction op
     &   p1, p2

      if (ntest>=100) call debug_header("assign_index", item%out)

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

      if (ntest>=100) then
         write(item%out,*) "Inital index strings in normal order"
         write(item%out,*) "T1: {", str1%str, "}"
         write(item%out,*) "T2: {", str2%str, "}"
         write(item%out,*) "Res: {", str1%str, "}{", str2%str, "}"
         write(item%out,*) "Contraction T1: ", str1%cnt_poss
         write(item%out,*) "Contraction T2: ", str2%cnt_poss
         write(item%out,*) "Inital factor: ", p_factor
      end if


      ! Rearrange intermediate index to matach previously declared inter
      ! index. This uses the itype array from the previous lines
      ! TODO: itype for two inters on same line
      if (item%inter(1)) then
         call arrange_inter_itype(item%rank1,item%rank2,str1,str2,
     &                           item%itype, item%itype2, item%label_t1)
         if (ntest>=100) then
            write(item%out,*) "After arragne_inter_itpye, inter(1)"
            write(item%out,*) "T1: {", str1%str, "}"
            write(item%out,*) "itpye: {", str1%itype, "}"
            write(item%out,*) "cnt_poss: ", str1%cnt_poss
            write(item%out,*) "previous itype: ", itype
         end if
      end if
      if (item%inter(2)) then
         call arrange_inter_itype(item%rank2,item%rank1,str2,str1,
     &                           item%itype, item%itype2, item%label_t2)
         if (ntest>=100) then
            write(item%out,*) "After arragne_inter_itpye, inter(2)"
            write(item%out,*) "T2: {", str1%str, "}"
            write(item%out,*) "itpye: {", str1%itype, "}"
            write(item%out,*) "cnt_poss: ", str1%cnt_poss
            write(item%out,*) "previous itype: ", itype
         end if
      end if


      ! Work out the factor due to permuation of contraction indicies
      if (.not. item%den(1) .and. .not. item%den(2)) then

         do i = 1, n_cnt
            if (mod(item%rank1-str1%cnt_poss(i)+str2%cnt_poss(i)-1,2)
     &          /=0)then
               p_factor = p_factor * -1.0d0
               if (ntest>=100) then
                  write(item%out,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
               end if
            end if
         end do

      else if (item%den(1) .and. item%rank1/=0) then

         do i = 1, n_cnt
            if (str1%cnt_poss(i)<=item%rank1/2) then
               ! Creation operator of density is contracted
               if (mod(item%rank1-str1%cnt_poss(i)+str2%cnt_poss(i)-1
     &                 -item%rank1/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>=100) then
                     write(item%out,*)"Update factor 1 (contraction of "
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
                  if (ntest>=100) then
                     write(item%out,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
                  end if
               end if

            end if
         end do

      else if (item%den(2) .and. item%rank2/=0) then

         do i = 1, n_cnt
            if (str2%cnt_poss(i)<=item%rank2/2) then
               ! Creation operator of density is contracted
               if (mod(item%rank2-str2%cnt_poss(i)+str1%cnt_poss(i)-1
     &                 -item%rank2/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>=100) then
                     write(item%out,*)"Update factor 1 (contraction of "
     &                             //"contraction indicies): ", p_factor
                  end if
               end if

            else if (str2%cnt_poss(i)>item%rank2/2) then
               ! Annhilation operator of density is contracted
               if (mod(item%rank1-str1%cnt_poss(i)+str2%cnt_poss(i)-1
     &                 -item%rank2/2,2)/=0)then
                  p_factor = p_factor * -1.0d0
                  if (ntest>=100) then
                     write(item%out,*)"Update factor 2 (contraction of "
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

      if (ntest>=100) then
         write(item%out,*) "Index strings in normal order"
         write(item%out,*) "Result string in arbitary order"
         write(item%out,*) "T1: {", str1%str, "}"
         write(item%out,*) "T2: {", str2%str, "}"
         write(item%out,*) "Res: {", str3%str, "}"
         write(item%out,*) "Contraction T1: ", str1%cnt_poss
         write(item%out,*) "Contraction T2: ", str2%cnt_poss
         write(item%out,*) "Inital factor: ", p_factor
      end if

      ! Rearrange the result string so it is in normal order (all
      ! creation operators to the left of the annhilation). This can
      ! also introduce a sign change.
      if (item%rank2/=0) then
         tstr = ""
         do i = 1, item%rank1/2
            shift = 1
            do j = 1, item%rank3

               if (str3%str(j) == str1%str(i)) then
                  tstr(shift:shift) = str3%str(j)
                  shift = shift + 1

                  do k = 1, item%rank3
                     if (str3%str(k) /= str3%str(j)) then
                        tstr(shift:shift) = str3%str(k)
                        shift = shift + 1
                     end if
                  end do

                  do k = 1, item%rank3
                     str3%str(k) = tstr(k:k)
                  end do

                  if (mod(item%rank3-j,2)==0) then
                     p_factor = p_factor * -1.0d0
                     if (ntest>=100) then
                        write(item%out,*) "Update factor 2 (rearrange",
     &                 " the result string to normal order): ", p_factor
                     end if
                  end if
                  exit
               end if

            end do
         end do


         do i = 1, item%rank2/2
            shift = 1
            do j = 1, item%rank3

               if (str3%str(j) == str2%str(i)) then
                  tstr(shift:shift) = str3%str(j)
                  shift = shift + 1

                  do k = 1, item%rank3
                     if (str3%str(k) /= str3%str(j)) then
                        tstr(shift:shift) = str3%str(k)
                        shift = shift + 1
                     end if
                  end do

                  do k = 1, item%rank3
                     str3%str(k) = tstr(k:k)
                  end do

                  if (mod(item%rank3-j,2)==0) then
                     p_factor = p_factor * -1.0d0
                     if (ntest>=100) then
                        write(item%out,*) "Update factor 2 (rearrange",
     &                 " the result string to normal order): ", p_factor
                     end if
                  end if
                  exit
               end if

            end do
         end do
      end if

      ! Due to how the R:ea residual is defined, {p a^+} instead of
      ! {a^+ p}, we need an extra minus to flip the normal ordered
      ! string.
      if (item%nops3(1)==0 .and. item%nops3(2)==1 .and.
     &    item%nops3(3)==1 .and. item%nops3(4)==0) then
         if (.not. item%inter(3)) then

            p_factor = p_factor * -1.0d0
            if (ntest>100) then
               write(item%out,*) "Update factor (R:ea)", p_factor
            end if
         end if
      end if

      if (ntest>=100) then
         write(item%out,*) "Result string in normal order"
         write(item%out,*) "Res: {", str3%str, "}"
      end if


      ! TODO: put this in reorder_integral
      if (item%rank1 == 2) then
         if (str1%itype(1)>str1%itype(2)) then
            tmp = str1%str(1)
            str1%str(1) = str1%str(2)
            str1%str(2) = tmp
         end if
      end if
      if (item%rank2 == 2) then
         if (str2%itype(1)>str2%itype(2)) then
            tmp = str2%str(1)
            str2%str(1) = str2%str(2)
            str2%str(2) = tmp
         end if
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

      if (ntest>=100) then
         if (p_factor>0.0d0) then
            tmp = '+'
         else
            tmp = '-'
         end if

         write(item%out,*) "Final indices"
         write(item%out,*) "---------------------------"
         write(item%out,*) "{",str3%str,"} ",tmp,"= {",str1%str,
     &                         "}{",str2%str,"}"

         write(item%out,*) "[",trim(s3),"] ",tmp,"= [",trim(s1),
     &                         "][",trim(s2),"]"
         write(item%out,*) "---------------------------"
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
      subroutine arrange_inter_itype(rank1, rank2, idx1, idx2, itype,
     &                               itype2, label)
*----------------------------------------------------------------------*
!     Rearrange intermeidate index according to the itype of the
!     previously declared intermediate
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   rank1,
     &   rank2,
     &   itype(INDEX_LEN),
     &   itype2(MAXINT,INDEX_LEN)
      type(index_str), intent(inout) ::
     &   idx1,
     &   idx2
      character(len=maxlen_bc_label), intent(in) ::
     &   label

      integer ::
     &   i, j,
     &   shift,
     &   ninter,
     &   titype(INDEX_LEN)
      character (len=INDEX_LEN) ::
     &   tstr
      logical ::
     &   same

      ! Find itype corresponding to numbered intermediate
      i = len(trim(label))
      read(label(i:i),'(1i)') ninter

      titype = 0
      do i = 1, rank1
         titype(i) = itype2(ninter, i)
      end do

      same = .true.
      do i = 1, rank1
         ! Check if intermediate already has correct index itype
         if (idx1%itype(i) /= itype2(ninter,i)) then
            same = .false.
         end if
      end do

      if (same) return


      do i = 1, rank1/2
         do j = 1, rank1/2
            if (itype2(ninter,i) == idx1%itype(j)) then
               tstr(i:i) = idx1%str(j)
               idx1%str(j) = ''
               idx1%itype(j) = 0
               exit
            end if
         end do
      end do

      do i = rank1/2 + 1, rank1
         do j = rank1/2+1, rank1
            if (itype2(ninter,i) == idx1%itype(j)) then
               tstr(i:i) = idx1%str(j)
               idx1%str(j) = ''
               idx1%itype(j) = 0
               exit
            end if
         end do
      end do

      do i = 1, rank1
         idx1%str(i) = tstr(i:i)
      end do

      do i = 1, rank1
         idx1%itype(i) = itype2(ninter,i)
      end do

      ! Update contraction index posistions
      shift = 1
      do i = 1, rank1
         do j = 1, rank2
            if (idx1%str(i) == idx2%str(j)) then
               idx1%cnt_poss(shift) = i
               shift = shift + 1
            end if
         end do
      end do



!      same = .true.
!      do i = 1, rank1
!         ! Check if intermediate already has correct index itype
!         if (idx1%itype(i) /= itype(i)) then
!            same = .false.
!         end if
!      end do
!
!      if (same) return
!
!
!      do i = 1, rank1/2
!         do j = 1, rank1/2
!            if (itype(i) == idx1%itype(j)) then
!               tstr(i:i) = idx1%str(j)
!               idx1%str(j) = ''
!               idx1%itype(j) = 0
!               exit
!            end if
!         end do
!      end do
!
!      do i = rank1/2 + 1, rank1
!         do j = rank1/2+1, rank1
!            if (itype(i) == idx1%itype(j)) then
!               tstr(i:i) = idx1%str(j)
!               idx1%str(j) = ''
!               idx1%itype(j) = 0
!               exit
!            end if
!         end do
!      end do
!
!      do i = 1, rank1
!         idx1%str(i) = tstr(i:i)
!      end do
!
!      do i = 1, rank1
!         idx1%itype(i) = itype(i)
!      end do
!
!      ! Update contraction index posistions
!      shift = 1
!      do i = 1, rank1
!         do j = 1, rank2
!            if (idx1%str(i) == idx2%str(j)) then
!               idx1%cnt_poss(shift) = i
!               shift = shift + 1
!            end if
!         end do
!      end do

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
            do i = 1, rank/2
!               if (list%plist(i)%pindex(1) >
!     &                                     list%plist(i)%pindex(2)) then
               if (list%plist(i)%nval(1) >=
     &                                     list%plist(i)%nval(2)) then

               if (list%plist(i)%nval(1) ==
     &                                     list%plist(i)%nval(2)) then
               if (list%plist(i)%pindex(1) >
     &                                     list%plist(i)%pindex(2)) then
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
      subroutine assign_spin(item, ntest)
*----------------------------------------------------------------------*
!     Assign spin to tensors, then sum remaining contraction indices
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &   item
      integer, intent(in) ::
     &   ntest

      integer ::
     &   i,j,k,l,m,n,z,            ! Loop index
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
               if (item%abba_line) then
                  ! abba
                  item%t_spin(3)%spin(1,1) = 1
                  item%t_spin(3)%spin(1,2) = 2
                  item%t_spin(3)%spin(2,1) = 2
                  item%t_spin(3)%spin(2,2) = 1
               else
                  ! abab
                  item%t_spin(3)%spin(1,1) = 1
                  item%t_spin(3)%spin(1,2) = 2
                  item%t_spin(3)%spin(2,1) = 1
                  item%t_spin(3)%spin(2,2) = 2
               end if
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
         if (ntest>=100) then
         call debug_header("assign_spin", item%out)
         call print_spin(item%t_spin(3)%spin,item%rank3,item%label_res,
     &                   item%out)
         call print_spin(item%t_spin(1)%spin,item%rank3,item%label_t1,
     &                   item%out)
         call print_spin(item%t_spin(2)%spin,item%rank3,item%label_t2,
     &                   item%out)
         end if

         call spin_index(item, ntest)

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

                  call assign_tensor_spin(item%idx3, item%idx1,
     &                 hr3, hr1, item%t_spin(3),
     &                 item%t_spin(1))
                  call assign_tensor_spin(item%idx3, item%idx2,
     &                 hr3, hr2, item%t_spin(3),
     &                 item%t_spin(2))

                  call spin_index(item, ntest)

               else
                  do k = 1, 2
                     item%t_spin(3)%spin(1,2) = k
                     do l = 1, 2
                        item%t_spin(3)%spin(2,2) = l
                        if (item%rank3==4) then
                           ! Assign spin of external indicies to T1 and T2
                           item%t_spin(1)%spin = 0
                           item%t_spin(2)%spin = 0

                           call assign_tensor_spin(item%idx3, item%idx1,
     &                          hr3, hr1, item%t_spin(3),
     &                          item%t_spin(1))
                           call assign_tensor_spin(item%idx3, item%idx2,
     &                          hr3, hr2, item%t_spin(3),
     &                          item%t_spin(2))

                           call spin_index(item, ntest)
                        else
                           do i = 1, 2
                              item%t_spin(3)%spin(1,3) = i
                              do j = 1, 2
                                 item%t_spin(3)%spin(2,3) = j
                                 ! Assign spin of external indicies to T1 and T2
                                 item%t_spin(1)%spin = 0
                                 item%t_spin(2)%spin = 0

                                 call assign_tensor_spin(item%idx3,
     &                               item%idx1,hr3, hr1, item%t_spin(3),
     &                               item%t_spin(1))
                                 call assign_tensor_spin(item%idx3,
     &                               item%idx2,hr3, hr2, item%t_spin(3),
     &                                item%t_spin(2))

                                 call spin_index(item, ntest)
                              end do
                           end do

                        end if
                     end do
                  end do
               end if
            end do
         end do


      else if (item%rank3==0) then
         call spin_index(item, ntest)
      end if

      end if

      if (ntest>=100) then
         call print_spin(item%t_spin(3)%spin,item%rank3,"Result",
     &                   item%out)
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine assign_tensor_spin(r_idx, t_idx, r_hrank, t_hrank,
     &                              rs, ts)
*----------------------------------------------------------------------*
!     Assign spin from a result tensor to a tensor tensor in the correct
!     possition
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=INDEX_LEN), intent(in) ::
     &   r_idx,                         ! Result index
     &   t_idx                          ! Tensor index
      integer, intent(in) ::
     &   r_hrank,                       ! Result half rank
     &   t_hrank                        ! Tensor half rank
      type(spin_info), intent(in) ::
     &   rs                             ! Result spin
      type(spin_info), intent(inout) ::
     &   ts                             ! Tensor spin

      integer ::
     &   i, j

      do j=1, r_hrank
         do i=1, t_hrank

            if (r_idx(j:j) == t_idx(i:i)) then
               ts%spin(1,i) = rs%spin(1,j)
            else if (r_idx(j:j)== t_idx(i+t_hrank:i+t_hrank)) then
               ts%spin(2,i) = rs%spin(1,j)
            end if

            if (r_idx(j+r_hrank:j+r_hrank) == t_idx(i:i)) then
               ts%spin(1,i) = rs%spin(2,j)
            else if (r_idx(j+r_hrank:j+r_hrank) ==
     &                                   t_idx(i+t_hrank:i+t_hrank))then
               ts%spin(2,i) = rs%spin(2,j)
            end if

         end do
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine spin_index(item, ntest)
*----------------------------------------------------------------------*
!     Find contraction index used in spin summation
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &   item
      integer, intent(in) ::
     &   ntest

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

      if (ntest>=100) call debug_header("spin_index", item%out)

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

      if (item%contri>0)  then
         allocate(poss(2,item%contri))


      if (ntest>=100) then
         call print_spin(item%t_spin(1)%spin,item%rank1,"T1",item%out)
         call print_spin(item%t_spin(2)%spin,item%rank2,"T2",item%out)
      end if

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

      if (ntest>=100) then
         if (item%contri>0) then
            write(item%out,*) "spin z1: ", size(item%t_spin(z1)%spin,2)
            write(item%out,*) "spin z2: ", size(item%t_spin(z2)%spin,2)
            write(item%out,*) "CONTRI: ", item%contri
            write(item%out,*) "SHIFT: ", shift
            write(item%out,*) "========================"
            do i = 1, 2
               do j = 1, item%contri
                  write(item%out,*) "poss: ", poss(i,j)%elements
               end do
            end do
            write(item%out,*) "========================"
         end if
      end if

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
            call save_spin_case(item,eloop,ntest)
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
                  call save_spin_case(item,eloop,ntest)
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
                        call save_spin_case(item,eloop,ntest)
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
                              call save_spin_case(item,eloop,ntest)
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
                                    call
     &                                  save_spin_case(item,eloop,ntest)
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
                                       call
     &                                  save_spin_case(item,eloop,ntest)
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
         if (ntest>=100) then
            write(item%out,*) "No contraction indices: ", item%contri
         end if

         !if (r1>0) item%t_spin(z1)%spin = item%t_spin(3)%spin
         !if (r2>0) item%t_spin(z2)%spin = item%t_spin(3)%spin

         ! TODO: factorise this, maybe combine with stuff above?
         if (item%rank1>0) then
            do i = 1, item%rank3/2
               do j = 1, item%rank1/2
                  if (item%idx3(i:i)==item%idx1(j:j)) then
                     item%t_spin(1)%spin(1,j) = item%t_spin(3)%spin(1,i)
                  end if
               end do
            end do
            do i = item%rank3/2+1, item%rank3
               do j = item%rank1/2+1, item%rank1
                  if (item%idx3(i:i)==item%idx1(j:j)) then
                     item%t_spin(1)%spin(2,j-item%rank1/2) =
     &                             item%t_spin(3)%spin(2,i-item%rank3/2)
                  end if
               end do
            end do
         end if

         if (item%rank2>0) then
            do i = 1, item%rank3/2
               do j = 1, item%rank2/2
                  if (item%idx3(i:i)==item%idx2(j:j)) then
                     item%t_spin(2)%spin(1,j) = item%t_spin(3)%spin(1,i)
                  end if
               end do
            end do
            do i = item%rank3/2+1, item%rank3
               do j = item%rank2/2+1, item%rank2
                  if (item%idx3(i:i)==item%idx2(j:j)) then
                     item%t_spin(2)%spin(2,j-item%rank2/2) =
     &                             item%t_spin(3)%spin(2,i-item%rank3/2)
                  end if
               end do
            end do
         end if


         call save_spin_case(item,eloop,ntest)
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
      subroutine save_spin_case(item,eloop,ntest)
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
      integer, intent(in) ::
     &   ntest

      integer ::
     &   i
      integer ::
     &   sum1a, sum1b, sum2a, sum2b

      if (ntest>=100) call debug_header("save_spin_case", item%out)

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


            item%all_spins(item%nspin_cases)%t_spin(1)%spin=
     &                                          item%t_spin(1)%spin
            item%all_spins(item%nspin_cases)%t_spin(2)%spin=
     &                                          item%t_spin(2)%spin
            item%all_spins(item%nspin_cases)%t_spin(3)%spin=
     &                                          item%t_spin(3)%spin
            item%nspin_cases = item%nspin_cases + 1


            if (ntest>=100) then
               write(item%out,*) "Saving spin case to item%t_spin:"
               call print_spin(item%t_spin(3)%spin,item%rank3,"Spin 3",
     &                         item%out)
               call print_spin(item%t_spin(1)%spin,item%rank1,"Spin 1",
     &                         item%out)
               call print_spin(item%t_spin(2)%spin,item%rank2,"Spin 2",
     &                         item%out)
               write(item%out,*) "Number of spin cases thus far: ",
     &                           item%nspin_cases - 1

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
      subroutine print_itf_contr(item)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &     item     ! Object which holds information necessary to print out an ITF algo line

      integer ::
     &   i,j,k,
     &   l1,l2,l3

      write(item%out,*)
      write(item%out,*) "================================="
      write(item%out,*) trim(item%label_res)//"["//
     & trim(item%idx3)//"]"//
     & " = "//trim(item%label_t1)//"["//
     & trim(item%idx1)//"]"//" "//
     & trim(item%label_t2)//"["//
     & trim(item%idx2)//"]"
      write(item%out,*) "FACTOR: ", item%fact
      write(item%out,*) "SPIN CASES: ", item%spin_cases-1

      l1 = item%rank1/2
      l2 = item%rank2/2
      l3 = item%rank3/2

      do i = 1, item%spin_cases-1
      write(item%out,*)
      write(item%out,*) "---------------------------------"
      do k = 2, 1, -1
      write(item%out,'(2a)',advance='no') "  "
      do j = 1, l3
      write(item%out,'(i1)',advance='no')
     &     item%all_spins(i)%t_spin(3)%spin(k,j)
      end do
      write(item%out,'(2a)',advance='no') "  "
      do j = 1, l1
      write(item%out,'(i1)',advance='no')
     &     item%all_spins(i)%t_spin(1)%spin(k,j)
      end do
      write(item%out,'(2a)',advance='no') "  "
      do j = 1, l2
      write(item%out,'(i1)',advance='no')
     &     item%all_spins(i)%t_spin(2)%spin(k,j)
      end do
      write(item%out,*)
      end do
      end do

      write(item%out,*) "================================="
      write(item%out,*)

      return
      end

*----------------------------------------------------------------------*
      subroutine print_symmetrise(item, ntest)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &   item
      integer, intent(in) ::
     &   ntest

      character(len=70) ::
     &   line
      character(len=INDEX_LEN) ::
     &   tindex
      character(len=MAXLEN_BC_LABEL) ::
     &   new

      ! Return if an intermediate or result rank less than 2
      if (item%inter(3) .or. item%rank3<=2) return

      ! Return if not a symmetric result
      if (.not. item%symmetric) return

      ! Return if INTPP block
      if (item%intpp) return

      if (ntest>=100) call debug_header("print_symmetrise", item%out)

      new = rename_tensor(item%old_name, item%rank3)

      line = '.'//trim(new)//'['//trim(item%idx3)//'] += '//
     &       trim(item%label_res)//'['//trim(item%idx3)//']'
      write(item%out,'(a)') trim(line)

      tindex = ' '
      tindex(1:1) = item%idx3(2:2)
      tindex(2:2) = item%idx3(1:1)
      tindex(3:3) = item%idx3(4:4)
      tindex(4:4) = item%idx3(3:3)

      line = '.'//trim(new)//'['//trim(item%idx3)//'] += '//
     &       trim(item%label_res)//'['//trimal(tindex)//']'
      write(item%out,'(a)') trim(line)

      return
      end


*----------------------------------------------------------------------*
      subroutine check_symmetric(contr_info, command, symmetric, nosym)
*----------------------------------------------------------------------*
!     Check if a tensor has permutational symmetry
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &   contr_info     ! Information about binary contraction
      integer, intent(in) ::
     &   command
      logical, intent(inout) ::
     &   symmetric,  ! Check if tensor has R[abij] = R[baji] symmetry
     &   nosym       ! Check if tensor is R[apiq]

      integer ::
     &   c(ngastp,2),
     &   e1(ngastp,2),
     &   e2(ngastp,2),
     &   nops(ngastp),
     &   i

      call itf_ops(contr_info, c, e1, e2, command)

      nops = sum(e1, dim=2) + sum(e2, dim=2)

      if (sum(nops)==0) then
         ! Scalars don't have symmetry
         symmetric = .false.
         nosym = .false.
         return
      end if

      symmetric = .true.
      do i = 1, ngastp
         if (mod(nops(i),2) /= 0) then
            symmetric = .false.
         end if
      end do

      nosym = .false.
      if (nops(1)==1 .and. nops(2)==1 .and. nops(3)==2) then
         if (.not. check_inter(contr_info%label_res)) then
            nosym = .true.
         end if
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine check_k4e(contr_info, command, k4e_line)
*----------------------------------------------------------------------*
!     Check if a tensor has permutational symmetry
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &   contr_info     ! Information about binary contraction
      integer, intent(in) ::
     &   command
      logical, intent(inout) ::
     &   k4e_line

      integer ::
     &   c(ngastp,2),
     &   e1(ngastp,2),
     &   e2(ngastp,2),
     &   nops(ngastp),
     &   i

      call itf_ops(contr_info, c, e1, e2, command)

      nops = sum(e1, dim=2) + sum(c, dim=2)

      k4e_line = .false.

      if (sum(nops)==0 .or. sum(nops)==2) then
         return
      end if

      if (nops(2)==4) then
         k4e_line = .true.
      else
         nops = sum(e2, dim=2) + sum(c, dim=2)
         if (nops(2)==4) then
            k4e_line = .true.
         end if
      end if

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
           call count_index(contr_info%occ_op1(1:,1:,i), e1)
         end do
      else
         ! Get occupation info
         do i = 1, contr_info%n_cnt
           call count_index(contr_info%occ_cnt(1:,1:,i), c)
         end do
         do i = 1, contr_info%nj_op1
           call count_index(contr_info%occ_ex1(1:,1:,i), e1)
         end do
         do i = 1, contr_info%nj_op2
           call count_index(contr_info%occ_ex2(1:,1:,i), e2)
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
