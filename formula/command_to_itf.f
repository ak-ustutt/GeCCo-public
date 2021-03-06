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
*
c     * moved to extra file as the compiler incorrectly quits with
c     * an error message (error #6633)
c*----------------------------------------------------------------------*
c      pure function rename_tensor(string, rank, nops, itf_names)
c*----------------------------------------------------------------------*
c!     Rename tensor according to taste
c!     This should be expanded to rename all tensors so they correspond
c!     with the ITF algo file
c*----------------------------------------------------------------------*
c
c      implicit none
c      include 'stdunit.h'
c      include 'opdim.h'
c      include 'def_contraction.h'
c      include 'def_itf_contr.h'
c
c      character(len=MAXLEN_BC_LABEL), intent(in) ::
c     &    string
c      integer, intent(in) ::
c     &    rank,          ! Rank of tensor
c     &    nops(ngastp)   ! Number of creation/annhilation ops per excitation space
c      type(tensor_names), intent(in) ::
c     &     itf_names            ! contains renaming information
c      character(len=MAXLEN_BC_LABEL) ::
c     &     rename_tensor
c
c      logical ::
c     &     found
c      integer ::
c     &     ii
c
c      found = .false.
c      do ii = 1, itf_names%nrename
c        if (trim(string)==trim(itf_names%gc_itf_rename(1,ii))) then
c          found = .true.
c          select case(itf_names%rename_type(ii))
c          case(RENAME_BASIC)
c            rename_tensor = trim(itf_names%gc_itf_rename(2,ii))
c          case(RENAME_HAMILTONIAN)
c            if (rank==2) then
c              rename_tensor='f'
c            else
c              if (nops(2)==3) then
c                rename_tensor='J' ! Currently use J:eeec for 3 external integrals
c              else
c                rename_tensor='K'
c              end if
c            end if
c          case(RENAME_ADD_RANK)
c            if (rank>18) then
c              !write(lulog,*) 'rank = ',rank
c              rename_tensor='ERROR_unexp_lg_rank'
c            end if
c            write(rename_tensor,'(a,i1)')
c     &           trim(itf_names%gc_itf_rename(2,ii)),rank/2
c          case default
c            !write(lulog,*) 'case ',itf_names%rename_type
c            rename_tensor='ERROR_unknown_case'
c          end select
c          exit
c        end if
c      end do
c
c      if (.not.found) rename_tensor = trim(string)
c
c      if (rename_tensor(1:1)=='_') rename_tensor(1:1)=''
c
c      end function

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
      if (index(label, "STIN")>0 .or.  ! <--- this is a bit dangerous:
     &    index(label, "LTIN")>0 .or.  !      the user might use a variable like 'FASTIN'
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
cc      if (index(label, "H")>0) then
      if (trim(label)=="H") then
         check_int=.true.
      else
         check_int=.false.
      end if

      end function


!*----------------------------------------------------------------------*
!      pure function check_den(label)
!*----------------------------------------------------------------------*
!!     Check if tensor is a density matrix
!*----------------------------------------------------------------------*
!
!      implicit none
!      include 'opdim.h'
!      include 'def_contraction.h'
!      include 'def_itf_contr.h'
!
!      character(len=MAXLEN_BC_LABEL), intent(in) ::
!     &     label
!
!      logical ::
!     &     check_den
!
!      ! Assume these are the names of intermediates
!      if (index(label, "GAM")>0) then
!         check_den=.true.
!      else
!         check_den=.false.
!      end if
!
!      end function


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
      else if (scan("pqrstuvwxyz",idx)>0) then
         itype = vals(3)
      !else if (scan("xyz",idx)>0) then
      !   itype = vals(4)
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
      write(logfile,'(1x,10i4)') (spins(2,i), i=1, rank/2)
      write(logfile,'(1x,10i4)') (spins(1,i), i=1, rank/2)
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

! for debugging, it is better to exit here
      call quit(1,'ITF','Trouble! Check bcontr.tmp')

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
      subroutine set_itype(item, itype, ntest)
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
     &   itype(MAXINT, INDEX_LEN)
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

      ! More advanced itype for handling n intermeidates which are
      ! defined and used anywhere within the complete diagram
      ! WARNING: This won't work for more than 9 intermediates...
      i = len(trim(item%label_res))
      read(item%label_res(i:i),*) ninter
      if (ninter==0)
     &          call line_error("0 ninter in set_itype",item)

      if (ninter.gt.MAXINT) then
        write(item%out,'(1x,"MAXINT = ",i5,"  ninter=",i10)')
     &       MAXINT, ninter
        call quit(1,'set_itype','dimension error')
      end if

      ! Clear all itype info from previous diagram
      if (i == 1) then
         itype = 0
      end if

      if (item%rank3==0) then
        ! nothing to do
      else if (item%rank3==2) then
         itype(ninter,1) = get_itype(item%idx3(1:1))
         itype(ninter,2) = get_itype(item%idx3(2:2))
      else if (item%rank3==4) then
         itype(ninter,1) = get_itype(item%idx3(1:1))
         itype(ninter,2) = get_itype(item%idx3(2:2))
         itype(ninter,3) = get_itype(item%idx3(4:4))
         itype(ninter,4) = get_itype(item%idx3(3:3))
      else if (item%rank3==6) then
         itype(ninter,1) = get_itype(item%idx3(1:1))
         itype(ninter,2) = get_itype(item%idx3(2:2))
         itype(ninter,3) = get_itype(item%idx3(3:3))
         itype(ninter,4) = get_itype(item%idx3(6:6))
         itype(ninter,5) = get_itype(item%idx3(5:5))
         itype(ninter,6) = get_itype(item%idx3(4:4))
      else if (item%rank3==8) then
         itype(ninter,1) = get_itype(item%idx3(1:1))
         itype(ninter,2) = get_itype(item%idx3(2:2))
         itype(ninter,3) = get_itype(item%idx3(3:3))
         itype(ninter,4) = get_itype(item%idx3(4:4))
         itype(ninter,5) = get_itype(item%idx3(8:8))
         itype(ninter,6) = get_itype(item%idx3(7:7))
         itype(ninter,7) = get_itype(item%idx3(6:6))
         itype(ninter,8) = get_itype(item%idx3(5:5))
      else if (item%rank3==10) then
         itype(ninter,1) = get_itype(item%idx3(1:1))
         itype(ninter,2) = get_itype(item%idx3(2:2))
         itype(ninter,3) = get_itype(item%idx3(3:3))
         itype(ninter,4) = get_itype(item%idx3(4:4))
         itype(ninter,5) = get_itype(item%idx3(5:5))
         itype(ninter,6) = get_itype(item%idx3(10:10))
         itype(ninter,7) = get_itype(item%idx3(9:9))
         itype(ninter,8) = get_itype(item%idx3(8:8))
         itype(ninter,9) = get_itype(item%idx3(7:7))
         itype(ninter,10)= get_itype(item%idx3(6:6))
       else
         write(item%out,*) 'rank3 = ',item%rank3
         write(item%out,*) 'extend me for missing intermediate rank'
         call quit(1,'set_itype',
     &       'extend me for missing intermediate rank')
      end if


      if (ntest>=100) then
         call debug_header("set_itype", item%out)
         write(item%out,'(1x,"Have just registered: idx3= ",a)')
     &        item%idx3
         do j = 1, MAXINT
            write(item%out,'(a6)',advance='no') "     {"
            do i = 1, item%rank3
               write(item%out,'(i0, 1x)', advance='no') itype(j, i)
            end do
            write(item%out,'(a1)') "}"
         end do
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine command_to_itf(contr_info, itin, itflog, command,
     &                          inter_itype,
     &                          itf_names,
     &                          counter, tasks, x_dict,
     &                          inter_spin_dict)
*----------------------------------------------------------------------*
!     Take GeCco binary contraction and produce ITF algo code.
!     Includes antisymmetry of residual equations and spin summation.
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(inout) ::
     &     contr_info           ! Information about binary contraction
      type(tensor_names), intent(in)   ::
     &     itf_names            ! contains renaming information
      logical, intent(in) ::
     &   itin,               ! Print ITIN lines or not
     &   tasks
      integer, intent(in) ::
     &   itflog,             ! Output file
     &   command             ! Type of formula item command, ie. contraction, copy etc.
      integer, intent(inout) ::
     &   counter(4)
      integer, intent(inout) ::
     &   inter_itype(MAXINT, INDEX_LEN)      ! Store itypes of intermediates between lines
      type(x_inter), intent(inout) ::
     &   x_dict(MAXX)
      type(inter_spin_cases), intent(inout) ::
     &   inter_spin_dict     ! Store intermediates spin cases between lines

      type(itf_contr) ::
     &   item,               ! ITF contraction object; holds all info about the ITF algo line
     &   pitem               ! Permutation ITF contraction object
      integer ::
     &   i, j,                 ! Loop index
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

      logical, parameter :: new = .true.

      if (.not.new) call warn('command_to_itf','running in old mode!')
      
      if (ntest1>=100 .or. ntest2>=100) then
         write(lulog,*) "DEBUG MODE: look for output in file bcontr.tmp"
         write(itflog,*) "FORMULA NUMBER: ", counter(1)
         call prt_bcontr(itflog,contr_info)
      end if

      ! Begin a special block which the python processor will pull out
      ! into its own code block
      intpp = .false.
      !if (contr_info%label_res=='INTpp') thenend
      !   intpp = .true.
      !   write(itflog,'(a)') "BEGIN_INTPP"
      !end if

      ! Mark begining of spin summed block
      if (.not. tasks) write(itflog,'(a5)') 'BEGIN'

      ! 1. Initalise itf_contr
      call itf_contr_init(contr_info,item,1,itin,command,itflog,
     &                    inter_itype,counter,tasks,new,ntest1)


      ! 2. Assign index / Determine sign
      call assign_index(item,contr_info,ntest2)


      ! 3. Determine if we need a permutation line and create new item
      !    for it
      !    Permutation line required for:
      !        a) Symmetric residuals with two pairs of external indicies
      !        a) Non-symmetric residuals with one pair of external indicies
      pline = .false.
      perm_case = 0
      do i = 1, ngastp
        do j = 1, 2
          if(contr_info%perm(i,j)) perm_case = perm_case + 1
        end do
      end do

      if (ntest3>=100) then
         call debug_header("Determine permutation", item%out)
         write(item%out,'(1x,"contr_info%perm was ",4l4,2x,4l4)')
     &        contr_info%perm(1:ngastp,2)        
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
c         if (.not. item%inter(3) .and. .not. item%product) then  ! <--- why not for tensor product?
         if (.not. item%inter(3) ) then  ! <--- why not for tensor product?
            pline = .true.

            call itf_contr_init(contr_info,pitem,1,itin,command,itflog,
     &                           inter_itype,counter,tasks,new,ntest5)
            call assign_index(pitem,contr_info,ntest5)

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
     &                       inter_itype,counter,tasks,new,ntest5)
         call assign_index(pitem,contr_info,ntest5)

         pitem%old_name = pitem%label_res
         pitem%label_res = item%label_res
         pitem%abba_line = .true.
      end if


      ! 6. Spin sum
      call assign_spin(item, ntest6)
      if (pline) call assign_spin(pitem, ntest6)


      ! 7. Loop over spin cases and print out each line
      call print_spin_cases(item, x_dict, inter_spin_dict, itf_names,
     &                      ntest7)
      if (pline) call print_spin_cases(pitem, x_dict, inter_spin_dict,
     &                                 itf_names, ntest7)


      ! 8. Print symmetrisation term
      if (.not. item%inter(3)) then
         if (.not. tasks) call print_symmetrise(item, itf_names, ntest8)
      end if


      ! 9. If an intermediate, set inter_itype for use in next line
      ! where the intermediate is created

      call set_itype(item, inter_itype, ntest9)


      if (ntest10>=100) then
         call print_itf_contr(item)
         if (pline) call print_itf_contr(pitem)
      end if


      ! Mark end of spin block
      if (.not. tasks) write(itflog,'(a)') "END"

      if (intpp) then
         write(itflog,'(a)') "END_INTPP"
      end if


      ! Update counter of quantaties
      ! Number of intermediates
      if (.not. item%inter(3)) then
         counter(2) = item%cntr(2) + 1
      else
         counter(2) = item%cntr(2)
      end if

      ! Update K4E counter
      if (item%k4e_line) then
         counter(3) = item%cntr(3) + 1
      end if

      ! Update X number
      counter(4) = item%cntr(4)

      ! Reset inter_spin_dict
      if (.not. item%inter(3)) then
         inter_spin_dict%ncase = 0
      end if

      ! Deallocate memroy used when construcitng item
      call itf_deinit(item)

      return
      end


*----------------------------------------------------------------------*
      subroutine itf_contr_init(contr_info,item,perm,itin,comm,loclog,
     &                          itype,counter,tasks,use_sign,ntest)
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
     &   loclog,        ! Output file
     &   itype(MAXINT,INDEX_LEN),
     &   counter(4),
     &   ntest
      logical, intent(in) ::
     &   itin,
     &   tasks,
     &   use_sign

      integer ::
     &   i
      character(len=60) ::
     &     tmp

      if (ntest>=100) then
        call debug_header("itf_contr_init", loclog)
      end if

      ! Assign output file
      item%out=loclog

      ! Set tasks file
      item%tasks = tasks

      ! Assign command type
      item%command=comm

      ! Initalise index strings
      item%idx1 = ''
      item%idx2 = ''
      item%idx3 = ''

      item%binary = .true.
      if (item%command==command_cp_intm .or.
     &     item%command==command_add_intm .or.
     &     item%command==command_add_reo ) then
         ! For [ADD] and [COPY]
         ! Not a binary contraction
         item%binary = .false.
      end if
      
      ! Get number of contraction and external indices on each tensor
      call itf_ops(contr_info, item%c, item%e1, item%e2, item%binary)

      ! Set ranks of tensors using matricies from itf_ops
      call itf_rank(item%e1, item%c, item%rank1, item%nops1, .false.)
      call itf_rank(item%e2, item%c, item%rank2, item%nops2, .false.)  ! should be fine for non-binary cmd.
      call itf_rank(item%e1, item%e2, item%rank3, item%nops3, .false.)

      if (ntest>=100) then
        write(item%out,'(1x,"determined rankX: ",3i6)')
     &       item%rank1, item%rank2, item%rank3
        write(item%out,'(1x,"determined nopsX: ",3i6)')
     &       item%nops1, item%nops2, item%nops3        
      end if
      
      ! Set number of contraction indicies
      item%contri = sum(sum(item%c, dim=1))

      ! Set external lines (used for permuations)
      ! UNUSED item%perm_case = contr_info%perm

      ! Determine factor from equivalent lines
      item%fact = 1.0d+0
      !! is now contained in fact_itf
      !!call itf_equiv_lines_factor(item%c, item%fact)

      if (ntest>=100 .and. abs(item%fact-1d0).gt.1d-12)
     &    write(item%out,'(1x,"eqv. lines factor was",f20.12)')item%fact

      ! Get any remaining factors from GeCCo
      if (use_sign) item%fact = item%fact * contr_info%fact_itf
      ! note: fact_itf contains the general contraction prefactor as defined in GeCCo and a possible sign
      ! change according to the sign conventions employed for ITF code generation

      ! this is old (to be deleted soon):
      if (.not.use_sign) item%fact = item%fact *abs(contr_info%fact_itf)
      !item%fact = item%fact * contr_info%fact_itf

      if (ntest>=100)
     &     write(item%out,'(1x,"factor set to: ",f12.6)') item%fact

      ! Assign labels
      item%label_t1=contr_info%label_op1

      ! Check if an intermediate
      item%inter(1) = check_inter(item%label_t1)
      if (.not. item%inter(1)) then
         if (.not. item%den(1)) then
            ! Check if an integral
            item%int(1) = check_int(item%label_t1)
         end if
      end if

      ! Operator 2 does not exist in [ADD] or [COPY]
      if (item%binary) then

         item%label_t2=contr_info%label_op2

         item%inter(2) = check_inter(item%label_t2)

         if (.not. item%inter(2)) then
            if (.not. item%den(2)) then
               item%int(2) = check_int(item%label_t2)
            end if
         end if
      else
         item%label_t2="" 
      end if

      ! Assign permutation number
      item%permute=perm

      item%label_res=contr_info%label_res
      item%inter(3) = check_inter(item%label_res)
      item%int(3) = .false.

      item%symm_res = .false. ! do not leave it undefined although never used
      call check_symmetric(contr_info, item%binary, item%symmetric,
     &                     item%nosym, item%out)

      ! If a residual, is it symmetric (R_{ab}^{ij] = R_{ba}^{ji})?
      if (.not. item%inter(3)) then
         call check_symmetric(contr_info, item%binary, item%symm_res,
     &                        item%nosym, item%out)

         ! Add factor to account for factor of two when the residual is
         ! symmetrised:
         ! R:eecc[abij] += G:eecc[abij]
         ! R:eecc[abij] += G:eecc[baji]
         if (item%symm_res .and. item%permute==0) then
            if (contr_info%label_res=='INTpp')
     &          call quit(1,'command_to_itf','does this happen?')
            ! dirty explicit reference
            if (contr_info%label_res/='INTpp') then
               write(item%out,*) "0.5 for ",trim(contr_info%label_res) 
               item%fact = item%fact * 0.5d+0
            end if
         end if
      end if


      ! Remove zeorth-body density from equations
      ! will be redundant soon; direct reference to name is problematic
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

      ! OLD: see comment beneath

      ! OLD: Assign factor --- use special ITF factor
      ! OLD: the ITF factor is closer to the value expected from standard
      ! OLD: diagram rules, but still some care has to be taken when translating
      ! OLD: to ITF tensors; e.g. the (HP;HP) integrals are stored by GeCCo as
      ! OLD: <aj||bi> while the standard sign for ring terms assumes <aj||ib>
      ! OLD -- commented out code:item%fact=contr_info%fact_itf

      ! All this is now set a few lines before ... the sign is now determined
      ! exactly for the chosen index ordering, pairwise shifts of corresponding
      ! indices will always be OK

      ! Assign index string. Tensor ranks and number
      ! of operators are also set here

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


      ! Check if a tensor product (i.e. no contraction)
      if (item%rank3==4 .and. item%rank1==2 .and. item%rank2==2) then
         item%product=.true.
      end if

      ! Number of spin cases associated with this line
      item%nspin_cases=1

      ! Set itype from previous line (can be 0)
      item%itype = itype

      ! Set K4E info
      item%k4e_line = .false. ! do not leave it undefined!
      if (item%int(1) .or. item%int(2))then
         call check_k4e(contr_info, item%command, item%k4e_line)
      end if


      ! Set vertex info, used to calculate overal sign in assign_index()
      allocate(item%vertex(contr_info%nj_op1 + contr_info%nj_op2))
      do i = 1, contr_info%nj_op1 + contr_info%nj_op2
         item%vertex(i) = -1 ! OBSOLETE contr_info%svertex_itf(i)
      end do

      !if (ntest>=0) write(item%out,*) "vertex ", item%vertex

      item%nj_op1 = contr_info%nj_op1
      item%nj_op2 = contr_info%nj_op2
      item%nj_res = contr_info%nj_res

      call  itf_vertex_ops(contr_info, item, comm)


      ! Set counter information that is outside of the loop
      item%cntr(1) = counter(1) ! Contraction number
      item%cntr(2) = counter(2) ! Intermediate number
      item%cntr(3) = counter(3) ! K4E number
      item%cntr(4) = counter(4) ! X intermediate number

      if (tasks) then
         ! Need to number intermediates in the task file
         write(tmp,*) item%cntr(2)

         if (item%inter(1)) then
            item%label_t1 = "STIN"//trim(adjustl(tmp))//
     &                       item%label_t1(6:9)
         end if

         if (item%inter(2)) then
            item%label_t2 = "STIN"//trim(adjustl(tmp))//
     &                       item%label_t2(6:9)
         end if

         if (item%inter(3)) then
            item%label_res = "STIN"//trim(adjustl(tmp))//
     &                       item%label_res(6:9)
         end if
      end if


      if (ntest>=100) then
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

      deallocate(item%vertex)

      ! TODO: probably don't need to check if associated
      if (associated(item%v1)) deallocate(item%v1)
      if (associated(item%v2)) deallocate(item%v2)
      if (associated(item%v3)) deallocate(item%v3)
      if (associated(item%vc1)) deallocate(item%vc1)
      if (associated(item%vc2)) deallocate(item%vc2)
      if (associated(item%vc3)) deallocate(item%vc3)

      deallocate(item%vnops1)
      deallocate(item%vnops2)
      deallocate(item%vnops3)

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
     &   perm(ngastp,2)
      integer, intent(in) ::
     &   ntest

      integer ::
     &   ex_itype,
     &   i, j,
     &   shift, rank_off
      character(len=1) ::
     &   ex_ind(2)
      character(len=INDEX_LEN) ::
     &   tmp1, tmp2

      if (ntest>=100) then
         call debug_header("create_permuation", item%out)
         write(item%out,'(1x,"on entry: idx1 = ",a)') item%idx1
         write(item%out,'(1x,"on entry: idx2 = ",a)') item%idx2
         write(item%out,'(1x,"on entry: idx3 = ",a)') item%idx3
      end if

      ! be save
      if (item%rank3/=4) call quit(1,'create permutation',
     &                   'I only work for result rank 4')

      tmp1 = item%idx1
      tmp2 = item%idx2

      ! Identify and swap external indicies
      do j = 1, 2
        do i = 1, ngastp
          if (perm(i,j)) then
            ex_itype = i
            exit
          end if
        end do
      end do

      ! check whether C or A should by symmetrised
      if (get_itype(item%idx3(1:1),.true.)==ex_itype.and.
     &    get_itype(item%idx3(2:2),.true.)==ex_itype) then
        rank_off = 0
      else if (get_itype(item%idx3(3:3),.true.)==ex_itype.and.
     &         get_itype(item%idx3(4:4),.true.)==ex_itype) then
        rank_off = 2
      else
        write(item%out,*) 'trouble in create_permutation:'
        write(item%out,*) 'idx3 = ',item%idx3,' ex_itype = ',ex_itype
        call quit(1,'create_permutation','I''am in trouble')
      end if

      shift = 1
      do i = rank_off+1, rank_off+item%rank3/2
         !if (get_itype(item%idx3(i:i),.true.)/=ex_itype) cycle

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
         write(item%out,*) "after create_permuation:"
         write(item%out,*) "external index: ", ex_ind
         write(item%out,*) "external itype: ", ex_itype
         write(item%out,*) "new factor: ", item%fact
         write(item%out,*) "new idx1: ",item%idx1
         write(item%out,*) "new idx2: ",item%idx2
         write(item%out,*) "    idx3: ",item%idx3
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine convert_to_abab_block(item, t_spin, idx1, idx2,
     &                                 idx3, fact)
*----------------------------------------------------------------------*
!     Convert integrals and amplitudes to abab spin blocks, this may also
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
     &   t_spin(3)
      character(len=INDEX_LEN), intent(inout) ::
     &   idx1, idx2, idx3
      real(8), intent(inout) ::
     &   fact

      character(len=1) ::
     &   tmp
      integer ::
     &   ops(4,2)


      if (item%rank1 > 2 .and. .not. item%inter(1)) then
         call reorder_spin_index(item, t_spin(1)%spin, item%rank1,
     &                           idx1, fact)
      end if

      if (item%rank2 > 2 .and. .not. item%inter(2)) then
         call reorder_spin_index(item, t_spin(2)%spin, item%rank2,
     &                           idx2, fact)
      end if

      ! For the R[apiq] abba case, we need to rearange
      if (item%rank3 > 2 .and. item%abba_line
     &    .and. .not. item%inter(3)) then
         call reorder_spin_index(item, t_spin(3)%spin, item%rank3,
     &                           idx3, fact)
      end if


      return
      end


*----------------------------------------------------------------------*
      subroutine reorder_spin_index(item, spin, rank, idx, fact)
*----------------------------------------------------------------------*
!     Reorder index according to spin.
!     All alpha indicies are moved to the left
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(in) ::
     &   item
      integer, intent(in) ::
     &   spin(2,*)            ! Spin info of a tensor
      integer, intent(in) ::
     &   rank
      character(len=INDEX_LEN), intent(inout) ::
     &   idx
      real(8), intent(inout) ::
     &   fact                 ! Overall factor

      character(len=INDEX_LEN) ::
     &   tmp
      integer ::
     &   a_shift,             ! Position of the next alpha index
     &   s,
     &   i, j, k, l,          ! Loop indices
     &   z


      ! Loop over covarient/contravarient indcies
      do l = 1, 2

         if (l == 1) then
            a_shift = 1
         else
            a_shift = 1 + rank/2
         end if

!         write(item%out,*) "Old index:"
!         do k = 1, rank
!            write(item%out,*) idx(k:k)
!         end do
!         write(item%out,*) "Old factor: ", fact
!         call print_spin(spin, rank, "reorder_spin_index", item%out)

         ! Loop over spin information
         do i = 1, rank/2

            if (l == 1) then
               z = i
            else
               z = i + rank/2
            end if

            ! Find a alpha spin index
            if (spin(l,i)==1) then

               ! Skip if alpha index already in position
               if (z == a_shift) then
                  a_shift = a_shift + 1
                  cycle
               end if

               ! Move index to the left
!               write(item%out,*) "Moving index: ", idx(z:z)
!               write(item%out,*) "Spin: ", spin(l,i)
               tmp(a_shift:a_shift) = idx(z:z)

               ! Calculate factor from permuting indices
               !write(item%out,*) "distance ", z - a_shift
               if (mod(z-a_shift,2)/=0)then
                  fact = fact * -1.0d0
               end if

               ! Arrange renaiming indicies into tmp
               s = 1
               do j = 1, rank
                  if (j == z) cycle
                  if (s == a_shift) then
                     tmp(s+1:s+1) = idx(j:j)
                     s = s + 2
                  else
                     tmp(s:s) = idx(j:j)
                     s = s + 1
                  end if
               end do

               ! Update a_shift for next alpha position
               a_shift = a_shift + 1

               ! Replace idx with tmp
               do k = 1, rank
                  idx(k:k) = tmp(k:k)
               end do

!               write(item%out,*) "New index:"
!               do k = 1, rank
!                  write(item%out,*) idx(k:k)
!               end do
!               write(item%out,*) "New factor: ", fact

            end if
         end do

      end do


      return
      end


*----------------------------------------------------------------------*
      subroutine create_exchange_string(rank,spin,str,nxstr,str_x)
*----------------------------------------------------------------------*
*     create the corresponding exchange strings, depending on spin
*     symmetry:
*     rank = 4, aaaa/bbbb case:  [abij] -> [baij]
*     rank = 6, aaaaaa/bbbbbb case: [abcijk] -> [bcaijk], [cabijk], etc.
*     rank = 6, aaaabb/bbbbaa case: [abcijk] -> [acbijk]
*
*     andreas, dec 2020 based on routines by j.a.black
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      logical, intent(in) :: spin ! true if pure aaaa... or bbbb... case
      integer, intent(in) :: rank ! rank of string
      integer, intent(inout) :: nxstr ! in:  max. dim of str_x
                                      ! out: number of gen. str_x
      character(len=INDEX_LEN), intent(in) :: str  ! input str
      character(len=INDEX_LEN), intent(out) :: str_x(nxstr) ! output str

      if (rank<4) then
        nxstr = 0
        return                  ! nothing to do
      end if
      if (rank>6) call quit(1,'create_exchange_string',
     &     'cannot handle rank > 6 at present!')

      if (rank==4.and.spin) then
        if (nxstr<1) call quit(1,'create_exchange_string',
     &                           'dimension error (rk 4)')
        str_x(1) = f_index(str,2)
        nxstr = 1
        return
      else if (rank==4.and..not.spin) then
        nxstr = 0
        return
      end if

      if (rank == 6) then
         if (spin) then
           if (nxstr<5) call quit(1,'create_exchange_string',
     &          'dimension error (rk 6 a)')
           str_x(1) = c_index(str,1)  ! cab
           str_x(2) = c_index(str,2)  ! bca
           str_x(3) = f_index(str,3)  ! acb
           str_x(4) = f_index(c_index(str,1),3) ! cba
           str_x(5) = f_index(c_index(str,2),3) ! bac
           nxstr = 5
           return
         else
           if (nxstr<1) call quit(1,'create_exchange_string',
     &          'dimension error (rk 6 b)')
           str_x(1) = f_index(str,3)  ! acb
           nxstr = 1
           return
         end if
      end if

      end subroutine


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

c      write(item%out,'(1x,"in prepare_symmetrise: ",2l4)')
c     &     item%symmetric, item%symm_res

      ! Return if INTPP block
      if (item%intpp) return  ! now deprecated

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
      subroutine print_spin_cases(item, x_dict, inter_spin_dict,
     &                            itf_names, ntest)
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
      type(x_inter), intent(inout) ::
     &   x_dict(MAXX)
      type(inter_spin_cases), intent(inout) ::
     &     inter_spin_dict
      type(tensor_names), intent(in) ::
     &     itf_names            ! contains renaming information
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

c      if (actual_spin_cases<=1) write(item%out,*) 'Only 0 spin cases?',
c     &     item%nspin_cases
c      if (actual_spin_cases<=1)
c     &     call warn('ITF:print_spin_cases','less than 2 spin cases')

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
     &                        x_dict, inter_spin_dict, itf_names, ntest)

      end do

      return
      end

*----------------------------------------------------------------------*
      subroutine print_itf_line(item, s1, s2, t_spin, x_dict,
     &                          inter_spin_dict, itf_names, ntest)
*----------------------------------------------------------------------*
!     Print line of ITF code
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'stdunit.h' ! for low-level debug
      include 'opdim.h'
      include 'mdef_operator_info.h' ! For def_formular_item.h
      include 'def_contraction.h'
      include 'def_formula_item.h' ! For command parameters
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &   item
      logical, intent(in) ::
     &   s1,s2
      type(spin_info2), intent(inout) ::
     &   t_spin(3)
      type(x_inter), intent(inout) ::
     &   x_dict(MAXX)
      type(inter_spin_cases), intent(inout) ::
     &     inter_spin_dict
      type(tensor_names), intent(in) ::
     &     itf_names            ! contains renaming information
      integer, intent(in) ::
     &   ntest

      integer, parameter :: maxxstr = 5
      character(len=MAXLEN_BC_LABEL) ::
     &   nres, nt1, nt2,         ! Name of tensors involved in the contraction
     &   nt1x(maxxstr), nt2x(maxxstr)
      character(len=INDEX_LEN) ::
     &   slabel1, slabel2, slabel3,
     &   new_idx1, new_idx2, new_idx3,
     &   new_idx1_x(maxxstr), new_idx2_x(maxxstr)
      character(len=264) ::
     &   itf_line,          ! Line of ITF code
     &   st1, st2,          ! Name of spin summed tensors + index
     &   tst1, tst2,          ! Name of spin summed tensors + index
     &   xst1, xst2           ! Name of spin summed tensors + index
      character(len=2) ::
     &   equal_op           ! ITF contraction operator; ie. +=, -=, :=
      character(len=MAXLEN_BC_LABEL) ::
     &   sfact,             ! String representation of factor
     &   sfact_star,        ! String representation of factor formatted for output
     &   k4e_no,            ! Counter of K4E tensors
     &   nx, nh
      integer ::
     &   i,j,k,nxstr1,nxstr2
      real(8) ::
     &   c_fact               ! Copy of orginal factor
      logical ::
     &   new_j,
     &   old,
     &   found

      character(len=MAXLEN_BC_LABEL), external ::
     &     rename_tensor

      new_idx1 = item%idx1
      new_idx2 = item%idx2
      new_idx3 = item%idx3
      new_idx1_x = ""
      new_idx2_x = ""
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


      ! For pure aaaa or bbbb cases create exchange string
      if (.not.item%inter(1)) then
        nxstr1 = maxxstr        ! assure subroutine that enough space is there
        call create_exchange_string(item%rank1,s1,new_idx1,
     &                              nxstr1,new_idx1_x)
      else
        nxstr1 = 0
      end if

      if (.not.item%inter(2)) then
        nxstr2 = maxxstr
        call create_exchange_string(item%rank2,s2,new_idx2,
     &                              nxstr2,new_idx2_x)
      else
        nxstr2 = 0
      end if

      if (ntest>=100) then
         write(item%out,*) "---------------------------"
         write(item%out,*) "After convert_to_abab_block"
         write(item%out,*) "Idx1: ", new_idx1, new_idx1_x(1:nxstr1)
         write(item%out,*) "Idx2: ", new_idx2, new_idx2_x(1:nxstr2)
         write(item%out,*) "Idx3: ", new_idx3
         write(item%out,*) "Factor: ", c_fact
         write(item%out,*) "---------------------------"
      end if

      ! If dealing with a result line, check the spin of any
      ! intermediates that contribute to it. If the intermediates begin
      ! with a beta spin, then the spin will be flipped,
      ! i.e. baab -> abba. These spin cases are the same, so we don't
      ! need to contruct both intermediates.
      !if (.not. item%inter(3)) then
         if (item%inter(1) .and. item%rank1 /= 0) then
            if (t_spin(1)%spin(1,1) == 2) then

               do i = 1, 2
                  do j = 1, item%rank1/2
                     if (t_spin(1)%spin(i,j) == 2) then
                        t_spin(1)%spin(i,j) = 1
                     else
                        t_spin(1)%spin(i,j) = 2
                     end if
                  end do
               end do

            end if
         end if
         if (item%inter(2) .and. item%rank2 /= 0) then
            if (t_spin(2)%spin(1,1) == 2) then

               do i = 1, 2
                  do j = 1, item%rank2/2
                     if (t_spin(2)%spin(i,j) == 1) then
                        t_spin(2)%spin(i,j) = 2
                     else
                        t_spin(2)%spin(i,j) = 1
                     end if
                  end do
               end do

            end if
         end if
      !end if


      ! Change names of specific tensors
      nres = rename_tensor(item%label_res,
     &                     item%rank3, item%nops3, itf_names)
      nt1  = rename_tensor(item%label_t1,
     &                     item%rank1, item%nops1, itf_names)
      nt2  = rename_tensor(item%label_t2,
     &                     item%rank2, item%nops2, itf_names)

      ! copy name to exchange contributions (some may get renamed)
      nt1x(1:nxstr1) = nt1
      nt2x(1:nxstr2) = nt2

      ! Add intermediate spin strings to names
      if (item%inter(1)) then
         call inter_spin_name(t_spin(1)%spin,item%rank1/2, slabel1)
         nt1 = trim(nt1)//trim(slabel1)
         call reorder_intermediate(item%rank1,new_idx1,item%nops1)
       else
         call reorder_to_slots(item%int(1),item%rank1,
     &                         new_idx1,
     &                         nt1,item%nops1,item%out)
         do i = 1, nxstr1
           call reorder_to_slots(item%int(1),item%rank1,
     &                         new_idx1_x(i),
     &                         nt1x(i),item%nops1,item%out)
         end do
      end if

      if (item%inter(2)) then
         call inter_spin_name(t_spin(2)%spin,item%rank2/2,slabel2)
         nt2 = trim(nt2)//trim(slabel2)
         call reorder_intermediate(item%rank2,new_idx2,item%nops2)
       else
         call reorder_to_slots(item%int(2),item%rank2,
     &                         new_idx2,
     &                         nt2,item%nops2,item%out)
         do i = 1, nxstr2
           call reorder_to_slots(item%int(2),item%rank2,
     &                         new_idx2_x(i),
     &                         nt2x(i),item%nops2,item%out)
         end do
      end if

      ! If a residual (or intermediate) depends on a non-existant intermediate spin case
      ! (ie. a non-totally spin conserving case = 0), then skip printing
      ! out the line. This may not be the case for ionisation/spin flip
! cases.
c     dbg
c      write(item%out,*) 'Names: ',trim(nt1),' ',trim(nt2)
c      if (trim(nt1)==' STIN0002abba') then
c        write(item%out,*) 'Now it becomes interesting'
c        write(item%out,*) 'inter_spin_dict:'
c        do i = 1,  inter_spin_dict%ncase
c          write(item%out,'(i4,a)') i,trim(inter_spin_dict%names(i))
c        end do
c      end if
c dbg
      if (item%inter(1))then
         if (inter_spin_dict%ncase > 0) then
            found = .false.
            do i = 1, inter_spin_dict%ncase
               if (nt1 == inter_spin_dict%names(i)) then
                  found = .true.
                  exit
               end if
            end do

            if (.not. found) return
         end if
      end if
      if (item%inter(2))then
         if (inter_spin_dict%ncase > 0) then
            found = .false.
            do i = 1, inter_spin_dict%ncase
               if (nt2 == inter_spin_dict%names(i)) then
                  found = .true.
                  exit
               end if
            end do

            if (.not. found) return
         end if
      end if

      if (item%inter(3)) then
         call inter_spin_name(t_spin(3)%spin,item%rank3/2,slabel3)
         nres = trim(nres)//trim(slabel3)
         call reorder_intermediate(item%rank3,new_idx3,item%nops3)

         ! Save spin cases name to dictionary
         ! When the residual line is printed, the intermediate will be
         ! checked against this list to see if it exists (ie. is 0 or not)
         inter_spin_dict%names(inter_spin_dict%ncase+1) = nres
         inter_spin_dict%ncase = inter_spin_dict%ncase + 1

      else
         call reorder_to_slots(item%int(3),item%rank3,
     &                         new_idx3,
     &                         nres,item%nops3,item%out)
      end if


c      ! Reorder integrals into fixed slot order
c      call reorder_integral(item%int(1),item%rank1,new_idx1,
c     &                      new_j,nt1,item%nops1)
c      call reorder_integral(item%int(2),item%rank2,new_idx2,
c     &                      new_j,nt2,item%nops2)
c
c      if (.not. item%inter(1) .and. .not. item%int(1)) then
c         call reorder_amp(item%rank1, new_idx1)
c      end if
c      if (.not. item%inter(2) .and. .not. item%int(2)) then
c         call reorder_amp(item%rank2, new_idx2)
c      end if
c      if (.not. item%inter(3) .and. .not. item%int(3)) then
c         call reorder_amp(item%rank3, new_idx3)
c       end if

      if (ntest>=100) then
         write(item%out,*) "---------------------------"
         write(item%out,*) "After reorder()"
         write(item%out,*) "Idx1: ", new_idx1, new_idx1_x(1:nxstr1)
         write(item%out,*) "Idx2: ", new_idx2, new_idx2_x(1:nxstr2)
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

c      ! Change tensor to spatial orbital quantity, unless it is an
c      ! intermediate
c      call spatial_string(st1,new_idx1,nt1,s1,item%inter(1),item%rank1,
c     &                1,item%binary,item%int(1),item%nops1,new_j,
c     &                item%out)
c      call spatial_string(st2,new_idx2,nt2,s2,item%inter(2),item%rank2,
c     &                2,item%binary,item%int(2),item%nops2,new_j,
c     &                item%out)

      ! create operands for ITF string, like K[abij] or (K[abij]-J[abji])
      call create_tensor_string(st1,nt1,new_idx1,
     &                              nt1x,new_idx1_x,nxstr1)
      if (item%binary) then
        call create_tensor_string(st2,nt2,new_idx2,
     &                              nt2x,new_idx2_x,nxstr2)
      else
        st2=''
      end if

      ! Create intermediate instead of brackets
      if (item%tasks) then

         if (s1 .and. .not. item%inter(1)) then

            ! Need to check if already declared this X intermediate
            ! Names an operator info are saved in x_dict between loops
            call check_x(item, x_dict, tst1, new_idx1, .true., old)

            if (.not. old) then
               ! Print out new X intermediate
               write(nh,*) item%cntr(4)
               nx = trimal('X'//trimal(nh))

               call create_Xtensor_lines(tst1,xst1,nx,
     &                       nt1,new_idx1,nt1x,new_idx1_x,nxstr1)

c               tst1='X'//trimal(nx)//'['//trim(new_idx1)//']'
c               call spatial_string2(xst1,new_idx1,nt1,s1,item%inter(1),
c     &                      item%rank1,
c     &                      1,item%binary,item%int(1),item%nops1,new_j,
c     &                      item%cntr(4),item%out)

               write(item%out,'(a)') trim(xst1)

               ! Increase X counter
               item%cntr(4) = item%cntr(4) + 1
            end if

         else
            tst1 = st1
         end if


         if (item%binary .and. s2 .and. .not. item%inter(2)) then

            call check_x(item, x_dict, tst2, new_idx2, .false., old)

            if (.not. old) then
               write(nh,*) item%cntr(4)
               nx = trimal('X'//trimal(nh))

               call create_Xtensor_lines(tst2,xst2,nx,
     &                       nt2,new_idx2,nt2x,new_idx2_x,nxstr2)

c               tst2='X'//trimal(nx)//'['//trim(new_idx2)//']'
c               call spatial_string2(xst2,new_idx2,nt2,s2,item%inter(2),
c     &                      item%rank2,
c     &                      2,item%binary,item%int(2),item%nops2,new_j,
c     &                      item%cntr(4),item%out)

               write(item%out,'(a)') trim(xst2)

               item%cntr(4) = item%cntr(4) + 1
            end if

         else
            tst2 = st2
         end if

      end if


      ! Add factor to scalar result (ie. the energy) cases (going to skip half the spin
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
            write(sfact,'(f24.15)') c_fact

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
      if (item%tasks) then
         itf_line='.'//trimal(nres)//
     &       '['//trim(new_idx3)//'] '//equal_op//' '//
     &       trim(sfact_star)//trimal(tst1)//' '//trimal(tst2)
      else
         itf_line='.'//trimal(nres)//
     &      '['//trim(new_idx3)//'] '//equal_op//' '//
     &      trim(sfact_star)//trimal(st1)//' '//trimal(st2)
      end if


      ! Print K4E line instead
      if (item%k4e_line) then
         write(k4e_no,*) item%cntr(3)

         itf_line='.'//trimal(nres)//
     &      '['//trim(new_idx3)//'] '//equal_op//' '//
c     &      trim(sfact_star)//'K4E'//trim(adjustl(k4e_no))//
c without number:
     &      trim(sfact_star)//'K4E'//
     &      '['//trim(new_idx3)//']'
      end if

      ! Print to output file
      write(item%out,'(a)') trim(itf_line)

      ! Increment number of printed spn cases
      item%spin_cases = item%spin_cases + 1

      return
      end


*----------------------------------------------------------------------*
      subroutine check_x(item, x_dict, tst, idx, t1, old)
*----------------------------------------------------------------------*
!     Check if already delcared an X intermediate, if not return
!     new=false and add info to x_dict
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &   item
      type(x_inter), intent(inout) ::
     &   x_dict(MAXX)
      character(len=264) ::
     &   tst
      character(len=INDEX_LEN), intent(in) ::
     &   idx
      logical, intent(in) ::
     &   t1
      logical, intent(inout) ::
     &   old

      integer ::
     &   i,j,k,
     &   e(ngastp,2), c(ngastp,2)
      character(len=MAXLEN_BC_LABEL) ::
     &   label
      character(len=25) ::
     &   nx

      if (t1) then
         label = item%label_t1
         e = item%e1
         c = item%c
      else
         label = item%label_t2
         e = item%e2
         do j = 1, ngastp
            c(j,1) = item%c(j,2)
            c(j,2) = item%c(j,1)
         end do
      end if

      old = .false.
      do i = 1, item%cntr(4)
         if (label==x_dict(i)%label) then

            old = .true.
            do j = 1, ngastp
               do k = 1, 2
                  ! Compare operator numbers
                  if (x_dict(i)%ops(j,k)/=e(j,k)+c(j,k)) then
                     old = .false.
                     exit
                  end if
               end do
            end do

            if (old) then
               ! We have declared this X inter before, so just use
               ! the old X number
               write(nx,*) x_dict(i)%n
               tst='X'//trimal(nx)//'['//trim(idx)//']'
               return
            end if

         end if
      end do

      if (.not. old) then
         ! If a new X intermediate, save info to x_dict
         x_dict(item%cntr(4))%n = item%cntr(4)
         x_dict(item%cntr(4))%label = label
         do j = 1, ngastp
            do k = 1, 2
               x_dict(item%cntr(4))%ops(j,k)=e(j,k)+c(j,k)
            end do
         end do
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine reorder_to_slots(integral,rank,idx,label,nops)
*----------------------------------------------------------------------*
*     analyse index string (on idx) and resort to match the actual
*     storage convention of the associated tensors
*     we only interchange indices according to the generally underlying
*     symmetries
*     .not.integral: assuming 1<>2 etc. particle exchange symmetry
*               we will also interchange C and A if there are more
*               external indices on A (deexcitation amplitudes)
*     integral: assuming 1<>2 and individual C<>A symmetry
*               map to special integrals, special orders if required
*
*     andreas, dec. 2020, based on routines by j.a.black
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
     &   integral
      character(len=MAXLEN_BC_LABEL), intent(inout) ::
     &   label

       integer ::
     &   itype(rank),
     &   ii, jj,
     &   itmp, exc, exa
      character(len=1) ::
     &   tmp
      character(len=INDEX_LEN) ::
     &   tstr
      logical ::
     &   all_same

      if (rank == 0 .or. rank == 2) return
      if (mod(rank,1)>0) call quit(1,'reorder_to_slots','odd rank?')

      ! get itypes
      do ii = 1, rank
        itype(ii) = get_itype(idx(ii:ii))
      end do

      ! round one: C<>A exchange?
      if (integral) then        ! for integrals, we assume C<>A symmetry for all pairs
        do ii = 1, rank/2
          if (itype(ii)>itype(ii+rank/2)) then
            ! Swap creation and annhilation
            tmp = idx(ii:ii)
            idx(ii:ii) = idx(ii+rank/2:ii+rank/2)
            idx(ii+rank/2:ii+rank/2) = tmp
            ! Update itype
            itmp = itype(ii)
            itype(ii) = itype(ii+rank/2)
            itype(ii+rank/2) = itmp
          end if
        end do
      else                      ! otherwise check on which part of the string more external indices are
        exc = 0
        exa = 0
        do ii = 1, rank/2
          if (itype(ii)==1)        exc = exc+1
          if (itype(ii+rank/2)==1) exa = exa+1
        end do
        if (exa>exc) then
          tstr(1:INDEX_LEN) = " "
          ! permute all C and A
          tstr(1:rank/2) = idx(rank/2+1:rank)
          tstr(rank/2+1:rank) = idx(1:rank/2)
          idx = tstr
          ! set itype from new
          do ii = 1, rank
            itype(ii) = get_itype(idx(ii:ii))
          end do
        end if
      end if

      ! check if all C and A indices or of same itype
      all_same = .true.
      do ii = 2, rank/2
        all_same = all_same.and.itype(ii)==itype(1)
      end do
      do ii= rank/2+2, rank
        all_same = all_same.and.itype(ii)==itype(rank/2+1)
      end do

      ! in this case we can exit
      if (all_same) return

      ! beyond this point, we currently can only handle rank 4 correctly
      ! rank 6+ requires proper sort algorithm
      if (rank>4) call quit(1,'reorder_to_slots',
     &                        'cannot handle this yet: '//trim(idx))

      ! explicit rank 4 code
      ! compare C, if C equal compare A
      if ( itype(1)>itype(2).or.
     &    (itype(1)==itype(2).and.itype(3)>itype(4))) then
        ! switch
        tmp = idx(1:1); idx(1:1) = idx(2:2); idx(2:2) = tmp
        tmp = idx(3:3); idx(3:3) = idx(4:4); idx(4:4) = tmp
        itmp = itype(1); itype(1) = itype(2); itype(2) = itmp
        itmp = itype(3); itype(3) = itype(4); itype(4) = itmp
      end if

      ! handle a few special cases
      if (integral) then
        ! all integrals that have an ext. orbital in position 1 and 3
        ! (this includes also 3-ext integrals) have a 1122 slot ordering
        ! also the cases ecac and acac
        if ((itype(1)==1.and.itype(3)==1).or.
     &      all(itype==(/1,3,2,3/)).or.
     &      all(itype==(/2,3,2,3/))) then
          ! switch (abcd) to (acbd)
          tmp = idx(2:2); idx(2:2) = idx(3:3); idx(3:3) = tmp
          itmp = itype(2); itype(2) = itype(3); itype(3) = itmp
          ! unless the case of 3-externals, we have to rename K to J
          if (nops(2).ne.3) label = 'J'
          ! Need to have special case of J:eeca not J:eeac
          if (all(itype==(/1,1,2,3/))) then
            tmp = idx(1:1); idx(1:1) = idx(2:2); idx(2:2) = tmp
            tmp = idx(3:3); idx(3:3) = idx(4:4); idx(4:4) = tmp
          end if
        end if
        ! Need to have special case of K:ccca not K:accc
        if (all(itype==(/2,3,3,3/))) then
          ! we can make a cyclic permutation
          tmp = idx(1:1)
          idx(1:1) = idx(2:2); idx(2:2) = idx(3:3)
          idx(3:3) = idx(4:4); idx(4:4) = tmp
        end if
      end if

      end subroutine
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

      symmetric = .true.  ! <- unused?
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
         do i = 1, rank/2-1   ! <- may fail for rank 6, e.g. for a c e
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
      subroutine create_tensor_string(st,nt,idx,nt_x,idx_x,nidxx)
*----------------------------------------------------------------------*
*     construct tensor
*        string "<nt>[<idx>]" if nidxx==0
*            or "(<nt>[<idx>] - <nt_x>[<idx_x])" if nidxx==1
*            or "(<nt>[<idx>] + ... + ... - ... - ... - ...)" if nidxx==5
*
*     andreas, dec 2020 based on routines by j.a.black
*----------------------------------------------------------------------*

      use itf_utils

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=264), intent(inout) ::
     &     st                   ! Name of spin summed tensors + index
      integer, intent(in) ::
     &     nidxx                ! Number of exchange contributions
      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &     nt                   ! Name of tensors involved in the contraction
      character(len=INDEX_LEN), intent(in) ::
     &     idx                  ! Index of tensor
      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &     nt_x(nidxx)          ! Name of tensors for exchange part
      character(len=INDEX_LEN), intent(in) ::
     &     idx_x(nidxx)         ! Index of tensor for exchange part


      ! hard coded ... could generalize if array of signs is provided
      if (nidxx==0) then
        st=trimal(nt)//'['//trim(idx)//']'
      else if (nidxx==1) then
        st='('//trimal(nt)//'['//trim(idx)//'] - '//
     &          trimal(nt_x(1))//'['//trim(idx_x(1))//'])'
      else if (nidxx==5) then
        st='('//trimal(nt)//'['//trim(idx)//'] + '//
     &          trimal(nt_x(1))//'['//trim(idx_x(1))//'] + '//
     &          trimal(nt_x(2))//'['//trim(idx_x(2))//'] - '//
     &          trimal(nt_x(3))//'['//trim(idx_x(3))//'] - '//
     &          trimal(nt_x(4))//'['//trim(idx_x(4))//'] - '//
     &          trimal(nt_x(5))//'['//trim(idx_x(5))//'])'
      else
        call quit(1,'create_tensor_string','not prepared for this case')
      end if

      end subroutine

*----------------------------------------------------------------------*
      subroutine create_Xtensor_lines(tst,st,nx,nt,idx,nt_x,idx_x,nidxx)
*----------------------------------------------------------------------*
*     construct line(s) defining X tensor, one for adressing
*       tst = "<nx>[<idx>]"
*     and one for defining
*       st  = ".<nx>[<idx>] += <nt>[<idx>] <newline>
*              .<nx>[<idx>] -= <nt_x>[<idx_x]"
*
*     andreas, dec 2020 based on routines by j.a.black
*----------------------------------------------------------------------*

      use itf_utils

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=264), intent(inout) ::
     &     tst,                 ! Name of intermediate w indices
     &     st                   ! Name of spin summed tensors + index
      integer, intent(in) ::
     &     nidxx                ! Number of exchange contributions
      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &     nx,                  ! Name of intermediate
     &     nt                   ! Name of tensors involved in the contraction
      character(len=INDEX_LEN), intent(in) ::
     &     idx                  ! Index of tensor
      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &     nt_x(nidxx)          ! Name of tensors for exchange part
      character(len=INDEX_LEN), intent(in) ::
     &     idx_x(nidxx)         ! Index of tensor for exchange part


      tst = trimal(nx)//'['//trim(idx)//']'

      ! hard coded ... could generalize if array of signs is provided
      if (nidxx==1) then
        st='.'//trimal(nx)//'['//trim(idx)//'] += '//
     &          trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &     '.'//trimal(nx)//'['//trim(idx)//'] -= '//
     &          trimal(nt_x(1))//'['//trim(idx_x(1))//']'
      else
        call quit(1,'create_Xtensor_lines','not prepared for this case')
      end if

      end subroutine


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
                        ! K:eaac - K:ecaa   <---- ?? exchange should be: - K:eaca
                        st='('//trimal(nt)//'['//trim(idx)//']'//' - '//
     &                      trimal(nt)//'['//
! ???     &                 f_index(idx,hrank,.false.,.false.,.false.,.true.)
     &                 f_index(idx,hrank,.true.)
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
      subroutine spatial_string2(st,idx,nt,spin,inter,rank,tensor,
     &                          binary,
     &                          integral,nops,j_int,nx,lulog)
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
     &   nx,
     &   lulog       ! Logfile

      integer ::
     &   hrank       ! Half rank
      character(len=25) ::
     &   x

      hrank = rank / 2

      write(x,*) nx
      st = ''

      if (spin .and. .not.inter) then
         ! Pure spin
         select case (rank)
            case (4)

               if (integral) then
                  if ((nops(1)==3 .or. nops(2)==3 .or. nops(3)==3)) then
                   ! Three-somthing integral (ie. K:eccc, K:accc, ...)

                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             'K'//'['//f_index(idx,hrank,.false.,.true.)//']'

                  else if (j_int) then
                     if (trim(nt)=='K') then
                        ! Need (K:eecc - J:eecc)

                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             'J'//'['//f_index(idx,hrank,.true.)//']'

                     else if (trim(nt)=='J') then

                        ! Need (J:eecc - K:eecc)
                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             'K'//'['//f_index(idx,hrank,.true.)//']'

                     end if

                  else if ((nops(1)==1.and.nops(2)==2.and.nops(3)==1))
     &              then
                    ! K:eeac - K:eeac
                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             trimal(nt)//'['//f_index(idx,hrank)//']'

                  else if ((nops(1)==1.and.nops(2)==1.and.nops(3)==2))
     &              then
                     if (get_itype(idx(2:2))==2 .and.
     &                   get_itype(idx(3:3))==2) then
                        ! K:eaac - K:ecaa

                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             trimal(nt)//'['//
     &             f_index(idx,hrank,.false.,.false.,.false.,.true.)//
     &             ']'

                     else

                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             trimal(nt)//'['//f_index(idx,hrank,.true.)//']'

                     end if

                  else

                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             trimal(nt)//'['//f_index(idx,hrank,.true.)//']'

                  end if

               else
                  ! Amplitude or density
                  if ((nops(1)==2.and.nops(2)==1.and.nops(3)==1)) then
                    ! T:eacc

                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             trimal(nt)//'['//f_index(idx,hrank,.true.)//']'

                  else if ((nops(1)==1.and.nops(2)==1.and.
     &                      nops(3)==2)) then
                    ! T:eaca or T:eaac

                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             trimal(nt)//'['//f_index(idx,hrank,.true.)//']'

                  else

                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             trimal(nt)//'['//f_index(idx,hrank)//']'

                  end if
               end if

            case (6)
               ! TODO: ignore for now...
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

                   st='.X'//trimal(x)//'['//trim(idx)//'] += '//
     &             trimal(nt)//'['//trim(idx)//']'//new_line('a')//
     &             '.X'//trimal(x)//'['//trim(idx)//'] -= '//
     &             trimal(nt)//'['//f_index(idx,hrank)//']'

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
      subroutine create_index_str2(idx, cnt, ex, c_shift, e_shift, rank,
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
     &   rank
      integer, intent(inout) ::
     &   c_shift(4),
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


      idx%str = ""
      idx%itype = 0

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

      c_shift = cnt_shift

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
      subroutine revert_index_str(str,nidx)
*----------------------------------------------------------------------*
*     revert squence (for transposed operator)
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(inout) :: str
      integer, intent(in) :: nidx

      integer :: idx
      type(index_str) :: scr_str
      
      call init_index_str(scr_str,nidx,1) ! cnt_poss is dummy only

      do idx = 1, nidx
        scr_str%str(idx) = str%str(idx)
        scr_str%itype(idx) = str%itype(idx)
      end do

      do idx = 1, nidx
        str%str(idx) = scr_str%str(nidx+1-idx)
        str%itype(idx) = scr_str%itype(nidx+1-idx)
      end do

      call deinit_index_str(scr_str)

      return
      end subroutine
      
*----------------------------------------------------------------------*
      subroutine set_index_str(str,index_info,nidx)
*----------------------------------------------------------------------*
*     use the information on index_info to create an index string
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(inout) :: str
      integer :: nidx, index_info(nidx)

      character, parameter ::
     &     label(10,4) = reshape(
     &                     (/'i','j','k','l','m','n','o',' ',' ',' ',
     &                       'a','b','c','d','e','f','g','h',' ',' ',
     &                       'p','q','r','s','t','u','v','w','x','y',
     &                       'P','Q','R','S','T','U',' ',' ',' ',' '/),
     &                     (/10,4/))
      integer, parameter ::
     &     shift_type(4) = (/3,1,2,4/)  ! 'canonical' values

      integer :: ii, type, idx

      do ii = 1, nidx
        type = index_info(ii) / 1000
        idx  = mod(index_info(ii),1000)
        str%str(ii) = label(idx,type)
        str%itype(ii) = shift_type(type)
      end do

      end subroutine

*----------------------------------------------------------------------*
      subroutine normalize_index_str(str,nidx,occ_op,nj_op)
*----------------------------------------------------------------------*
*     for operators that orignally consisted of several vertices
*     we have to postprocess the index, such that all indices originally
*     associated with creations are on the left
*     we assume that this routine is called before cnt_poss is set!
*
*     andreas, dec 2020
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(index_str), intent(inout) :: str
      integer :: nidx, nj_op, occ_op(ngastp,2,nj_op)

      integer :: idx, nidx_c1, nidx_a1, nidx_v2, nidx_c3, nidx_a3,
     &           nidx0, ij, off_a, off_c
      integer :: nidx_c(nj_op), nidx_a(nj_op)

      type(index_str) :: scr_str
      integer :: ireo(nidx)

      call init_index_str(scr_str,nidx,1) ! cnt_poss is dummy only


c      write(lulog,*) 'renormalize:'
c      call wrt_occ_n(lulog,occ_op,nj_op)
c      write(lulog,'(1x,10a1)') str%str(1:nidx)

      if (nj_op.eq.2) then
!     leave creations of first vertex in the front and move
!     annihilations of first vertex behind second vertex -
!     without changing their sequence
        nidx_c1 = sum(occ_op(1:ngastp,1,1))
        nidx_a1 = sum(occ_op(1:ngastp,2,1))
        nidx_v2 = sum(occ_op(1:ngastp,1:2,2))
        do idx = 1, nidx_c1
          ireo(idx) = idx  ! C1 stay at same place
        end do
        do idx = nidx_c1+1, nidx_c1+nidx_a1
          ireo(idx) = idx+nidx_v2  ! A1 now behind C2A2
        end do
        do idx = nidx_c1+nidx_a1+1, nidx
          ireo(idx) = idx-nidx_a1  ! C2A2 now directly after C1 (shift by length of A1)
        end do
      else if (nj_op.eq.3) then
!     3 vertices:
!     move all A from v1 behind v3 and all C from v3 between v1 and v2
        nidx_c1 = sum(occ_op(1:ngastp,1,1))
        nidx_a1 = sum(occ_op(1:ngastp,2,1))
        nidx_v2 = sum(occ_op(1:ngastp,1:2,2))
        nidx_c3 = sum(occ_op(1:ngastp,1,3))
        nidx_a3 = sum(occ_op(1:ngastp,2,3))
        do idx = 1, nidx_c1
          ireo(idx) = idx  ! C1 in same place
        end do
        nidx0 = nidx_c1
        do idx = nidx0+1, nidx0+nidx_a1
          ireo(idx) = idx+nidx_c3+nidx_v2+nidx_a3 ! A1 now at the end
        end do
        nidx0 = nidx0+nidx_a1
        do idx = nidx0+1, nidx0+nidx_v2
          ireo(idx) = idx-nidx_a1+nidx_c3 ! C2A2 now after C1C3 (A1 is now away, but C3 is there)
        end do
        nidx0 = nidx0+nidx_v2
        do idx = nidx0+1, nidx0+nidx_c3
          ireo(idx) = idx-nidx_a1-nidx_v2 ! C3 is now after C1 (so A1 and C2A2 have moved away)
        end do
        nidx0 = nidx0+nidx_c3
        do idx = nidx0+1, nidx
          ireo(idx) = idx-nidx_a1 ! A3 is now directly after C2A2 (A1 has moved away)
        end do
      else
        ! general strategy: compatible to 2-vertex case
        ! but would give different order for 3-vertex case (can break pairing of middle vertex)
        ! should not influence numerical result, efficiency of final code has to be checked
        do ij = 1, nj_op
          nidx_c(ij) = sum(occ_op(1:ngastp,1,ij))
          nidx_a(ij) = sum(occ_op(1:ngastp,2,ij))
        end do
        nidx0 = 0
        off_c = 0
        off_a = sum(nidx_c(1:nj_op))+sum(nidx_a(2:nj_op))
        do ij = 1, nj_op
          do idx = 1, nidx_c(ij)
            ireo(nidx0+idx) = off_c + idx
          end do
          nidx0 = nidx0+nidx_c(ij)
          do idx = 1, nidx_a(ij)
            ireo(nidx0+idx) = off_a + idx
          end do
          nidx0 = nidx0+nidx_a(ij)
          off_c = off_c+nidx_c(ij)
          if (ij<nj_op) off_a = off_a-nidx_a(ij+1)
        end do
      end if

      do idx = 1, nidx
        scr_str%str(ireo(idx)) = str%str(idx)
        scr_str%itype(ireo(idx)) = str%itype(idx)
      end do

      str%str = scr_str%str
      str%itype = scr_str%itype

c      write(lulog,'(1x,10a1)') str%str(1:nidx)

      call deinit_index_str(scr_str)

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine assign_index(item,contr_info,ntest)
*----------------------------------------------------------------------*
!     Assign an ITF index string to each tensor in a line
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &     item                 ! ITF binary contraction
      type(binary_contr), intent(in) ::
     &     contr_info
      integer, intent(in) ::
     &   ntest

      integer ::
     &   n_cnt,         ! Number of contraction operators
     &   i, idx, jdx, ioff, nidx, ncnt
      character(len=INDEX_LEN) ::
     &   s1, s2, s3  ! Tmp ITF index strings
      character(len=1) ::
     &   tmp       ! Scratch space to store index letter
      type(index_str) ::
     &   str1,          ! Stores 'normal ordered' index string for T1
     &   str2,          ! Stores 'normal ordered' index string for T2
     &   str3           ! Stores 'normal ordered' index string for Res

      if (ntest>=100) call debug_header("assign_index", item%out)


      ! Set number of contraction indicies
      n_cnt = item%contri

      ! Allocate index_str objects
      call init_index_str(str1, item%rank1, n_cnt)
      call init_index_str(str2, item%rank2, n_cnt)
      call init_index_str(str3, item%rank3, n_cnt)

c     dbg
c      write(item%out,'(1x,"ranks:",3i6)')
c     &     item%rank1, item%rank2, item%rank3
c      write(item%out,'(1x,"index info: ",3i6)')
c     &     contr_info%itf_index_info(1:3)
c      write(item%out,'(1x,"o1:  ",10i6)')
c     &     contr_info%itf_index_info(3+1:3+contr_info%itf_index_info(1))
c     dbg

      ioff = 3 ! first three entries are dimensions
      nidx = contr_info%itf_index_info(1)
      if (nidx.ne.item%rank1) then
        write(item%out,'(1x,"nidx/rank1: ",i4,"/",i4)') nidx,item%rank1
        call quit(1,'assign_index','rank mismatch (1)')
      end if
      call set_index_str(str1,
     &     contr_info%itf_index_info(ioff+1:ioff+nidx), nidx)
      if (contr_info%nj_op1.gt.1)
     &     call normalize_index_str(str1,nidx,
     &     contr_info%occ_op1,contr_info%nj_op1)
      if (contr_info%tra_op1) call revert_index_str(str1,nidx) ! invert sequence if transposition was requested
      
      ioff = ioff+nidx  ! shift offset by indices of tensor 1
      nidx = contr_info%itf_index_info(2)
      if (nidx.ne.item%rank2)
     &     call quit(1,'assign_index','rank mismatch (2)')
      call set_index_str(str2,
     &     contr_info%itf_index_info(ioff+1:ioff+nidx), nidx)
      if (contr_info%nj_op2.gt.1)
     &     call normalize_index_str(str2,nidx,
     &     contr_info%occ_op2,contr_info%nj_op2)
      if (contr_info%tra_op2) call revert_index_str(str2,nidx) ! invert sequence if transposition was requested

      ioff = ioff+nidx  ! shift offset by indices of tensor 2
      nidx = contr_info%itf_index_info(3)
      if (nidx.ne.item%rank3)
     &     call quit(1,'assign_index','rank mismatch (3)')
      call set_index_str(str3,
     &     contr_info%itf_index_info(ioff+1:ioff+nidx), nidx)
      if (contr_info%nj_res.gt.1)
     &     call normalize_index_str(str3,nidx,
     &     contr_info%occ_res,contr_info%nj_res)
      if (contr_info%tra_res) call revert_index_str(str3,nidx) ! invert sequence if transposition was requested

! get cnt_poss from index matching
      ncnt = 0
      do idx = 1, item%rank1
        do jdx = 1, item%rank2
          if (str1%str(idx)==str2%str(jdx)) then
            ncnt = ncnt+1
            str1%cnt_poss(ncnt) = idx
            str2%cnt_poss(ncnt) = jdx
          end if
        end do
      end do

      if (ntest>=100) then
         write(item%out,*) "After set_index_str"
         write(item%out,*) "T1: {", str1%str, "}"
         write(item%out,'(3x,a,10i4)') "itpye:    ", str1%itype
         write(item%out,'(3x,a,10i4)')   "cnt_poss: ", str1%cnt_poss
         write(item%out,*) "T2: {", str2%str, "}"
         write(item%out,'(3x,a,10i4)') "itpye:    ", str2%itype
         write(item%out,'(3x,a,10i4)')   "cnt_poss: ", str2%cnt_poss
      end if



      ! Rearrange intermediate index to match previously declared inter
      ! index. This uses the itype array from the previous lines
      if (item%inter(1)) then
         call arrange_inter_itype(item%rank1,item%rank2,str1,str2,
     &                           item%itype, item%label_t1, item%out)
         if (ntest>=100) then
            write(item%out,*) "After arrange_inter_itype, inter(1)"
            write(item%out,*) "T1: {", str1%str, "}"
            write(item%out,*) "itpye: {", str1%itype, "}"
            write(item%out,*) "cnt_poss: ", str1%cnt_poss
            write(item%out,*) "previous itype: ", item%itype
         end if
      end if
      if (item%inter(2)) then
         call arrange_inter_itype(item%rank2,item%rank1,str2,str1,
     &                           item%itype, item%label_t2, item%out)
         if (ntest>=100) then
            write(item%out,*) "After arrange_inter_itype, inter(2)"
            write(item%out,*) "T2: {", str2%str, "}"
            write(item%out,*) "itpye: {", str2%itype, "}"
            write(item%out,*) "cnt_poss: ", str2%cnt_poss
            write(item%out,*) "previous itype: ", item%itype
         end if
      end if

      if (ntest>=100) then
         write(item%out,*) "Index strings in normal order"
         write(item%out,*) "Result string in arbitrary order"
         write(item%out,*) "T1: {", str1%str, "}"
         write(item%out,*) "T2: {", str2%str, "}"
         write(item%out,*) "Res: {", str3%str, "}"
         write(item%out,*) "Contraction T1: ", str1%cnt_poss
         write(item%out,*) "Contraction T2: ", str2%cnt_poss
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
         tmp = '+'
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

      call deinit_index_str(str1)
      call deinit_index_str(str2)
      call deinit_index_str(str3)

      return
      end


*----------------------------------------------------------------------*
      subroutine arrange_inter_itype(rank1, rank2, idx1, idx2, itype,
     &                               label,itflog)
*----------------------------------------------------------------------*
!     rearrange intermeidate index according to the itype of the
!     previously declared intermediate
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   rank1,
     &   rank2,
     &   itype(maxint, index_len),
     &   itflog
      type(index_str), intent(inout) ::
     &   idx1,
     &   idx2
      character(len=maxlen_bc_label), intent(in) ::
     &   label

      integer ::
     &   i, j,
     &   shift,
     &   ninter, nfound
      character (len=index_len) ::
     &   tstr
      logical ::
     &   same

      ! find itype corresponding to numbered intermediate
      i = len(trim(label))
      read(label(i:i),*) ninter

      same = .true.
      do i = 1, rank1
         ! check if intermediate already has correct index itype
         if (idx1%itype(i) /= itype(ninter,i)) then
            same = .false.
         end if
      end do

      if (same) return

      nfound = 0
      do i = 1, rank1/2
         do j = 1, rank1/2
            if (itype(ninter,i) == idx1%itype(j)) then
               tstr(i:i) = idx1%str(j)
               !idx1%str(j) = ''
               !idx1%itype(j) = 0
               nfound = nfound+1
               exit
            end if
         end do
      end do

      do i = rank1/2 + 1, rank1
         do j = rank1/2+1, rank1
            if (itype(ninter,i) == idx1%itype(j)) then
               tstr(i:i) = idx1%str(j)
               !idx1%str(j) = ''
               !idx1%itype(j) = 0
               nfound = nfound+1
               exit
            end if
         end do
      end do

      if (nfound.ne.rank1) then
        write(itflog,'(1x,"error in matching previous intermediate:")')
        write(itflog,'(1x,"label = ",a)') trim(label)
        write(itflog,'(1x,"idx1%itype      = ",12i3)')
     &       idx1%itype(1:rank1)
        write(itflog,'(1x,"itype(ninter,:) = ",12i3)')
     &       itype(ninter,1:rank1)

        call quit(1,'arrange_inter_itype',
     &                'error in matching previous intermediate')
      end if

      do i = 1, rank1
         idx1%str(i) = tstr(i:i)
      end do

      do i = 1, rank1
         idx1%itype(i) = itype(ninter,i)
      end do

      ! update contraction index posistions
      shift = 1
      do i = 1, rank1
         do j = 1, rank2
            if (idx1%str(i) == idx2%str(j)) then
               idx1%cnt_poss(shift) = i
               shift = shift + 1
            end if
         end do
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine arrange_inter_itype2(rank1, rank2, idx1, idx2, itype,
     &                                nop, label)
*----------------------------------------------------------------------*
!     rearrange intermeidate index according to the itype of the
!     previously declared intermediate
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &   rank1,
     &   rank2,
     &   itype(maxint, index_len),
     &   nop
      type(index_str), intent(inout) ::
     &   idx1,
     &   idx2
      character(len=maxlen_bc_label), intent(in) ::
     &   label

      integer ::
     &   i, j,
     &   shift,
     &   ninter
      character (len=index_len) ::
     &   tstr
      logical ::
     &   same

      if (nop==0) then
         write(15,*) "ZERO"
         return
      end if

      if (nop<rank1) then
         write(15,*) "TROUBLE"
         return
      end if

      if (nop==rank1) then
         write(15,*) "GOOD"
         return
      end if

      return

      ! find itype corresponding to numbered intermediate
      i = len(trim(label))
      read(label(i:i),*) ninter

      same = .true.
      do i = 1, rank1
         ! check if intermediate already has correct index itype
         if (idx1%itype(i) /= itype(ninter,i)) then
            same = .false.
         end if
      end do

      if (same) return


      do i = 1, rank1/2
         do j = 1, rank1/2
            if (itype(ninter,i) == idx1%itype(j)) then
               tstr(i:i) = idx1%str(j)
               idx1%str(j) = ''
               idx1%itype(j) = 0
               exit
            end if
         end do
      end do

      do i = rank1/2 + 1, rank1
         do j = rank1/2+1, rank1
            if (itype(ninter,i) == idx1%itype(j)) then
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
         idx1%itype(i) = itype(ninter,i)
      end do

      ! update contraction index posistions
      shift = 1
      do i = 1, rank1
         do j = 1, rank2
            if (idx1%str(i) == idx2%str(j)) then
               idx1%cnt_poss(shift) = i
               shift = shift + 1
            end if
         end do
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

      if (ntest>=100) then
        write(item%out,'(1x,"Top of assign_spin")')
        write(item%out,'(1x,"idx1 = ",a)') item%idx1
        write(item%out,'(1x,"idx2 = ",a)') item%idx2
        write(item%out,'(1x,"idx3 = ",a)') item%idx3
      end if

      ! Assign spin to indicies on the result tensor
      ! If not an intermediate
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
               call quit(1,'assign_spin',
     &             'I think I am not ready for rank 6 results')
               do i=1, 2
                  ! aaaaaa
                  item%t_spin(3)%spin(1,i) = 1
                  item%t_spin(3)%spin(2,i) = 1
               end do
               item%t_spin(3)%spin(1,3) = 1
               item%t_spin(3)%spin(2,3) = 1
            case default
               call line_error("Uncovered case for result tensor rank",
     &                         item)
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
         call print_spin(item%t_spin(1)%spin,item%rank1,item%label_t1,
     &                   item%out)
         call print_spin(item%t_spin(2)%spin,item%rank2,item%label_t2,
     &                   item%out)
         end if

         call spin_index(item, ntest)

      else
         ! Else the result is an intermeidate, spin sum all possible
         ! spin cases

      if (item%rank3>0) then
         ! This first index only runs from 1 and not 2. This is because
         ! abba == baab, so we don't need all the spin cases where beta
         ! is the first spin index. In print_itf_line(), the
         ! intermeidates that have a beta as the first spin undergo a
         ! spin flip (baba -> abab). This reduces the number of intermediates and
         ! contractions, as well as speeding up the python
         ! post-processing
         do n = 1, 1
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
                                 if (item%rank3>6)
     &                                call quit(1,
     &                                'command_to_ift>assign_spin',
     &                                'adapt for rank>6')
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
     &   i, j, k, l, m, n, o, p, q, r, shift,
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
         write(item%out,'(1x,"idx1 = ",a)') item%idx1
         write(item%out,'(1x,"idx2 = ",a)') item%idx2
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
                call save_spin_case(item,eloop,ntest)
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
                  call save_spin_case(item,eloop,ntest)
                 end if
                 if (shift > 6) then
                  do o = 1, 2
                   i1 = poss(1,7)%elements(1)
                   i2 = poss(1,7)%elements(2)
                   i3 = poss(2,7)%elements(1)
                   i4 = poss(2,7)%elements(2)
                   item%t_spin(z1)%spin(i1, i2) = o
                   item%t_spin(z2)%spin(i3, i4) = o
                   if (shift <= 7) then
                    call save_spin_case(item,eloop,ntest)
                   end if
                   if (shift > 7) then
                    do p = 1, 2
                     i1 = poss(1,8)%elements(1)
                     i2 = poss(1,8)%elements(2)
                     i3 = poss(2,8)%elements(1)
                     i4 = poss(2,8)%elements(2)
                     item%t_spin(z1)%spin(i1, i2) = p
                     item%t_spin(z2)%spin(i3, i4) = p
                     if (shift <= 8) then
                      if (item%rank3 == 0 .and. i == 2) exit
                      call save_spin_case(item,eloop,ntest)
                     end if
                     if (shift > 8) then
                      do q = 1, 2
                       i1 = poss(1,9)%elements(1)
                       i2 = poss(1,9)%elements(2)
                       i3 = poss(2,9)%elements(1)
                       i4 = poss(2,9)%elements(2)
                       item%t_spin(z1)%spin(i1, i2) = q
                       item%t_spin(z2)%spin(i3, i4) = q
                       if (shift <= 9) then
                        call save_spin_case(item,eloop,ntest)
                       end if
                       if (shift > 9) then
                        do r = 1, 2
                         i1 = poss(1,10)%elements(1)
                         i2 = poss(1,10)%elements(2)
                         i3 = poss(2,10)%elements(1)
                         i4 = poss(2,10)%elements(2)
                         item%t_spin(z1)%spin(i1, i2) = r
                         item%t_spin(z2)%spin(i3, i4) = r
                         if (shift <= 10) then
                          if (item%rank3 == 0 .and. i == 2) exit
                          call save_spin_case(item,eloop,ntest)
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

            if (item%nspin_cases .gt. MAX_SPIN_CASES)
     &           call quit(1,'save_spin_case','extend max dimensions')

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

      include 'stdunit.h'
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
      subroutine print_symmetrise(item, itf_names, ntest)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(inout) ::
     &     item
      type(tensor_names) ::
     &     itf_names    ! contains renaming information
      integer, intent(in) ::
     &   ntest

      character(len=70) ::
     &   line
      character(len=INDEX_LEN) ::
     &   tindex, new_idx
      character(len=MAXLEN_BC_LABEL) ::
     &     new

      character(len=MAXLEN_BC_LABEL), external ::
     &     rename_tensor

      ! Return if an intermediate or result rank less than 2
      if (item%inter(3) .or. item%rank3<=2) return

      ! Return if not a symmetric result
      if (.not. item%symmetric) return

      ! Return if INTPP block
      if (item%intpp) return

      if (ntest>=100) call debug_header("print_symmetrise", item%out)

      new = rename_tensor(item%old_name,
     &                    item%rank3, item%nops3, itf_names)

      new_idx = item%idx3
      call reorder_to_slots(item%int(3),item%rank3,
     &                      new_idx,
     &                      new,item%nops3,item%out)
      
      line = '.'//trim(new)//'['//trim(new_idx)//'] += '//
     &       trim(item%label_res)//'['//trim(new_idx)//']'
      write(item%out,'(a)') trim(line)

      tindex = ' '
      tindex(1:1) = new_idx(2:2)
      tindex(2:2) = new_idx(1:1)
      tindex(3:3) = new_idx(4:4)
      tindex(4:4) = new_idx(3:3)

      line = '.'//trim(new)//'['//trim(new_idx)//'] += '//
     &       trim(item%label_res)//'['//trimal(tindex)//']'
      write(item%out,'(a)') trim(line)

      return
      end


*----------------------------------------------------------------------*
      subroutine check_symmetric(contr_info, binary, symmetric, nosym,
     &      lu)
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
     &   lu
      logical, intent(in) ::
     &   binary
      logical, intent(inout) ::
     &   symmetric,  ! Check if tensor has R[abij] = R[baji] symmetry
     &   nosym       ! Check if tensor is R[apiq]

      integer ::
     &   c(ngastp,2),
     &   e1(ngastp,2),
     &   e2(ngastp,2),
     &   nops(ngastp,2),
     &   ii, jj, ij, nj

      call itf_ops(contr_info, c, e1, e2, binary)

      nops = e1 + e2
c      nops = sum(e1, dim=2) + sum(e2, dim=2)
c     dbg
c      write(lu,'(1x,"nops:",4i4)') nops 
c     dbg
      
      if (sum(nops)==0) then
         ! Scalars don't have symmetry
         symmetric = .false.
         nosym = .false.
         return
      end if

      !if (sum(nops)>4) then
      !  call quit(1,'command_to_itf','review me for 3-index results')
      !end if


      ! new: do the check based on the actual info on contr_info
      nj = contr_info%nj_res
      symmetric = .true.
      do ij = 1, nj
        do jj = 1, 2
          do ii = 1, ngastp
            if (contr_info%occ_res(ii,jj,ij).eq.1) then
              symmetric = .false.
            end if
          end do
        end do
      end do
      
c      symmetric = .true.
c      do jj = 1, 2
c        do ii = 1, ngastp
c          if (mod(nops(ii,jj),2) /= 0) then
c            symmetric = .false.
c          end if
c        end do
c      end do

      nosym = .false.
c      if (sum(nops)>2 .and. nops(3,1)==1 .and. nops(3,2)==1) then
c         if (.not. check_inter(contr_info%label_res)) then
c            nosym = .true.
c         end if
c     end if
      if (.not.check_inter(contr_info%label_res)) then
        if (sum(nops)==4) then
          nosym = .true.
          do jj = 1, 2
            do ii = 1, ngastp
              if (nops(ii,jj)==2) nosym=.false. 
            end do
          end do 
        end if
        
        if (sum(nops)>4) then
          call quit(1,'command_to_itf->check_symmetric',
     &         'extend me for ops with > 4 indices')
        end if
      end if
        
c     dbg
c      write(lu,*) 'symmetric, nosym: ',symmetric, nosym
c     dbg
      
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
      subroutine itf_ops(contr_info, c, e1, e2, binary)
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

      logical, intent(in) ::
     &   binary

      integer ::
     &   i

      c=0
      e1=0
      e2=0

      if (.not.binary) then
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
      subroutine itf_vertex_ops(contr_info, item, command)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'mdef_operator_info.h' ! For def_formular_item.h
      include 'def_contraction.h'
      include 'def_formula_item.h' ! For command parameters
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &   contr_info     ! Information about binary contraction
      type(itf_contr), intent(inout) ::
     &   item
      integer, intent(in) ::
     &   command
      integer, pointer ::
     &  op1(:,:,:) => null(),
     &  op2(:,:,:) => null(),
     &  res(:,:,:) => null()

      integer ::
     &   i, j, k

      ! Read in for external indicies
      ! Calculate contraction indicies + which vertex

      allocate(item%v1(contr_info%nj_op1,ngastp,2))
      allocate(item%v2(contr_info%nj_op2,ngastp,2))
      allocate(item%v3(contr_info%nj_res,ngastp,2))
      allocate(item%vc1(contr_info%nj_op1,ngastp,2))
      allocate(item%vc2(contr_info%nj_op2,ngastp,2))
      allocate(item%vc3(contr_info%nj_res,ngastp,2))

      allocate(item%vnops1(contr_info%nj_op1))
      allocate(item%vnops2(contr_info%nj_op2))
      allocate(item%vnops3(contr_info%nj_res))

      allocate(op1(contr_info%nj_op1,ngastp,2))
      allocate(op2(contr_info%nj_op2,ngastp,2))

      item%v1 = 0
      item%v2 = 0
      item%v3 = 0
      item%vc1 = 0
      item%vc2 = 0
      item%vc3 = 0

      op1 = 0
      op2 = 0

      if (command==command_cp_intm .or.
     &    command==command_add_intm .or.
     &    command==command_add_reo) then
         do i = 1, contr_info%nj_op1
           call count_index(contr_info%occ_op1(1:,1:,i), item%v1(i,:,:))
           item%vnops1(i) = sum(sum(item%v1(i,:,:),dim=1))
         end do
      else
         ! Get external indicies
         do i = 1, contr_info%nj_op1
           call count_index(contr_info%occ_ex1(1:,1:,i), item%v1(i,:,:))
         end do
         do i = 1, contr_info%nj_op2
           call count_index(contr_info%occ_ex2(1:,1:,i),item%v2(i,:,:))
         end do

         ! Get complete index
         do i = 1, contr_info%nj_op1
           call count_index(contr_info%occ_op1(1:,1:,i), op1(i,:,:))
         end do
         do i = 1, contr_info%nj_op2
           call count_index(contr_info%occ_op2(1:,1:,i), op2(i,:,:))
         end do

         do i = 1, contr_info%nj_op1
            do j = 1, ngastp
               do k = 1, 2
                  item%vc1(i,j,k) = op1(i,j,k) - item%v1(i,j,k)
               end do
            end do
         end do

         do i = 1, contr_info%nj_op2
            do j = 1, ngastp
               do k = 1, 2
                  item%vc2(i,j,k) = op2(i,j,k) - item%v2(i,j,k)
               end do
            end do
         end do

         do i = 1, contr_info%nj_op1
            item%vnops1(i) = sum(sum(item%v1(i,:,:),dim=1)) +
     &                       sum(sum(item%vc1(i,:,:),dim=1))
            !write(item%out,*) "vnops1 ", item%vnops1(i)
         end do
         do i = 1, contr_info%nj_op2
            item%vnops2(i) = sum(sum(item%v2(i,:,:),dim=1)) +
     &                       sum(sum(item%vc2(i,:,:),dim=1))
            !write(item%out,*) "vnops2 ", item%vnops2(i)
          end do

          do i = 1, contr_info%nj_res
            do j = 1, ngastp
              do k = 1, 2
                item%v3(i,j,k) =
     &               contr_info%occ_res(j,k,i)
              end do
            end do
          end do
         !write(item%out,*) "External"
         !do k = 1, contr_info%nj_op1
         !   do j = 1, 2
         !      if (j==1) write(item%out,'(a6)',advance='no') "/"
         !      if (j==2) write(item%out,'(a6)',advance='no') "\"
         !      do i = 1, ngastp
         !         write(item%out,'(i0, 1x)', advance='no')item%v1(k,i,j)
         !      end do
         !      if (j==1) write(item%out,'(a1)',advance='no') "\"
         !      if (j==2) write(item%out,'(a1)',advance='no') "/"
         !      write(item%out,*) "      "
         !   end do
         !   write(item%out,*)
         !end do

         !write(item%out,*) "Contraction"
         !do k = 1, contr_info%nj_op1
         !   do j = 1, 2
         !      if (j==1) write(item%out,'(a6)',advance='no') "/"
         !      if (j==2) write(item%out,'(a6)',advance='no') "\"
         !      do i = 1, ngastp
         !         write(item%out,'(i0, 1x)',advance='no')item%vc1(k,i,j)
         !      end do
         !      if (j==1) write(item%out,'(a1)',advance='no') "\"
         !      if (j==2) write(item%out,'(a1)',advance='no') "/"
         !      write(item%out,*) "      "
         !   end do
         !   write(item%out,*)
         !end do

      end if

      deallocate(op1)
      deallocate(op2)

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
