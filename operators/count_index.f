*----------------------------------------------------------------------*
      module itf_utils
*----------------------------------------------------------------------*
!     Contins functions used throughout the code
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
      else if (trim(string).eq.'T2g' .or. trim(string).eq.'T2') then
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
      pure function t_index(index)
*----------------------------------------------------------------------*
!     Transpose ITF index string, abij => abji
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=index_len), intent(in) ::
     &     index       ! ITF index string
      character(len=index_len) ::
     &     t_index      ! Transpose of ITF index string

      character(len=1) ::
     &     tmp

      t_index=index
      tmp=index(1:1)
      t_index(1:1)=index(2:2)
      t_index(2:2)=tmp
      
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

      character(len=100), intent(in) ::
     &     error 
      type(itf_contr), intent(in) ::
     &     item

      write(item%logfile,*) "ERROR: "//error
      write(item%logfile,*) "Result: ", item%label_res, item%idx3
      write(item%logfile,*) "Tensor1: ", item%label_t1, item%idx1
      write(item%logfile,*) "Tensor2: ", item%label_t2, item%idx2

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

      ! Annhilation operators
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

      type(binary_contr), intent(in) ::
     &     contr_info      ! Inofrmation about binary contraction
      integer, intent(in) ::
     &     itflog,         ! Output file
     &     command         ! Type of formula item command, ie. contraction, copy etc.

      type(itf_contr) ::
     &     itf_item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &    perm_array(4),   ! Info of permutation factors
     &    i                ! Loop index
      logical ::
     &    inter            ! True if result is an intermediate

      ! Initalise permutation factors to 0 == no permutation
      perm_array=0
      ! Determine if result needs permuting
      !call check_inter(contr_info%label_res,inter)
      inter = check_inter(contr_info%label_res)
      if (.not.inter) then
         call permute_tensors(contr_info,perm_array,itflog)
      end if

      if (command==command_add_intm .or. command==command_cp_intm) then
         ! For [ADD] and [COPY] cases
         call itf_contr_init(contr_info,itf_item,perm_array(1),
     &                       command,itflog)
         call print_itf_line(itf_item,.false.,.false.)
      else
         ! For other binary contractions
         if (perm_array(1)==0) then
            ! No permutations
            call itf_contr_init(contr_info,itf_item,perm_array(1),
     &                          command,itflog)
            call assign_spin(itf_item)
         else
            do i=1, size(perm_array)
               ! Loop over permuation cases and send seperatley to
               ! assign_spin
               call itf_contr_init(contr_info,itf_item,perm_array(i),
     &                             command,itflog)
               call assign_spin(itf_item)
               if (perm_array(i+1)==0) exit
            end do
         end if
      end if
      
      return
      end

*----------------------------------------------------------------------*
      subroutine intermediate_to_itf(contr_info,itflog,command,
     &                               spin_inters)
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
     &     contr_info      ! Inofrmation about binary contraction
      integer, intent(in) ::
     &     itflog,         ! Output file
     &     command         ! Type of formula item command, ie. contraction, copy etc.
      type(spin_cases), dimension(4) ::
     &     spin_inters

      type(itf_contr) ::
     &     itf_item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &    perm_array(4),   ! Info of permutation factors
     &    i

      perm_array=0

      call itf_contr_init(contr_info,itf_item,perm_array(1),
     &                    command,itflog)

      ! Allocate space to store information about intermediates and
      ! their spin cases. Only allocate 2 objects as there can only be
      ! at most two intermediates on a line
      allocate(itf_item%inter_spins(2))

      ! Do not want to print out the lines while gatheing info about
      ! intermediates
      itf_item%print_line = .false.
      call assign_spin(itf_item)
      itf_item%print_line = .true.

      if (itf_item%ninter == 0) call line_error("Couldn't find
     &                              intermediate", itf_item)

      ! Copy information back to array in print_itf()
      do i = 1, itf_item%ninter
         spin_inters(i)=itf_item%inter_spins(i)
      end do

      deallocate(itf_item%inter_spins)

      return
      end

*----------------------------------------------------------------------*
      subroutine intermediate2_to_itf(contr_info,itflog,command,
     &                               spin_case)
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
     &     contr_info      ! Inofrmation about binary contraction
      integer, intent(in) ::
     &     itflog,         ! Output file
     &     command         ! Type of formula item command, ie. contraction, copy etc.
      integer, intent(in) ::
     &     spin_case(4)

      type(itf_contr) ::
     &     itf_item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &    perm_array(4),   ! Info of permutation factors
     &    i,j
      character(len=4) ::
     &    spin_name

      perm_array=0

      call itf_contr_init(contr_info,itf_item,perm_array(1),
     &                    command,itflog)

      itf_item%spin_case = spin_case

      ! Change intermediate name to reflect spin case
      j = 1
      spin_name = ''
      do i=1, 4
         if (spin_case(i) == 1) then
            spin_name(j:j) = 'a'
            j = j + 1
         else if (spin_case(i) == 2) then
            spin_name(j:j) = 'b'
            j = j + 1
         end if
      end do

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
      ! TODO: is there a better way of avoiding this problem?
      if (item%permute>1) then
         if (item%inter(1)) then
            do i = 1, item%rank1
               if (item%inter1(i:i)=='a') then
                  item%inter1(i:i)='b'
               else if (item%inter1(i:i)=='b') then
                  item%inter1(i:i)='a'
               end if
            end do
         end if

         if (item%inter(2)) then
            do i = 1, item%rank2
               if (item%inter2(i:i)=='a') then
                  item%inter2(i:i)='b'
               else if (item%inter2(i:i)=='b') then
                  item%inter2(i:i)='a'
               end if
            end do
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

      ! Need to include factor if permuting
      if (item%permute>1) then
         if (equal_op=='+=') then
            equal_op='-='
         else if (equal_op=='-=') then
            equal_op='+='
         end if
      end if

       ! Construct complete itf algo line from the above parts
       itf_line='.'//trimal(nres)//
     &     '['//trim(item%idx3)//'] '//equal_op//' '//
     &     trim(sfact_star)//trimal(st1)//' '//trimal(st2)

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
      subroutine assign_index(contr_info,item)
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
     &     c(4,2),       ! Occupations of contraction index
     &     e1(4,2),      ! Occupations of external index 1
     &     e2(4,2)      ! Occupations of external index 2
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
      ! x_array(1:4) = creation operators (par/val/hol/f12)
      ! x_array(5:8) = annhilation operators (par/val/hol/f12)
     &     c_array,     ! Contraction index array
     &     e1_array,    ! External index of operator 1 array
     &     e2_array,    ! External index of operator 2 array
     &     t1_array,
     &     t2_array
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
              idx_type(1)=1
          end if
          if (contr_info%label_op2.eq.trim(tensor_ham(i))) then
              idx_type(2)=1
          end if
          if (contr_info%label_res.eq.trim(tensor_ham(i))) then
              idx_type(3)=1
          end if
      end do

      c=0
      e1=0
      e2=0

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
      

      ! Order in ITF usually follows: apij
      ! Defualt [ccaa] as in the case of T[abij]

      ! Assign e1 (external indicies of t1)
      do i=1, e1(2,1)
          c1(i:)=par(i)
      end do
      e1_array(1)=c1
      do i=1, e1(3,1)
          c2(i:)=val(i)
      end do
      e1_array(2)=c2
      do i=1, e1(1,1)
          c3(i:)=hol(i)
      end do
      e1_array(3)=c3

      ! Need to to be shifted to not match assignment of creations above
      do i=1, e1(2,2)
          a1(i:)=par(i+e1(2,1))
      end do
      e1_array(5)=a1
      do i=1, e1(3,2)
          a2(i:)=val(i+e1(3,1))
      end do
      e1_array(6)=a2
      do i=1, e1(1,2)
          a3(i:)=hol(i+e1(1,1))
      end do
      e1_array(7)=a3

      c1='        '
      c2='        '
      c3='        '
      a1='        '
      a2='        '
      a3='        '

      ! Shifted so as not to match e1 index
      do i=1, e2(2,1)
          c1(i:)=par(i+e1(2,1)+e1(2,2))
      end do
      e2_array(1)=c1
      do i=1, e2(3,1)
          c2(i:)=val(i+e1(3,1)+e1(3,2))
      end do
      e2_array(2)=c2
      do i=1, e2(1,1)
          c3(i:)=hol(i+e1(1,1)+e1(1,2))
      end do
      e2_array(3)=c3

      ! Shifted so as not to match e1 index or above creations
      do i=1, e2(2,2)
          a1(i:)=par(i+e2(2,1)+e1(2,1)+e1(2,2))
      end do
      e2_array(5)=a1
      do i=1, e2(3,2)
          a2(i:)=val(i+e2(3,1)+e1(3,1)+e1(3,2))
      end do
      e2_array(6)=a2
      do i=1, e2(1,2)
          a3(i:)=hol(i+e2(1,1)+e1(1,1)+e1(1,2))
      end do
      e2_array(7)=a3
      
      c1='        '
      c2='        '
      c3='        '
      a1='        '
      a2='        '
      a3='        '

!      ! Permute indicies to get antisymm tensors
!      if (item%permute==0) then
!         ! No permutations
!         do i=1, size(t1_array)
!            t1_array(i)=e1_array(i)
!            t2_array(i)=e2_array(i)
!         end do
!      else if (item%permute==1) then
!         ! No permutations, but get a factor
!         do i=1, size(t1_array)
!            t1_array(i)=e1_array(i)
!            t2_array(i)=e2_array(i)
!         end do
!      else if (item%permute==2) then
!         ! Permute creations
!         do i=1, size(t1_array)/2
!            t1_array(i)=e2_array(i)
!            t1_array(i+4)=e1_array(i+4)
!            t2_array(i)=e1_array(i)
!            t2_array(i+4)=e2_array(i+4)
!         end do
!      else if (item%permute==3) then
!         ! Permute annhilations
!         do i=1, size(t1_array)/2
!            t1_array(i)=e1_array(i)
!            t1_array(i+4)=e2_array(i+4)
!            t2_array(i)=e2_array(i)
!            t2_array(i+4)=e1_array(i+4)
!         end do
!      else if (item%permute==4) then
!         ! Permute creations and annhilations
!         do i=1, size(t1_array)
!            t1_array(i)=e2_array(i)
!            t2_array(i)=e1_array(i)
!         end do
!      end if

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
            t1_array(i)=e2_array(i)
            t1_array(i+4)=e1_array(i+4)
            t2_array(i)=e1_array(i)
            t2_array(i+4)=e2_array(i+4)
         end do
      else if (item%permute==3) then
         ! Permute annhilations
         do i=1, size(t1_array)/2
            t1_array(i)=e1_array(i)
            t1_array(i+4)=e2_array(i+4)
            t2_array(i)=e2_array(i)
            t2_array(i+4)=e1_array(i+4)
         end do
      else if (item%permute==4) then
         ! Permute creations and annhilations
         do i=1, size(t1_array)
            t1_array(i)=e2_array(i)
            t2_array(i)=e1_array(i)
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
      item%idx1=trimal(t1_array(1))//trimal(c_array(1))//
     &          trimal(t1_array(5))//trimal(c_array(5))//
     &          trimal(t1_array(3))//trimal(c_array(3))//
     &          trimal(t1_array(7))//trimal(c_array(7))//
     &          trimal(t1_array(2))//trimal(c_array(2))//
     &          trimal(t1_array(6))//trimal(c_array(6))
      case default
      ! [apij] (aacc), ie. T[abij]
      item%idx1=trimal(t1_array(1))//trimal(c_array(1))//
     &          trimal(t1_array(2))//trimal(c_array(2))//
     &          trimal(t1_array(3))//trimal(c_array(3))//
     &          trimal(t1_array(5))//trimal(c_array(5))//
     &          trimal(t1_array(6))//trimal(c_array(6))//
     &          trimal(t1_array(7))//trimal(c_array(7))
      end select

      ! Operator 2
      ! c_array annhilations correspond to t2 creations and vice versa
      select case(idx_type(2))
      case(ham)
      item%idx2=trimal(t1_array(1))//trimal(c_array(1))//
     &          trimal(t1_array(5))//trimal(c_array(5))//
     &          trimal(t1_array(3))//trimal(c_array(3))//
     &          trimal(t1_array(7))//trimal(c_array(7))//
     &          trimal(t1_array(2))//trimal(c_array(2))//
     &          trimal(t1_array(6))//trimal(c_array(6))
      case default
      item%idx2=trimal(t2_array(1))//trimal(c_array(5))//
     &          trimal(t2_array(2))//trimal(c_array(6))//
     &          trimal(t2_array(3))//trimal(c_array(7))//
     &          trimal(t2_array(5))//trimal(c_array(1))//
     &          trimal(t2_array(6))//trimal(c_array(2))//
     &          trimal(t2_array(7))//trimal(c_array(3))
      end select

      ! Result
      select case(idx_type(3))
      case(ham)
      item%idx3=trimal(t1_array(1))//trimal(c_array(1))//
     &          trimal(t1_array(5))//trimal(c_array(5))//
     &          trimal(t1_array(3))//trimal(c_array(3))//
     &          trimal(t1_array(7))//trimal(c_array(7))//
     &          trimal(t1_array(2))//trimal(c_array(2))//
     &          trimal(t1_array(6))//trimal(c_array(6))
      case default
      item%idx3=trimal(e1_array(1))//trimal(e2_array(1))//
     &          trimal(e1_array(2))//trimal(e2_array(2))//
     &          trimal(e1_array(3))//trimal(e2_array(3))//
     &          trimal(e1_array(5))//trimal(e2_array(5))//
     &          trimal(e1_array(6))//trimal(e2_array(6))//
     &          trimal(e1_array(7))//trimal(e2_array(7))
      end select

!      call index_convention(t1_array,c_array,idx_type(1),
!     &                      item%idx1,.true.,.true.)
!      call index_convention(t2_array,c_array,idx_type(2),
!     &                      item%idx2,.true.,.false.)
!      call index_convention(e1_array,e2_array,idx_type(3),
!     &                      item%idx3,.true.,.true.)

      return
      end

*----------------------------------------------------------------------*
      subroutine index_convention(arr1,arr2,ttype,item_idx,t1,t2)
*----------------------------------------------------------------------*
!     Arrange index letters into correct order according to convention 
*----------------------------------------------------------------------*

      use itf_utils
      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=8), dimension(8), intent(in) ::
     &     arr1,     ! Contraction index array
     &     arr2
      integer, intent(in) ::
     &     ttype    ! Info about index convention
      logical, intent(in) ::
     &     t1, t2
      character(len=index_len), intent(inout) ::
     &     item_idx

      integer, parameter ::
     &     t_amp = 0,       ! [apij] (aacc)
     &     ham   = 1        ! [abip]
      integer, dimension(6) ::
     &     a1 = (/ 1, 2, 3, 5, 6, 7 /),
     &     a2 = (/ 5, 6, 7, 1, 2, 3 /),
     &     num

      if (t1) then
         num = t1
      else
         num = t2
      end if

      if (t2) then
         num = t1
      else
         num = t2
      end if

      select case(ttype)
      case(ham)
      item_idx=trimal(arr1(num(4)))//trimal(arr2(num(4)))//
     &         trimal(arr1(num(5)))//trimal(arr2(num(5)))//
     &         trimal(arr1(num(6)))//trimal(arr2(num(6)))//
     &         trimal(arr1(num(1)))//trimal(arr2(num(1)))//
     &         trimal(arr1(num(2)))//trimal(arr2(num(2)))//
     &         trimal(arr1(num(3)))//trimal(arr2(num(3)))
      case default
      item_idx=trimal(arr1(num(1)))//trimal(arr2(num(1)))//
     &         trimal(arr1(num(2)))//trimal(arr2(num(2)))//
     &         trimal(arr1(num(3)))//trimal(arr2(num(3)))//
     &         trimal(arr1(num(4)))//trimal(arr2(num(4)))//
     &         trimal(arr1(num(5)))//trimal(arr2(num(5)))//
     &         trimal(arr1(num(6)))//trimal(arr2(num(6)))
      end select

      return
      end


*----------------------------------------------------------------------*
      subroutine permute_tensors(contr_info,perm_array,lulog)
*----------------------------------------------------------------------*
!     
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(binary_contr), intent(in) ::
     &     contr_info   ! Information about binary contraction
      integer, intent(inout) ::
     &     perm_array(4)
      integer, intent(in) ::
     &     lulog

      integer ::
     &     e1(4,2),      ! Occupations of external index 1
     &     e2(4,2),      ! Occupations of external index 2
     &     c(4,2)
      integer ::
     &     i,
     &     sum_c1,sum_c2,sum_a1,sum_a2,
     &     shift

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

      !write(lulog,*) "e1+e2 ", sum(sum(e1,dim=1))+sum(sum(e2,dim=1))
      !write(lulog,*) "e1 ", e1
      !write(lulog,*) "e2 ", e2

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

         !write(lulog,*) "sum_c1: ", sum_c1
         !write(lulog,*) "sum_c2: ", sum_c1
         !write(lulog,*) "sum_a1: ", sum_c1
         !write(lulog,*) "sum_c1: ", sum_c1

         shift=1
         ! If sum==2, then both indicies come from same operator, therefore
         ! it doesn't need anti-symm
         if (sum_c1/=2 .and. sum_c2/=2) then
            if (sum_c1+sum_c2==2) then
               !write(lulog,*) "permute creations!"
               perm_array(shift)=1
               shift=shift+1
               perm_array(shift)=2
               shift=shift+1
            end if
         end if

         if (sum_a1/=2 .and. sum_a2/=2) then
            if (sum_a1+sum_a2==2) then
               !write(lulog,*) "permute annhilations!"
               if (shift>2) then
                  ! P(ab)P(ij)
                  !perm_array(shift)=3
                  !shift=shift+1
                  ! TODO: overwirtie 2 for now - mess!
                  perm_array(shift-1)=4
                  shift=shift+1
               else
                  perm_array(shift)=1
                  shift=shift+1
                  perm_array(shift)=3
                  shift=shift+1
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
!     Asign index to covarient and contravarient groups
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
!    Assign spin to tensors, then sum remaining contraction indicies 
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
        write(item%logfile,*) "END"
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
                  item%inter_spins(ishift)%name=item%label_t1

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
                  item%inter_spins(ishift)%name=item%label_t2

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

!               do i=1, 4
!                  write(item%logfile,*) "SC ",
!     &            item%inter_spins(1)%cases(i,shift)
!               end do

               write(item%logfile,*) "Name: ", item%inter_spins(1)%name
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
!     Initalise ITF contraction object
*----------------------------------------------------------------------*

      use itf_utils
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h' ! For def_formular_item.h
      include 'def_contraction.h'
      include 'def_formula_item.h' ! For command parameters
      include 'def_itf_contr.h'
      
      type(binary_contr), intent(in) ::
     &     contr_info   ! Inofrmation about binary contraction
      type(itf_contr), intent(inout) ::
     &     itf_item     ! Object which holds information necessary to print out an ITF algo line
      integer, intent(in) ::
     &     perm,        ! Permuation information
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

      ! Assign factor
      itf_item%fact=contr_info%fact

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

      return
      end

*----------------------------------------------------------------------*
      subroutine itf_rank(idx,nrank)
*----------------------------------------------------------------------*
!    Initalise 
*----------------------------------------------------------------------*

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

!*----------------------------------------------------------------------*
!      subroutine check_inter(label,intermediate)
!*----------------------------------------------------------------------*
!!    Check if tensor is an intermediate
!*----------------------------------------------------------------------*
!
!      implicit none
!      include 'opdim.h'
!      include 'def_contraction.h'
!      include 'def_itf_contr.h'
!      
!      character(len=maxlen_bc_label), intent(inout) ::
!     &     label
!      logical, intent(inout) ::
!     &     intermediate
!
!      ! If first character of the label is '_' then the tensor is
!      ! defined as an intermedate
!      ! The leading '-' is replaced with a whitespace
!      if (label(1:1).eq.'_') then
!         label(1:1)=' '
!         intermediate=.true.
!      else
!         intermediate=.false.
!      end if
!
!      ! Intermediates are also called like below
!      if (scan(label, "STIN")>0 .or.
!     &    scan(label, "LTIN")>0) then
!         intermediate=.true.
!      else
!         intermediate=.false.
!      end if
!
!      return
!      end
