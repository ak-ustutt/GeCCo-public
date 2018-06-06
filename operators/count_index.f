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

      nops(1,1)=nops(1,1) + iocc(1,1)
      nops(2,1)=nops(2,1) + iocc(2,1)
      nops(3,1)=nops(3,1) + iocc(3,1)
      nops(4,1)=nops(4,1) + iocc(4,1)
      nops(1,2)=nops(1,2) + iocc(1,2)
      nops(2,2)=nops(2,2) + iocc(2,2)
      nops(3,2)=nops(3,2) + iocc(3,2)
      nops(4,2)=nops(4,2) + iocc(4,2)

      return
      end

*----------------------------------------------------------------------*
      subroutine remove_whitespace(string)
*----------------------------------------------------------------------*
!     Remove whitespace inbetween string
*----------------------------------------------------------------------*

      implicit none

      character(len=8), intent(inout) ::
     &     string   ! IFT index string with whitespaces
      character(len=8) ::
     &     tmp      ! Tempory holder of index
      integer ::
     &     i,j      ! Loop index

      tmp='    '
      j=0

      do i=1, len(string)
        if(string(i:i).ne.' ') then
          j=j+1
          tmp(j:j)=string(i:i)
        endif
      enddo
 
      string=tmp

      return
      end

*----------------------------------------------------------------------*
      pure function rename_tensor(string)
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
      character(len=maxlen_bc_label) ::
     &    rename_tensor

      if (trim(string).eq.'OMG') then
          rename_tensor='R'
      else
          rename_tensor=trim(string)
      end if

      end function

*----------------------------------------------------------------------*
      subroutine print_itf_line(item,s1,s2)
*----------------------------------------------------------------------*
!     Print line of ITF code
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(itf_contr), intent(in) ::
     &     item
      logical, intent(in) ::
     &     s1,s2

      character(len=maxlen_bc_label) ::
     &     nres, nt1, nt2,          ! Name of tensors involved in the contraction
     &     ires, it1, it2,          ! Name of tensors + index
     &     sres, st1, st2           ! Name of spin summed tensors + index
      character(len=5) ::
     &     s_int                 ! Intermdiate tensor number
      character(len=50) ::
     &     itf_line,             ! Line of ITF code
     &     d_line,             ! Line of ITF code
     &     l_line,             ! Line of ITF code
     &     a_line,             ! Line of ITF code
     &     s_line              ! Line of ITF code
      character(len=25) ::
     &     sfact='                         ',             ! String representation of factor
     &     sfact_star='                         '              ! String representation of factor formatted for output
      character(len=maxlen_bc_label) ::
     &     rename_tensor

      ! Remove these...
      nres=item%label_res
      nt1=item%label_t1
      nt2=item%label_t2

      ! Change names of specific tensors
      nres=rename_tensor(nres)
      
      ! Spin summ
      if (s1) then
         ! Pure spin
         st1='('//trim(adjustl(nt1))//'['//trim(item%idx1)//']'//
     &       '-'//trim(adjustl(nt1))//'['//trim(item%idx1)//']'//')'
      else
         st1=trim(adjustl(nt1))//'['//trim(item%idx1)//']'
      end if

      if (s2) then
         ! Pure spin
         st2='('//trim(adjustl(nt2))//'['//trim(item%idx2)//']'//
     &       '-'//trim(adjustl(nt2))//'['//trim(item%idx2)//']'//')'
      else
         st2=trim(adjustl(nt2))//'['//trim(item%idx2)//']'
      end if

      ! Convert factor to string, ignore if 1.0 or -1.0
      if (item%fact.ne.1.0) then
          if (item%fact.ne.-1.0) then
              write(sfact,*) item%fact
              sfact_star=' '//trim(sfact)//'*'
          end if
      end if

      if (item%fact.lt.0.0) then
          itf_line='.'//trim(adjustl(nres))//'['//trim(item%idx3)//
     &        '] -= '//trim(sfact_star)//
     &        trim(adjustl(st1))//
     &        trim(adjustl(st2))
      else
          itf_line='.'//trim(adjustl(nres))//'['//trim(item%idx3)//
     &        '] += '//trim(sfact_star)//
     &        trim(adjustl(st1))//
     &        trim(adjustl(st2))
      end if

      write(item%logfile,*) trim(itf_line)

      return
      end

*----------------------------------------------------------------------*
      subroutine assign_index(contr_info,item)
*----------------------------------------------------------------------*
!     
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'
      
      type(binary_contr), intent(in) ::
     &     contr_info   ! Information about binary contraction
      type(itf_contr), intent(inout) ::
     &     item
      integer ::
     &     t1(4,2),      ! Occupations of operator 1
     &     c(4,2),       ! Occupations of contraction index
     &     e1(4,2),      ! Occupations of external index 1
     &     e2(4,2)       ! Occupations of external index 2
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
     &     t1_array,    ! Operator 1 array
     &     c_array,     ! Contraction index array
     &     e1_array,    ! External index of operator 1 array
     &     e2_array     ! External index of operator 2 array
      integer, dimension(3) ::
     &     idx_type     ! Info about index convention
      integer, parameter ::
     &     t_amp = 0,       ! [apij] (aacc)
     &     ham   = 1        ! [abip]
      character(len=20), dimension(4) ::
     &     tensor_ham=(/ 'H', 'INT_D', 'INT_HH', 'INT_PP' /)  ! Tensor names to use ham index convention


      !istr1='        '
      !istr2='        '
      !istr3='        '

      ! Set index convention
      idx_type=(/ 0, 0, 0 /)
      do i=1, len(tensor_ham)
          if (contr_info%label_op1.eq.trim(tensor_ham(i))) then
              ! Use default convention for now
              idx_type(1)=0
          end if
          if (contr_info%label_op2.eq.trim(tensor_ham(i))) then
              idx_type(2)=0
          end if
          if (contr_info%label_res.eq.trim(tensor_ham(i))) then
              idx_type(3)=0
          end if
      end do

      t1=0
      c=0
      e1=0
      e2=0

      ! Get occuation info
      do i = 1, contr_info%nj_op1
        call count_index(i,
     &     contr_info%occ_op1(1:,1:,i),
     &     contr_info%rst_op1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,t1)
      end do
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

      ! Assign t1 (tensor 1)
      ! Creations first
      do i=1, t1(2,1)
          c1(i:)=par(i)
          ! Below doesn't work...
          !t1_array(1)(i:)=par(i)
      end do
      t1_array(1)=c1
      do i=1, t1(3,1)
          c2(i:)=val(i)
      end do
      t1_array(2)=c2
      do i=1, t1(1,1)
          c3(i:)=hol(i)
      end do
      t1_array(3)=c3

      ! Annhilators second
      do i=1, t1(2,2)
          a1(i:)=par(i+t1(2,1))
      end do
      t1_array(5)=a1
      do i=1, t1(3,2)
          a2(i:)=val(i+t1(3,1))
      end do
      t1_array(6)=a2
      do i=1, t1(1,2)
          a3(i:)=hol(i+t1(1,1))
      end do
      t1_array(7)=a3

      c1='        '
      c2='        '
      c3='        '
      a1='        '
      a2='        '
      a3='        '


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

      ! Need to to be shifted so to match assignment of t1 above
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
      
      !write(lulog,*) e1(1,1)
      !write(lulog,*) e1(2,1)
      !write(lulog,*) e1(3,1)
      !write(lulog,*) e1(1,2)
      !write(lulog,*) e1(2,2)
      !write(lulog,*) e1(3,2)
      
      ! Assign c (contracted by)
      ! These need to be shifted, so as not to match e1
      do i=1, c(2,1)
          c1(i:)=par(i+e1(2,1)+e1(2,2))
      end do
      c_array(1)=c1
      do i=1, c(3,1)
          c2(i:)=val(i+e1(3,1)+e1(3,2))
      end do
      c_array(2)=c2
      do i=1, c(1,1)
          c3(i:)=hol(i+e1(1,1)+e1(1,2))
      end do
      c_array(3)=c3

      do i=1, c(2,2)
          a1(i:)=par(i+c(2,1)+e1(2,1)+e1(2,2))
      end do
      c_array(5)=a1
      do i=1, c(3,2)
          a2(i:)=val(i+c(3,1)+e1(3,1)+e1(3,2))
      end do
      c_array(6)=a2
      do i=1, c(1,2)
          a3(i:)=hol(i+c(1,1)+e1(1,1)+e1(1,2))
      end do
      c_array(7)=a3

      c1='        '
      c2='        '
      c3='        '
      a1='        '
      a2='        '
      a3='        '

!      do i=1, c(2,2)
!          a1(i:)=par(i+t1(2,1))
!      end do
!      c_array(5)=a1
!      do i=1, c(3,2)
!          a2(i:)=val(i+t1(3,1))
!      end do
!      c_array(6)=a2
!      do i=1, c(1,2)
!          a3(i:)=hol(i+t1(1,1))
!      end do
!      c_array(7)=a3


      ! Assign e2 (external indicies of t2)
      ! Needs to be shifted so it's different from e1
      ! Needs to be shifted so it's different from c
      do i=1, e2(2,1)
          c1(i:)=par(i+e1(2,1)+e1(2,2)+c(2,1)+c(2,2))
      end do
      e2_array(1)=c1
      do i=1, e2(3,1)
          c2(i:)=val(i+e1(3,1)+e1(3,2)+c(3,1)+c(3,2))
      end do
      e2_array(2)=c2
      do i=1, e2(1,1)
          c3(i:)=hol(i+e1(1,1)+e1(1,2)+c(1,1)+c(1,2))
      end do
      e2_array(3)=c3

      do i=1, e2(2,2)
          a1(i:)=par(i+e2(2,1)+e1(2,1)+e1(2,2)+c(2,1)+c(2,2))
      end do
      e2_array(5)=a1
      do i=1, e2(3,2)
          a2(i:)=val(i+e2(3,1)+e1(3,1)+e1(3,2)+c(3,1)+c(3,2))
      end do
      e2_array(6)=a2
      do i=1, e2(1,2)
          a3(i:)=hol(i+e2(1,1)+e1(1,1)+e1(1,2)+c(1,1)+c(1,2))
      end do
      e2_array(7)=a3
      
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
      item%idx1=trim(adjustl(t1_array(1)))//trim(adjustl(t1_array(5)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(7)))//
     &          trim(adjustl(t1_array(2)))//trim(adjustl(t1_array(6)))
      case default
      ! [apij] (aacc), ie. T[abij]
      item%idx1=trim(adjustl(t1_array(1)))//trim(adjustl(t1_array(2)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(5)))//
     &          trim(adjustl(t1_array(6)))//trim(adjustl(t1_array(7)))
      end select

      ! Operator 2
      ! c_array annhilations correspond to t2 creations and vice versa
      select case(idx_type(2))
      case(ham)
       item%idx2=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(5)))//
     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(3)))//
     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(7)))//
     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(6)))
      case default
      item%idx2=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(5)))//
     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(6)))//
     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(7)))//
     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(3)))
      end select

!      write(item%logfile,*) "Permute: ", item%permute
!      select case(item%permute)
!      case(1)
!      ! [apij] (aacc), ie. T[abij]
!      item%idx1=trim(adjustl(t1_array(1)))//trim(adjustl(t1_array(2)))//
!     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(5)))//
!     &          trim(adjustl(t1_array(6)))//trim(adjustl(t1_array(7)))
!
!      item%idx2=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(5)))//
!     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(6)))//
!     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(7)))//
!     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(1)))//
!     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(2)))//
!     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(3)))
!      case(2)
!      ! Permute creations
!       item%idx1=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(1)))//
!     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(2)))//
!     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(3)))//
!     &          trim(adjustl(e1_array(5)))//trim(adjustl(c_array(5)))//
!     &          trim(adjustl(e1_array(6)))//trim(adjustl(c_array(6)))//
!     &          trim(adjustl(e1_array(7)))//trim(adjustl(c_array(7)))
!
!      item%idx2=trim(adjustl(e1_array(1)))//trim(adjustl(c_array(5)))//
!     &          trim(adjustl(e1_array(2)))//trim(adjustl(c_array(6)))//
!     &          trim(adjustl(e1_array(3)))//trim(adjustl(c_array(7)))//
!     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(1)))//
!     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(2)))//
!     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(3)))
!      case(3)
!       item%idx1=trim(adjustl(e1_array(1)))//trim(adjustl(c_array(1)))//
!     &          trim(adjustl(e1_array(2)))//trim(adjustl(c_array(2)))//
!     &          trim(adjustl(e1_array(3)))//trim(adjustl(c_array(3)))//
!     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(5)))//
!     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(6)))//
!     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(7)))
!
!      item%idx2=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(5)))//
!     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(6)))//
!     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(7)))//
!     &          trim(adjustl(e1_array(5)))//trim(adjustl(c_array(1)))//
!     &          trim(adjustl(e1_array(6)))//trim(adjustl(c_array(2)))//
!     &          trim(adjustl(e1_array(7)))//trim(adjustl(c_array(3)))
!      case(4)
!       item%idx1=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(1)))//
!     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(2)))//
!     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(3)))//
!     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(5)))//
!     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(6)))//
!     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(7)))
!
!      item%idx2=trim(adjustl(e1_array(1)))//trim(adjustl(c_array(5)))//
!     &          trim(adjustl(e1_array(2)))//trim(adjustl(c_array(6)))//
!     &          trim(adjustl(e1_array(3)))//trim(adjustl(c_array(7)))//
!     &          trim(adjustl(e1_array(5)))//trim(adjustl(c_array(1)))//
!     &          trim(adjustl(e1_array(6)))//trim(adjustl(c_array(2)))//
!     &          trim(adjustl(e1_array(7)))//trim(adjustl(c_array(3)))
!      case default
!      write(item%logfile,*) "Error, could not determine permutation"
!      item%idx1=trim(adjustl(t1_array(1)))//trim(adjustl(t1_array(2)))//
!     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(5)))//
!     &          trim(adjustl(t1_array(6)))//trim(adjustl(t1_array(7)))
!       item%idx2=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(1)))//
!     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(5)))//
!     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(3)))//
!     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(7)))//
!     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(2)))//
!     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(6)))
!      end select

      

      ! Result
      select case(idx_type(3))
      case(ham)
      item%idx3=trim(adjustl(e1_array(1)))//trim(adjustl(e2_array(1)))//
     &          trim(adjustl(e1_array(5)))//trim(adjustl(e2_array(5)))//
     &          trim(adjustl(e1_array(3)))//trim(adjustl(e2_array(3)))//
     &          trim(adjustl(e1_array(7)))//trim(adjustl(e2_array(7)))//
     &          trim(adjustl(e1_array(2)))//trim(adjustl(e2_array(2)))//
     &          trim(adjustl(e1_array(6)))//trim(adjustl(e2_array(6)))
      case default
      item%idx3=trim(adjustl(e1_array(1)))//trim(adjustl(e2_array(1)))//
     &          trim(adjustl(e1_array(2)))//trim(adjustl(e2_array(2)))//
     &          trim(adjustl(e1_array(3)))//trim(adjustl(e2_array(3)))//
     &          trim(adjustl(e1_array(5)))//trim(adjustl(e2_array(5)))//
     &          trim(adjustl(e1_array(6)))//trim(adjustl(e2_array(6)))//
     &          trim(adjustl(e1_array(7)))//trim(adjustl(e2_array(7)))
      end select

      !call permute_tensors(e1,e2,c,item)

      return
      end

*----------------------------------------------------------------------*
      subroutine permute_tensors(e1,e2,c,item)
*----------------------------------------------------------------------*
!     
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      integer, intent(in) ::
     &     e1(4,2),      ! Occupations of external index 1
     &     e2(4,2),      ! Occupations of external index 2
     &     c(4,2)
      type(itf_contr), intent(inout) ::
     &     item
      integer ::
     &     i,
     &     sum_c1,sum_c2,sum_a1,sum_a2
      logical ::
     &     logi ! Debug delete

      ! Check if not antisym over different verticies

      ! For rank 4. Rank 2,6 and 0 don't need antisymetrising
      if (len(trim(item%idx3))==2 .or. len(trim(item%idx3))==6
     &    .or. len(trim(item%idx3))==0) then
         return
      end if

      ! Don't care about tensor products now
      if (len(trim(item%idx3))==4 .and. len(trim(item%idx1))==2 .and.
     &    len(trim(item%idx2))==2) then
         return
      end if

      sum_c1=0
      sum_c2=0
      sum_a1=0
      sum_a2=0

      do i=1, 4
         sum_c1=sum_c1+e1(i,1)
         sum_c2=sum_c2+e2(i,1)
         sum_a1=sum_a1+e1(i,2)
         sum_a2=sum_a2+e2(i,2)
      end do

      if (sum_c1/=2 .and. sum_c2/=2) then
         if (sum_c1+sum_c2==2) then
            write(item%logfile,*) "permute creations! 0.5*(1-P)"
         end if
      end if

      if (sum_a1/=2 .and. sum_a2/=2) then
         if (sum_a1+sum_a2==2) then
            write(item%logfile,*) "permute annhilations! 0.5*(1-P)"
         end if
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine permute_tensors2(contr_info,perm_array,lulog)
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

      if ((sum(sum(e1,dim=1))+sum(sum(e2,dim=1)))==2) then
         return
      else if ((sum(sum(e1,dim=1))+sum(sum(e2,dim=1)))==0) then
         return
      else if ((sum(sum(e1,dim=1))+sum(sum(e2,dim=1)))==6) then
         return
      end if

      if (e1(1,1)+e2(1,1)==2 .and. e1(3,2)+e2(3,2)==2 .or.
     &    e1(2,1)+e2(2,1)==2 .and. e1(3,2)+e2(3,2)==2 .or.
     &    e1(2,1)+e2(2,1)==2 .and. e1(1,2)+e2(1,2)==2) then
         
         write(lulog,*) "hello ", e1(1,1)+e2(1,1), e1(3,2)+e2(3,2)
         write(lulog,*) "hello ", e1(2,1)+e2(2,1), e1(3,2)+e2(3,2)
         write(lulog,*) "hello ", e1(2,1)+e2(2,1), e1(1,2)+e2(1,2)

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

         shift=1
         ! If sum==2, then both indicies come from same operator, therefore
         ! it doesn't need anti-symm
         if (sum_c1/=2 .and. sum_c2/=2) then
            if (sum_c1+sum_c2==2) then
               write(lulog,*) "permute creations! 0.5*(1-P)"
               perm_array(shift)=1
               shift=shift+1
               perm_array(shift)=2
               shift=shift+1
            end if
         end if

         if (sum_a1/=2 .and. sum_a2/=2) then
            if (sum_a1+sum_a2==2) then
               write(lulog,*) "permute annhilations! 0.5*(1-P)"
               if (shift>2) then
                  ! (1-P)*(1-P)
                  perm_array(shift)=3
                  shift=shift+1
                  perm_array(shift)=4
                  shift=shift+1
               else
                  perm_array(shift)=1
                  shift=shift+1
                  perm_array(shift)=3
                  shift=shift+1
               end if
            end if
         end if

         write(lulog,*) "perm_array ", perm_array

      else
         return
      end if

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
     &     zero_a,zero_b
      logical ::
     &     logi ! Debug delete

      ! Don't care about tensor products now
      if (len(trim(item%idx3))==4 .and. len(trim(item%idx1))==2 .and.
     &    len(trim(item%idx2))==2) then
         call print_itf_line(item,.false.,.false.)
         return
      else if (item%inter(3)) then
         ! Don't spin sum intermediate results
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

      ! Check the spin case
      if (len(trim(item%idx1))==4 .and. len(trim(item%idx2))==4) then
         t1a(1)=item%idx1(1:1)
         t1a(2)=item%idx1(2:2)
         t1b(1)=item%idx1(3:3)
         t1b(2)=item%idx1(4:4)

         t2a(1)=item%idx2(1:1)
         t2a(2)=item%idx2(2:2)
         t2b(1)=item%idx2(3:3)
         t2b(2)=item%idx2(4:4)
      else if (len(trim(item%idx1))==4 .and. 
     &         len(trim(item%idx2))==2) then
         ! Takes care of rank 4 and rank 2 result
         t1a(1)=item%idx1(1:1)
         t1a(2)=item%idx1(2:2)
         t1b(1)=item%idx1(3:3)
         t1b(2)=item%idx1(4:4)

         t2a(1)=item%idx2(1:1)
         t2b(1)=item%idx2(2:2)
      else if (len(trim(item%idx1))==2 .and.
     &         len(trim(item%idx2))==4) then
         ! Swap t1 and t2 if t1=rank 2. Only going to sum over
         ! contraction indices of t1.
         ! Takes care of rank 4 and rank 2 result
         t2a(1)=item%idx1(1:1)
         t2b(1)=item%idx1(2:2)

         t1a(1)=item%idx2(1:1)
         t1a(2)=item%idx2(2:2)
         t1b(1)=item%idx2(3:3)
         t1b(2)=item%idx2(4:4)
      else if (len(trim(item%idx1))==4 .and.
     &         len(trim(item%idx2))==0) then
         t1a(1)=item%idx1(1:1)
         t1a(2)=item%idx1(2:2)
         t1b(1)=item%idx1(3:3)
         t1b(2)=item%idx1(4:4)
      else if (len(trim(item%idx1))==0 .and.
     &         len(trim(item%idx2))==4) then
         ! Swap t1 and t2 if t1=rank 0
         t1a(1)=item%idx2(1:1)
         t1a(2)=item%idx2(2:2)
         t1b(1)=item%idx2(3:3)
         t1b(2)=item%idx2(4:4)
      else if (len(trim(item%idx1))==2 .and. len(trim(item%idx2))==2 
     &         .and. len(trim(item%idx3))==2) then
         ! Don't need to spin sum
         call print_itf_line(item,.false.,.false.)
         return
      else if (len(trim(item%idx1))==0 .and. len(trim(item%idx2))==2
     &         .and. len(trim(item%idx3))==2) then
         call print_itf_line(item,.false.,.false.)
         return
      else if (len(trim(item%idx1))==2 .and. len(trim(item%idx2))==0 
     &         .and. len(trim(item%idx3))==2) then
         call print_itf_line(item,.false.,.false.)
         return
      else if (len(trim(item%idx1))==2 .and. len(trim(item%idx2))==2
     &         .and. len(trim(item%idx3))==0) then
         call print_itf_line(item,.false.,.false.)
         return
      else if (len(trim(item%idx1))==0 .and. len(trim(item%idx2))==0
     &         .and. len(trim(item%idx3))==0) then
         call print_itf_line(item,.false.,.false.)
         return
      else if (len(trim(item%idx1))==6 .or.
     &         len(trim(item%idx2))==4) then
         call print_itf_line(item,.false.,.false.)
         return
      else
         write(item%logfile,*) "Error, can't determine spin case!"
         call print_itf_line(item,.false.,.false.)
         return
      end if

      if (len(trim(item%idx3))==4) then
         r1a(1)=item%idx3(1:1)
         r1a(2)=item%idx3(2:2)
         r1b(1)=item%idx3(3:3)
         r1b(2)=item%idx3(4:4)
      else if (len(trim(item%idx3))==2) then
         r1a(1)=item%idx3(1:1)
         r1b(1)=item%idx3(2:2)
      end if

      ! Working in index groups, set abab (1212) index to individual tensor
      ! index groups.
      do i=1, 2
         do j=1, 2
            if (r1a(j)==t1a(i) .and. r1a(j)/='') then
               s1a(i)=j
            else if (r1a(j)==t1b(i) .and. r1a(j)/='') then
               s1b(i)=j
            end if
    
            if (r1b(j)==t1a(i) .and. r1b(j)/='') then
               s1a(i)=j
            else if (r1b(j)==t1b(i) .and. r1b(j)/='') then
               s1b(i)=j
            end if
         end do
      end do

      do i=1, 2
         do j=1, 2
            if (r1a(j)==t2a(i) .and. r1a(j)/='') then
               s2a(i)=j
            else if (r1a(j)==t2b(i) .and. r1a(j)/='') then
               s2b(i)=j
            end if
    
            if (r1b(j)==t2a(i) .and. r1b(j)/='') then
               s2a(i)=j
            else if (r1b(j)==t2b(i) .and. r1b(j)/='') then
               s2b(i)=j
            end if
         end do
      end do

      ! Index assigned from result
      write(item%logfile,*) "s1b: ", s1b
      write(item%logfile,*) "s1a: ", s1a
      write(item%logfile,*)
      write(item%logfile,*) "s2b: ", s2b
      write(item%logfile,*) "s2a: ", s2a
      write(item%logfile,*)
      write(item%logfile,*)

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
      type(itf_contr) ::
     &     item

      integer ::
     &     i,j,k,l,m,n,o,p,q,r,s,
     &     second_idx,third_idx,fourth_idx

      second_idx=0
      third_idx=0
      fourth_idx=0
      ! Check for first unassigned index
      do i=1, size(s1a)
       if (s1a(i)==0) then
        do j=1, 2
         s1a(i)=j
         call find_idx(s2a,s2b,t1a,t2a,t2b,i,j)

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
     &                                    s2b,logi,item)
                    end do
                   end if
                  end do
                 else if (.not. any(s1b==0) .and. fourth_idx==0) then
                  ! Print result, sum over three indicies
                  call print_spin_case(s1b,s1a,s2a,s2b,logi,
     &                                 item)
                 end if

                end do
               end if
              end do
             else if (.not. any(s1b==0) .and. third_idx==0) then 
              ! Print result, sum over two indicies
              call print_spin_case(s1a,s1b,s2a,s2b,.not.logi,
     &                             item)
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
     &                            item)
            end do
           end if
          end do
         else if (.not. any(s1a==0) .and. second_idx==0) then
          ! Print result, sum over one indicies
          call print_spin_case(s1a,s1b,s2a,s2b,.not.logi,
     &                         item)
         end if
        end do ! Loop over a/b for first index
       end if ! Check for first 0 index
      end do

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

      return
      end

*----------------------------------------------------------------------*
      subroutine print_spin_case(s1a,s1b,s2a,s2b,logi,item)
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
      type(itf_contr), intent(inout) ::
     &     item
      logical ::
     &     s1,       ! True if tensor 1 is mixed spin
     &     s2        ! True if tensor 2 is mixed spin

       s1=.true.
       s2=.true.
       ! Pick out specific spin cases here
       if (sum(s1a)==sum(s1b) .and.
     &     sum(s2a)==sum(s2b)) then
         if (modulo(sum(s1a)+sum(s1b),2)==0 .and.
     &       modulo(sum(s2a)+sum(s2b),2)==0) then

            ! Doesn't work for rank 6 tensors yet...
            if (item%rank1==2 .or. item%rank1==0) then
               s1=.false.
            else if (s1a(1)/=s1a(2)) then
               s1=.false.
            else if (s1a(1)==s1a(2)) then
               ! Pure spin
               s1=.true.
            end if

            if (item%rank2==2 .or. item%rank2==0) then
               s2=.false.
            else if (s2a(1)/=s2a(2)) then
               s2=.false.
            else if (s2a(1)==s2a(2)) then
               s2=.true.
            end if

            if (logi) then
               write(item%logfile,*) "s1b: ", s1b
               write(item%logfile,*) "s1a: ", s1a
               write(item%logfile,*)
               write(item%logfile,*) "s2b: ", s2b
               write(item%logfile,*) "s2a: ", s2a
               write(item%logfile,*)
               write(item%logfile,*)
            else
               write(item%logfile,*) "s1b: ", s1b
               write(item%logfile,*) "s1a: ", s1a
               write(item%logfile,*)
               write(item%logfile,*) "s2b: ", s2b
               write(item%logfile,*) "s2a: ", s2a
               write(item%logfile,*)
               write(item%logfile,*)
            end if

            call print_itf_line(item,s1,s2)
         end if
       end if

       return
       end

*----------------------------------------------------------------------*
      subroutine find_idx(s2a,s2b,t1,t2a,t2b,idx,j)
*----------------------------------------------------------------------*
!    Find corresponding index in second tenor 
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
      subroutine itf_contr_init(contr_info,itf_item,perm,lulog)
*----------------------------------------------------------------------*
!    Initalise 
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'
      
      type(binary_contr), intent(in) ::
     &     contr_info   ! Inofrmation about binary contraction
      type(itf_contr), intent(inout) ::
     &     itf_item
      integer, intent(in) ::
     &     perm,
     &     lulog
      ! This should be in a 'constructor'

      ! Assign output file
      itf_item%logfile=lulog

      ! Assign permutation number
      if (perm==0) then
         itf_item%permute=1
      else
         itf_item%permute=perm
      end if

      ! Assign labels
      itf_item%label_t1=contr_info%label_op1
      itf_item%label_t2=contr_info%label_op2
      itf_item%label_res=contr_info%label_res

      ! Assign rank
      call itf_rank(itf_item%idx1,itf_item%rank1)
      call itf_rank(itf_item%idx2,itf_item%rank2)
      call itf_rank(itf_item%idx3,itf_item%rank3)

      ! Assign factor
      itf_item%fact=contr_info%fact

      ! Check if an intermediate
      call check_inter(itf_item%label_t1,itf_item%inter(1))
      call check_inter(itf_item%label_t2,itf_item%inter(2))
      call check_inter(itf_item%label_res,itf_item%inter(3))

      ! Assign index string
      call assign_index(contr_info,itf_item)

      return
      end

*----------------------------------------------------------------------*
      subroutine itf_rank(idx,nrank)
*----------------------------------------------------------------------*
!    Initalise 
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=index_len), intent(in) ::
     &     idx
      integer, intent(inout) ::
     &     nrank

      nrank=len(trim(adjustl(idx)))

      return
      end

*----------------------------------------------------------------------*
      subroutine check_inter(label,intermediate)
*----------------------------------------------------------------------*
!    Initalise 
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'
      
      character(len=maxlen_bc_label), intent(inout) ::
     &     label
      logical, intent(inout) ::
     &     intermediate

      ! If first character of the label is '_' then the tensor is
      ! defined as an intermedate
      ! The leading '-' is replaced with a whitespace
      if (label(1:1).eq.'_') then
         label(1:1)=' '
         intermediate=.true.
      else
         intermediate=.false.
      end if

      return
      end
