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

      nops(1,1)=iocc(1,1)
      nops(2,1)=iocc(2,1)
      nops(3,1)=iocc(3,1)
      nops(4,1)=iocc(4,1)
      nops(1,2)=iocc(1,2)
      nops(2,2)=iocc(2,2)
      nops(3,2)=iocc(3,2)
      nops(4,2)=iocc(4,2)

      return
      end

*----------------------------------------------------------------------*
      subroutine array_to_string(ans_arr, t1_arr, t2_arr, idx1, idx2,
     &                           idx3)
*----------------------------------------------------------------------*
!     Copy index array to index string
*----------------------------------------------------------------------*

      implicit none

      character, intent(in) ::
     &     ans_arr(2,2),
     &     t1_arr(2,2),
     &     t2_arr(2,2)
      character(len=4), intent(inout) ::
     &     idx1,
     &     idx2,
     &     idx3
      integer ::
     &     i,j      ! Loop index

      do i=1, 2
        do j=1, 2
          if (t1_arr(j,i).ne.'>') then
            idx1(j+(2*i-2):)=t1_arr(j,i)
          end if
          if (t2_arr(j,i).ne.'>') then
            idx2(j+(2*i-2):)=t2_arr(j,i)
          end if
          if (ans_arr(j,i).ne.'>') then
            idx3(j+(2*i-2):)=ans_arr(j,i)
          end if
        end do
      end do

      call remove_whitespace(idx1)
      call remove_whitespace(idx2)
      call remove_whitespace(idx3)

      return
      end

*----------------------------------------------------------------------*
      subroutine remove_whitespace(string)
*----------------------------------------------------------------------*
!     Remove whitespace inbetween ITF index string
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
      subroutine clear_index(str1, str2, str3)
*----------------------------------------------------------------------*
!     Clear ITF index array and string
*----------------------------------------------------------------------*

      implicit none

      character(len=8), intent(inout) ::
     &     str1, str2, str3

      str1='        '
      str2='        '
      str3='        '

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
      subroutine spin_sum_index(itf1, itf2, itf3, lulog)
*----------------------------------------------------------------------*
!     Spin sum index and produce resulting binary contractions
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_tensor.h'

      type(itf_tensor), intent(in) ::
     &    itf1, itf2, itf3          ! Op1, op2, res
      integer, intent(in) ::
     &    lulog
      integer ::
     &    i

      if (itf3%rank.eq.2) then
          ! Just one case, all alpha == all beta
          do i=1, len(itf1%idx)
          end do
      end if
      
      return
      end

*----------------------------------------------------------------------*
      subroutine construct_tensor(res, op1, op2, istr1, istr2, istr3,
     &                            itf1, itf2, itf3, factor, lulog)
*----------------------------------------------------------------------*
!     Form itf_tensor objects
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_tensor.h'

      character(len=maxlen_bc_label), intent(in) ::
     &     res, op1, op2           ! Name of tensors involved in the contraction
      character(len=index_len), intent(in) ::
     &    istr1, istr2, istr3
      type(itf_tensor), intent(inout) ::
     &    itf1, itf2, itf3         ! Op1, op2, res
      real(8) ::
     &    factor
      integer, intent(in) ::
     &    lulog
      
      itf1%name=op1
      itf2%name=op2
      itf3%name=res
      
      itf1%idx=istr1
      itf2%idx=istr2
      itf3%idx=istr3

      itf1%rank=len(trim(istr1))
      itf2%rank=len(trim(istr2))
      itf3%rank=len(trim(istr3))

      itf1%fact=factor
      itf2%fact=factor
      itf3%fact=1.0

      return
      end

*----------------------------------------------------------------------*
      subroutine recursive_formular_item(fl_item,lulog)
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      type(formula_item), intent(in), target ::
     &     fl_item
      integer, intent(in) ::
     &     lulog
      type(formula_item), pointer ::
     &     next_item,  ! Next formula_item
     &     next2_item  ! Next formula_item

      next_item=>fl_item%next
      write(lulog,*) "WHAT: ", next_item%command
      next2_item=>next_item%next
      write(lulog,*) "WHAT2: ", next2_item%command

      return
      end

*----------------------------------------------------------------------*
      subroutine print_itf_line(res, t1, t2, fact, idx1, idx2,
     &                          idx3, inter, lulog)
*----------------------------------------------------------------------*
!     Print line of ITF code
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      character(len=maxlen_bc_label), intent(in) ::
     &     res, t1, t2           ! Name of tensors involved in the contraction
      character(len=4), intent(in) ::
     &     idx1, idx2, idx3      ! Index strings
      integer, intent(in) ::
     &     inter,                ! Intermediate number involved in contraction
     &     lulog                 ! File to write to
      real(8), intent(in) ::
     &     fact                ! Factor associated with binary contraction

      character(len=maxlen_bc_label) ::
     &     nres, nt1, nt2,          ! Name of tensors involved in the contraction
     &     ires, it1, it2           ! Name of tensors + index
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
      logical ::
     &     interm_res=.false.,    ! True if result involves intermediate
     &     interm_t1=.false.,      ! True if contraction involves intermediate
     &     interm_t2=.false.       ! True if contraction involves intermediate
      character(len=maxlen_bc_label) ::
     &     rename_tensor

      ! Remove leading '_' from intermediate label
      nres=res
      nt1=t1
      nt2=t2

      if (nres(1:1).eq.'_') then
          nres(1:1)=' '
          ! Check if result is an intermediate
          interm_res=.true.
      end if
      if (nt1(1:1).eq.'_') then
          nt1(1:1)=' '
          interm_t1=.true.
      end if
      if (nt2(1:1).eq.'_') then
          nt2(1:1)=' '
          interm_t2=.true.
      end if

      ! Change names of specific tensors
      nres=rename_tensor(nres)

      ires=trim(adjustl(nres))//'['//trim(idx3)//']'
      it1=trim(adjustl(nt1))//'['//trim(idx1)//']'
      it2=trim(adjustl(nt2))//'['//trim(idx2)//']'

      if (interm_res) then
          ! Need to alloc intermediate
          a_line='alloc '//trim(ires)
          l_line='load '//trim(it1)//', '//trim(it2)
          d_line='drop '//trim(it2)//', '//trim(it1)
          s_line='store '//trim(ires)
      else
          if (interm_t1) then
              ! Intermediate already in memory from previous alloc
              ! Just drop it at the end, don't store it
              l_line='load '//trim(it2)
              d_line='drop '//trim(it2)//', '//trim(it1)
          elseif (interm_t2) then
              l_line='load '//trim(it1)
              d_line='drop '//trim(it2)//', '//trim(it1)
          else
              l_line='load '//trim(it1)//', '//trim(it2)
              d_line='drop '//trim(it2)//', '//trim(it1)
          end if
      end if

      ! Convert factor to string, ignore if 1.0 or -1.0
      if (fact.ne.1.0) then
          if (fact.ne.-1.0) then
              write(sfact,*) fact
              sfact_star=' '//trim(sfact)//'*'
          end if
      end if

      ! Print the ITF line
      if (interm_res) then
          write(lulog,*) a_line
      end if

      write(lulog,*) l_line

      if (fact.lt.0.0) then
          itf_line='.'//trim(adjustl(nres))//'['//trim(idx3)//'] -= '//
     &        trim(sfact_star)//
     &        trim(adjustl(nt1))//'['//trim(idx1)//'] '//
     &        trim(adjustl(nt2))//'['//trim(idx2)//']'
      else
          itf_line='.'//trim(adjustl(nres))//'['//trim(idx3)//'] += '//
     &        trim(sfact_star)//
     &        trim(adjustl(nt1))//'['//trim(idx1)//'] '//
     &        trim(adjustl(nt2))//'['//trim(idx2)//']'
      end if

      write(lulog,*) trim(itf_line)

      write(lulog,*) d_line

      interm_res=.false.
      interm_t1=.false.
      interm_t2=.false.

      return
      end

*----------------------------------------------------------------------*
      subroutine assign_index(contr_info,istr1,istr2,istr3)
*----------------------------------------------------------------------*
!     
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      
      type(binary_contr), intent(in) ::
     &     contr_info   ! Inofrmation about binary contraction
      character(len=8), intent(inout)  ::
     &     istr1,       ! Operator 1 index
     &     istr2,       ! Operator 2 index
     &     istr3        ! Result index
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
      character(len=4) ::
     &     c1='    ',    ! Workspace to assign index before array
     &     c2='    ',
     &     c3='    ',
     &     c4='    ',
     &     a1='    ',
     &     a2='    ',
     &     a3='    ',
     &     a4='    '
      character(len=4), dimension(8) ::
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


      istr1='        '
      istr2='        '
      istr3='        '

      ! Set index convention
      idx_type=(/ 0, 0, 0 /)
      do i=1, len(tensor_ham)
          if (contr_info%label_op1.eq.trim(tensor_ham(i))) then
              idx_type(1)=1
          end if
          if (contr_info%label_op2.eq.trim(tensor_ham(i))) then
              idx_type(2)=1
          end if
          if (contr_info%label_res.eq.trim(tensor_ham(i))) then
              idx_type(3)=1
          end if
      end do

      ! Get occuation info
      do i = 1, 1 ! Just the first block
        call count_index(i,
     &     contr_info%occ_op1(1:,1:,i),
     &     contr_info%rst_op1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,t1)
        call count_index(i,
     &     contr_info%occ_cnt(1:,1:,i),
     &     contr_info%rst_cnt(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,c)
        call count_index(i,
     &     contr_info%occ_ex1(1:,1:,i),
     &     contr_info%rst_ex1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e1)
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

      c1='    '
      c2='    '
      c3='    '
      a1='    '
      a2='    '
      a3='    '
      
      
      ! Assign c (contracted by)
      ! These will be assigned as above
      do i=1, c(2,1)
          c1(i:)=par(i)
      end do
      c_array(1)=c1
      do i=1, c(3,1)
          c2(i:)=val(i)
      end do
      c_array(2)=c2
      do i=1, c(1,1)
          c3(i:)=hol(i)
      end do
      c_array(3)=c3

      ! Need to to be shifted so to match assignment of t1 above
      do i=1, c(2,2)
          a1(i:)=par(i+t1(2,1))
      end do
      c_array(5)=a1
      do i=1, c(3,2)
          a2(i:)=val(i+t1(3,1))
      end do
      c_array(6)=a2
      do i=1, c(1,2)
          a3(i:)=hol(i+t1(1,1))
      end do
      c_array(7)=a3

      c1='    '
      c2='    '
      c3='    '
      a1='    '
      a2='    '
      a3='    '


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
          a1(i:)=par(i+t1(2,1))
      end do
      e1_array(5)=a1
      do i=1, e1(3,2)
          a2(i:)=val(i+t1(3,1))
      end do
      e1_array(6)=a2
      do i=1, e1(1,2)
          a3(i:)=hol(i+t1(1,1))
      end do
      e1_array(7)=a3

      c1='    '
      c2='    '
      c3='    '
      a1='    '
      a2='    '
      a3='    '


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
      
      c1='    '
      c2='    '
      c3='    '
      a1='    '
      a2='    '
      a3='    '


      ! Construct final index strings
      ! Operator 1
      select case(idx_type(1))
      case(ham)
          ! Hamiltonian/integral convention
          istr1=trim(adjustl(t1_array(1)))//trim(adjustl(t1_array(5)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(7)))//
     &          trim(adjustl(t1_array(2)))//trim(adjustl(t1_array(6)))
      case default
          ! [apij] (aacc), ie. T[abij]
          istr1=trim(adjustl(t1_array(1)))//trim(adjustl(t1_array(2)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(5)))//
     &          trim(adjustl(t1_array(6)))//trim(adjustl(t1_array(7)))
      end select


      ! Operator 2
      ! c_array annhilations correspond to t2 creations and vice versa
      select case(idx_type(2))
      case(ham)
          istr2=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(5)))//
     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(3)))//
     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(7)))//
     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(6)))
      case default
          istr2=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(5)))//
     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(6)))//
     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(7)))//
     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(3)))
      end select
      

      ! Result
      select case(idx_type(3))
      case(ham)
          istr3=trim(adjustl(e1_array(1)))//trim(adjustl(e2_array(1)))//
     &          trim(adjustl(e1_array(5)))//trim(adjustl(e2_array(5)))//
     &          trim(adjustl(e1_array(3)))//trim(adjustl(e2_array(3)))//
     &          trim(adjustl(e1_array(7)))//trim(adjustl(e2_array(7)))//
     &          trim(adjustl(e1_array(2)))//trim(adjustl(e2_array(2)))//
     &          trim(adjustl(e1_array(6)))//trim(adjustl(e2_array(6)))
      case default
          istr3=trim(adjustl(e1_array(1)))//trim(adjustl(e2_array(1)))//
     &          trim(adjustl(e1_array(2)))//trim(adjustl(e2_array(2)))//
     &          trim(adjustl(e1_array(3)))//trim(adjustl(e2_array(3)))//
     &          trim(adjustl(e1_array(5)))//trim(adjustl(e2_array(5)))//
     &          trim(adjustl(e1_array(6)))//trim(adjustl(e2_array(6)))//
     &          trim(adjustl(e1_array(7)))//trim(adjustl(e2_array(7)))
      end select
      
      return
      end

*----------------------------------------------------------------------*
      subroutine itf_tensor_init(contr_info,itf1,itf2,itf3)
*----------------------------------------------------------------------*
!     
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_tensor.h'
      
      type(binary_contr), intent(in) ::
     &     contr_info   ! Inofrmation about binary contraction
      type(itf_tensor), intent(inout)  ::
     &     itf1,       ! Operator 1 index
     &     itf2,       ! Operator 2 index
     &     itf3        ! Result index
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
      character(len=4) ::
     &     c1='    ',    ! Workspace to assign index before array
     &     c2='    ',
     &     c3='    ',
     &     c4='    ',
     &     a1='    ',
     &     a2='    ',
     &     a3='    ',
     &     a4='    '
      character(len=4), dimension(8) ::
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


      itf1%idx='        '
      itf1%idx='        '
      itf1%idx='        '

      ! Set index convention
      idx_type=(/ 0, 0, 0 /)
      do i=1, len(tensor_ham)
          if (contr_info%label_op1.eq.trim(tensor_ham(i))) then
              idx_type(1)=1
          end if
          if (contr_info%label_op2.eq.trim(tensor_ham(i))) then
              idx_type(2)=1
          end if
          if (contr_info%label_res.eq.trim(tensor_ham(i))) then
              idx_type(3)=1
          end if
      end do

      ! Get occuation info
      do i = 1, 1 ! Just the first block
        call count_index(i,
     &     contr_info%occ_op1(1:,1:,i),
     &     contr_info%rst_op1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,t1)
        call count_index(i,
     &     contr_info%occ_cnt(1:,1:,i),
     &     contr_info%rst_cnt(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,c)
        call count_index(i,
     &     contr_info%occ_ex1(1:,1:,i),
     &     contr_info%rst_ex1(1:,1:,1:,1:,1:,i),
     &     contr_info%ngas,contr_info%nspin,e1)
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

      itf1%c_ops=trim(c1)//trim(c2)//trim(c3)
      itf1%a_ops=trim(a1)//trim(a2)//trim(a3)

      c1='    '
      c2='    '
      c3='    '
      a1='    '
      a2='    '
      a3='    '
      
      
      ! Assign c (contracted by)
      ! These will be assigned as above
      do i=1, c(2,1)
          c1(i:)=par(i)
      end do
      c_array(1)=c1
      do i=1, c(3,1)
          c2(i:)=val(i)
      end do
      c_array(2)=c2
      do i=1, c(1,1)
          c3(i:)=hol(i)
      end do
      c_array(3)=c3

      ! Need to to be shifted so to match assignment of t1 above
      do i=1, c(2,2)
          a1(i:)=par(i+t1(2,1))
      end do
      c_array(5)=a1
      do i=1, c(3,2)
          a2(i:)=val(i+t1(3,1))
      end do
      c_array(6)=a2
      do i=1, c(1,2)
          a3(i:)=hol(i+t1(1,1))
      end do
      c_array(7)=a3

      c1='    '
      c2='    '
      c3='    '
      a1='    '
      a2='    '
      a3='    '


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
          a1(i:)=par(i+t1(2,1))
      end do
      e1_array(5)=a1
      do i=1, e1(3,2)
          a2(i:)=val(i+t1(3,1))
      end do
      e1_array(6)=a2
      do i=1, e1(1,2)
          a3(i:)=hol(i+t1(1,1))
      end do
      e1_array(7)=a3

      c1='    '
      c2='    '
      c3='    '
      a1='    '
      a2='    '
      a3='    '


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
      
      c1='    '
      c2='    '
      c3='    '
      a1='    '
      a2='    '
      a3='    '


      ! Construct final index strings
      ! Operator 1
      select case(idx_type(1))
      case(ham)
       ! Hamiltonian/integral convention
       itf1%idx=trim(adjustl(t1_array(1)))//trim(adjustl(t1_array(5)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(7)))//
     &          trim(adjustl(t1_array(2)))//trim(adjustl(t1_array(6)))
      case default
       ! [apij] (aacc), ie. T[abij]
       itf1%idx=trim(adjustl(t1_array(1)))//trim(adjustl(t1_array(2)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(5)))//
     &          trim(adjustl(t1_array(6)))//trim(adjustl(t1_array(7)))
      end select


      ! Operator 2
      ! c_array annhilations correspond to t2 creations and vice versa
      select case(idx_type(2))
      case(ham)
       itf2%idx=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(5)))//
     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(3)))//
     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(7)))//
     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(6)))
      case default
       itf2%idx=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(5)))//
     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(6)))//
     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(7)))//
     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(3)))
      end select
      
      itf2%c_ops=trim(adjustl(e2_array(1)))//trim(adjustl(c_array(5)))//
     &           trim(adjustl(e2_array(2)))//trim(adjustl(c_array(6)))//
     &           trim(adjustl(e2_array(3)))//trim(adjustl(c_array(7)))

      itf2%a_ops=trim(adjustl(e2_array(5)))//trim(adjustl(c_array(1)))//
     &           trim(adjustl(e2_array(6)))//trim(adjustl(c_array(2)))//
     &           trim(adjustl(e2_array(7)))//trim(adjustl(c_array(3)))

      ! Result
      select case(idx_type(3))
      case(ham)
       itf3%idx=trim(adjustl(e1_array(1)))//trim(adjustl(e2_array(1)))//
     &          trim(adjustl(e1_array(5)))//trim(adjustl(e2_array(5)))//
     &          trim(adjustl(e1_array(3)))//trim(adjustl(e2_array(3)))//
     &          trim(adjustl(e1_array(7)))//trim(adjustl(e2_array(7)))//
     &          trim(adjustl(e1_array(2)))//trim(adjustl(e2_array(2)))//
     &          trim(adjustl(e1_array(6)))//trim(adjustl(e2_array(6)))
      case default
       itf3%idx=trim(adjustl(e1_array(1)))//trim(adjustl(e2_array(1)))//
     &          trim(adjustl(e1_array(2)))//trim(adjustl(e2_array(2)))//
     &          trim(adjustl(e1_array(3)))//trim(adjustl(e2_array(3)))//
     &          trim(adjustl(e1_array(5)))//trim(adjustl(e2_array(5)))//
     &          trim(adjustl(e1_array(6)))//trim(adjustl(e2_array(6)))//
     &          trim(adjustl(e1_array(7)))//trim(adjustl(e2_array(7)))
      end select

      itf3%c_ops=trim(adjustl(e1_array(1)))//
     &           trim(adjustl(e2_array(1)))//
     &          trim(adjustl(e1_array(2)))//trim(adjustl(e2_array(2)))//
     &          trim(adjustl(e1_array(3)))//trim(adjustl(e2_array(3)))

      itf3%a_ops=trim(adjustl(e1_array(5)))//
     &          trim(adjustl(e2_array(5)))//
     &          trim(adjustl(e1_array(6)))//trim(adjustl(e2_array(6)))//
     &          trim(adjustl(e1_array(7)))//trim(adjustl(e2_array(7)))

      itf1%name=contr_info%label_op1
      itf2%name=contr_info%label_op2
      itf3%name=contr_info%label_res
      
      itf1%rank=len(trim(itf1%idx))
      itf2%rank=len(trim(itf2%idx))
      itf3%rank=len(trim(itf3%idx))

      itf1%fact=contr_info%fact
      itf2%fact=contr_info%fact
      itf3%fact=1

      return
      end
