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
      subroutine spin_sum_index(istr1, istr2, istr3, lulog)
*----------------------------------------------------------------------*
!     Spin sum index and produce resulting binary contractions
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'

      character(len=8), intent(in) ::
     &    istr1, istr2, istr3          ! Op1, op2, res
      integer, intent(in) ::
     &    lulog
      integer ::
     &    r1, r2, r3,               ! rank
     &    e1, e2, e3, e4,              ! external indicies
     &    c1, c2,                   ! contraction indicies
     &    i,j,k,l,m,n,
     &    extern

      r1=len(trim(istr1))
      r2=len(trim(istr2))
      r3=4
      extern=2

      if (r3.eq.4) then
         e1=1
         do i=0, 1
            e2=1
            e1=e1+i
            do j=0, 1
               e3=1
               e2=e2+j
               do k=0, 1
                  e4=1
                  e3=e3+k
                  do l=0, 1
                     e4=e4+l
                     !write(lulog,*) "hello", e1, e2, e3, e4
                     if (e1==1 .and. e2==1 .and. e3==2 .and. e4==2) then
                        c1=1
                        do m=0, 1
                           c2=1
                           c1=c1+m
                           do n=0, 1
                              c2=c2+n
                              select case(extern)
                                 case(0)
                                    if (mod(e1+e2+e3+e4,2)==0) then
                                       write(lulog,*) "A"
                                       write(lulog,*) "B ", e1, e2, e3,
     &                                                 e4
                                    end if
                                 case(2)
                                    if (mod(e1+e2+c1+c2,2)==0) then
                                       write(lulog,*) "A ", e1, e2,
     &                                                c1,
     &                                                c2
                                       write(lulog,*) "B ", e3, e4, c1,
     &                                                c2
                                    end if
                               end select
                           end do
                        end do
                     end if
                  end do
               end do
            end do
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
      character(len=8), intent(in) ::
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
!      if (interm_res) then
!          write(lulog,*) a_line
!      end if
!
!      write(lulog,*) l_line

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

!      write(lulog,*) d_line

      interm_res=.false.
      interm_t1=.false.
      interm_t2=.false.

      return
      end

*----------------------------------------------------------------------*
      subroutine assign_index(contr_info,istr1,istr2,istr3,lulog)
*----------------------------------------------------------------------*
!     
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      
      type(binary_contr), intent(in) ::
     &     contr_info   ! Inofrmation about binary contraction
      integer, intent(in) ::
     &     lulog        ! Debug delete
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


      istr1='        '
      istr2='        '
      istr3='        '

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


      call permute_tensors(e1,e2,c,istr1,istr2,istr3,lulog)

      return
      end

*----------------------------------------------------------------------*
      subroutine permute_tensors(e1,e2,c,istr1,istr2,istr3,lulog)
*----------------------------------------------------------------------*
!     
*----------------------------------------------------------------------*

      implicit none

      character(len=8), intent(in)  ::
     &     istr1,       ! Operator 1 index
     &     istr2,       ! Operator 2 index
     &     istr3        ! Result index
      integer, intent(in) ::
     &     e1(4,2),      ! Occupations of external index 1
     &     e2(4,2),      ! Occupations of external index 2
     &     c(4,2)
      integer, intent(in) ::
     &     lulog
      character(len=1) ::
     &     t1a(2),   ! First group of first tensor
     &     t1b(2),   ! Second group of first tensor
     &     t2a(2),
     &     t2b(2),
     &     r1a(2),
     &     r1b(2)
      integer ::
     &     i,j,k,l,m,n,o,p,q,
     &     sum_c1,sum_c2,sum_a1,sum_a2,
     &     s1a(2),
     &     s1b(2),
     &     s2a(2),
     &     s2b(2),
     &     second_idx,third_idx,
     &     zero_a,zero_b

      ! Check if not antisym over different verticies

      ! For rank 4. Rank 2 and 0 don't need antisymetrising
      if (len(trim(istr3))==0 .or.
     &    len(trim(istr3))==6) then
         return
      end if

      ! Don't care about tensor products now
      if (len(trim(istr3))==4 .and. len(trim(istr1))==2 .and.
     &    len(trim(istr2))==2) then
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

      if (len(trim(istr1))==4 .and. len(trim(istr2))==4) then
         t1a(1)=istr1(1:1)
         t1a(2)=istr1(2:2)
         t1b(1)=istr1(3:3)
         t1b(2)=istr1(4:4)

         t2a(1)=istr2(1:1)
         t2a(2)=istr2(2:2)
         t2b(1)=istr2(3:3)
         t2b(2)=istr2(4:4)
      else if (len(trim(istr1))==4 .and. len(trim(istr2))==2) then
         ! Takes care of rank 4 and rank 2 result
         t1a(1)=istr1(1:1)
         t1a(2)=istr1(2:2)
         t1b(1)=istr1(3:3)
         t1b(2)=istr1(4:4)

         t2a(1)=istr2(1:1)
         t2b(1)=istr2(2:2)
      else if (len(trim(istr1))==2 .and. len(trim(istr2))==4) then
         ! Swap t1 and t2 if t1=rank 2. Only going to sum over
         ! contraction indices of t1.
         ! Takes care of rank 4 and rank 2 result
         t2a(1)=istr1(1:1)
         t2b(1)=istr1(2:2)

         t1a(1)=istr2(1:1)
         t1a(2)=istr2(2:2)
         t1b(1)=istr2(3:3)
         t1b(2)=istr2(4:4)
      else if (len(trim(istr1))==4 .and. len(trim(istr2))==0) then
         t1a(1)=istr1(1:1)
         t1a(2)=istr1(2:2)
         t1b(1)=istr1(3:3)
         t1b(2)=istr1(4:4)
      else if (len(trim(istr1))==0 .and. len(trim(istr2))==4) then
         ! Swap t1 and t2 if t1=rank 0
         t1a(1)=istr2(1:1)
         t1a(2)=istr2(2:2)
         t1b(1)=istr2(3:3)
         t1b(2)=istr2(4:4)
      else if (len(trim(istr1))==2 .and. len(trim(istr2))==2 
     &         .and. len(trim(istr3))==2) then
         ! Don't need to spin sum
         return
      else if (len(trim(istr1))==0 .and. len(trim(istr2))==2
     &         .and. len(trim(istr3))==2) then
         ! Don't need to spin sum
         return
      else if (len(trim(istr1))==2 .and. len(trim(istr2))==0 
     &         .and. len(trim(istr3))==2) then
         ! Don't need to spin sum
         return
      end if

      if (len(trim(istr3))==4) then
         r1a(1)=istr3(1:1)
         r1a(2)=istr3(2:2)
         r1b(1)=istr3(3:3)
         r1b(2)=istr3(4:4)
      else if (len(trim(istr3))==2) then
         r1a(1)=istr3(1:1)
         r1b(1)=istr3(2:2)
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
            write(lulog,*) "permute creations! 0.5*(1-P)"
         end if
      end if

      if (sum_a1/=2 .and. sum_a2/=2) then
         if (sum_a1+sum_a2==2) then
            write(lulog,*) "permute annhilations! 0.5*(1-P)"
         end if
      end if
      

      ! Sum over contraction creation

      ! Working in index groups, set abab (1212) index to individual tensor
      ! index groups. Ordering of spins doesn't matter, only overall Sz.
      do i=1, 2
         if (r1a(1)==t1a(i) .and. r1a(1)/='') then
            s1a(i)=1
         end if
         if (r1a(1)==t1b(i) .and. r1a(1)/='') then
            s1b(i)=1
         end if
    
         if (r1a(2)==t1a(i) .and. r1a(2)/='') then
            s1a(i)=2
         end if
         if (r1a(2)==t1b(i) .and. r1a(2)/='') then
            s1b(i)=2
         end if
    
    
         if (r1b(1)==t1a(i) .and. r1b(1)/='') then
            s1a(i)=1
         end if
         if (r1b(1)==t1b(i) .and. r1b(1)/='') then
            s1b(i)=1
         end if
    
         if (r1b(2)==t1a(i) .and. r1b(2)/='') then
            s1a(i)=2
         end if
         if (r1b(2)==t1b(i) .and. r1b(2)/='') then
            s1b(i)=2
         end if
      end do

      do i=1, 2
         if (r1a(1)==t2a(i) .and. r1a(1)/='') then
            s2a(i)=1
         end if
         if (r1a(1)==t2b(i) .and. r1a(1)/='') then
            s2b(i)=1
         end if
    
         if (r1a(2)==t2a(i) .and. r1a(2)/='') then
            s2a(i)=2
         end if
         if (r1a(2)==t2b(i) .and. r1a(2)/='') then
            s2b(i)=2
         end if
    
    
         if (r1b(1)==t2a(i) .and. r1b(1)/='') then
            s2a(i)=1
         end if
         if (r1b(1)==t2b(i) .and. r1b(1)/='') then
            s2b(i)=1
         end if
    
         if (r1b(2)==t2a(i) .and. r1b(2)/='') then
            s2a(i)=2
         end if
         if (r1b(2)==t2b(i) .and. r1b(2)/='') then
            s2b(i)=2
         end if
      end do

      ! Index assigned from result
      write(lulog,*) "s1b: ", s1b
      write(lulog,*) "s1a: ", s1a
      write(lulog,*)
      write(lulog,*) "s2b: ", s2b
      write(lulog,*) "s2a: ", s2a
      write(lulog,*)
      write(lulog,*)

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
      second_idx=0
      third_idx=0
      ! Check for unassigned indices in first index group
      do i=1, size(s1a)
         if (s1a(i)==0) then
            do j=1, 2
               s1a(i)=j
               ! Find index in second tensor
               do m=1, size(t2a)
                  if (t2a(m)==t1a(i)) then
                     s2a(m)=j
                  else if (t2b(m)==t1a(i)) then
                     s2b(m)=j
                  end if
               end do
               ! Check if second index is 0
               do k=i, size(s1a)
                  if (s1a(k)==0 .or. second_idx>0
     &                .and. k==second_idx) then
                     ! If this is second time, s1a(k)/=0, so need to
                     ! remeber where contraction index was
                     ! Doubt this will work for rank 6 tensors...
                     if (second_idx==0) then
                        second_idx=k
                     end if
                     do l=1, 2
                        s1a(second_idx)=l
                        ! Find index in second tensor
                        do n=1, size(t2a)
                           if (t2a(n)==t1a(second_idx)) then
                              s2a(n)=l
                           else if (t2b(n)==t1a(second_idx)) then
                              s2b(n)=l
                           end if
                        end do
                        ! Check if contraction indicies in second index
                        ! group
                        if (any(s1b==0) .or. third_idx>0) then
                           do o=1, size(s1b)
                              if (s1b(o)==0 .or. third_idx>0
     &                            .and. o==third_idx) then
                                 if (third_idx==0) then
                                    third_idx=o
                                 end if
                                 do p=1, 2
                                    s1b(third_idx)=p
                                    do q=1, size(t2a)
                                       if (t2a(q)==t1b(third_idx)) then
                                          s2a(q)=p
                                       else if (t2b(q)==t1b(third_idx))
     &                                 then
                                          s2b(q)=p
                                       end if
                                    end do
                                    write(lulog,*) "s1b: ", s1b
                                    write(lulog,*) "s1a: ", s1a
                                    write(lulog,*)
                                    write(lulog,*) "s2b: ", s2b
                                    write(lulog,*) "s2a: ", s2a
                                    write(lulog,*)
                                    write(lulog,*)
                                 end do
                              end if
                           end do
                        else   
                           write(lulog,*) "s1b: ", s1b
                           write(lulog,*) "s1a: ", s1a
                           write(lulog,*)
                           write(lulog,*) "s2b: ", s2b
                           write(lulog,*) "s2a: ", s2a
                           write(lulog,*)
                           write(lulog,*)
                        end if
                     end do
                  else if (zero_a==1 .and. zero_b==1) then
                     ! Case where one contraction over s1a and s1b
                     do o=1, size(s1b)
                        if (s1b(o)==0 .or. third_idx>0
     &                      .and. o==third_idx) then
                           if (third_idx==0) then
                              third_idx=o
                           end if
                           do p=1, 2
                              s1b(third_idx)=p
                              do q=1, size(t2a)
                                 if (t2a(q)==t1b(third_idx)) then
                                    s2a(q)=p
                                 else if (t2b(q)==t1b(third_idx))
     &                           then
                                    s2b(q)=p
                                 end if
                              end do
                              write(lulog,*) "s1b: ", s1b
                              write(lulog,*) "s1a: ", s1a
                              write(lulog,*)
                              write(lulog,*) "s2b: ", s2b
                              write(lulog,*) "s2a: ", s2a
                              write(lulog,*)
                              write(lulog,*)
                           end do
                        end if
                     end do
                  else if (any(s1a/=0) .and. k==2) then
                     ! No more contraction indicies, print out result
                     ! Only prints if k==2, ie. at end of k loop, won't
                     ! work for rank 6 tensors
                     write(lulog,*) "s1b: ", s1b
                     write(lulog,*) "s1a: ", s1a
                     write(lulog,*)
                     write(lulog,*) "s2b: ", s2b
                     write(lulog,*) "s2a: ", s2a
                     write(lulog,*)
                     write(lulog,*)
                  end if
               end do ! Check if second index is 0
            end do ! Loop over a/b for first index
         end if ! Check for first 0 index
      end do

      else
      second_idx=0
      third_idx=0
      ! Check for unassigned indices in second index group
      ! This assumes both contraction indices are in the second index
      ! group, above loops check the remaing two cases
      do i=1, size(s1b)
         if (s1b(i)==0) then
            do j=1, 2
               s1b(i)=j
               ! Find index in second tensor
               do m=1, size(t2a)
                  if (t2a(m)==t1b(i)) then
                     s2a(m)=j
                  else if (t2b(m)==t1b(i)) then
                     s2b(m)=j
                  end if
               end do
               ! Check if second index is 0
               do k=i, size(s1b)
                  if (s1b(k)==0 .or. second_idx>0
     &                .and. k==second_idx) then
                     ! If this is second time, s1b(k)/=0, so need to
                     ! remeber where contraction index was
                     ! Doubt this will work for rank 6 tensors...
                     if (second_idx==0) then
                        second_idx=k
                     end if
                     do l=1, 2
                        s1b(second_idx)=l
                        ! Find index in second tensor
                        do n=1, size(t2a)
                           if (t2a(n)==t1b(second_idx)) then
                              s2a(n)=l
                           else if (t2b(n)==t1b(second_idx)) then
                              s2b(n)=l
                           end if
                        end do
                        ! Check if contraction indicies in second index
                        ! group
                        if (any(s1a==0) .or. third_idx>0) then
                           do o=1, size(s1a)
                              if (s1a(o)==0 .or. third_idx>0
     &                            .and. o==third_idx) then
                                 if (third_idx==0) then
                                    third_idx=o
                                 end if
                                 do p=1, 2
                                    s1a(third_idx)=p
                                    do q=1, size(t2a)
                                       if (t2a(q)==t1a(third_idx)) then
                                          s2a(q)=p
                                       else if (t2b(q)==t1a(third_idx))
     &                                 then
                                          s2b(q)=p
                                       end if
                                    end do
                                    write(lulog,*) "s1b: ", s1b
                                    write(lulog,*) "s1a: ", s1a
                                    write(lulog,*)
                                    write(lulog,*) "s2b: ", s2b
                                    write(lulog,*) "s2a: ", s2a
                                    write(lulog,*)
                                    write(lulog,*)
                                 end do
                              end if
                           end do
                        else
                           write(lulog,*) "s1b: ", s1b
                           write(lulog,*) "s1a: ", s1a
                           write(lulog,*)
                           write(lulog,*) "s2b: ", s2b
                           write(lulog,*) "s2a: ", s2a
                           write(lulog,*)
                           write(lulog,*)
                        end if
                     end do
!                  else if (any(s1a==0) .and. .not. (any(s1b==0))
!     &                     .or. third_idx>0 .and.
!     &                     .not. (any(s1b==0))) then
                   else if (zero_a==1 .and. zero_b==1) then
                     ! Case where one contraction over s1a and s1b
                     do o=1, size(s1a)
                        if (s1a(o)==0 .or. third_idx>0
     &                      .and. o==third_idx) then
                           if (third_idx==0) then
                              third_idx=o
                           end if
                           do p=1, 2
                              s1a(third_idx)=p
                              do q=1, size(t2a)
                                 if (t2a(q)==t1a(third_idx)) then
                                    s2a(q)=p
                                 else if (t2b(q)==t1a(third_idx))
     &                           then
                                    s2b(q)=p
                                 end if
                              end do
                              write(lulog,*) "hello"
                              write(lulog,*) "s1b: ", s1b
                              write(lulog,*) "s1a: ", s1a
                              write(lulog,*)
                              write(lulog,*) "s2b: ", s2b
                              write(lulog,*) "s2a: ", s2a
                              write(lulog,*)
                              write(lulog,*)
                           end do
                        end if
                     end do
                  else if (any(s1b/=0) .and. k==2) then
                     ! No more contraction indicies, print out result
                     ! Only prints if k==2, ie. at end of k loop, won't
                     ! work for rank 6 tensors
                     write(lulog,*) "s1b: ", s1b
                     write(lulog,*) "s1a: ", s1a
                     write(lulog,*)
                     write(lulog,*) "s2b: ", s2b
                     write(lulog,*) "s2a: ", s2a
                     write(lulog,*)
                     write(lulog,*)
                  end if
               end do ! Check if second index is 0
            end do ! Loop over a/b for first index
         end if ! Check for first 0 index
      end do

      end if

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
     &     t1_array,    ! Operator 1 array
     &     c_array,     ! Contraction index array
     &     e1_array,    ! External index of operator 1 array
     &     e2_array     ! External index of operator 2 array
      integer, parameter ::
     &     t_amp = 0,       ! [apij] (aacc)
     &     ham   = 1        ! [abip]
      character(len=20), dimension(4) ::
     &     tensor_ham=(/ 'H', 'INT_D', 'INT_HH', 'INT_PP' /)  ! Tensor names to use ham index convention


      itf1%idx='        '
      itf1%idx='        '
      itf1%idx='        '

      ! Tensor number
      itf1%t_numb=1
      itf2%t_numb=2
      itf3%t_numb=3

      ! Set index convention
      ! For now, don't bother
!      idx_type=(/ 0, 0, 0 /)
!      do i=1, len(tensor_ham)
!          if (contr_info%label_op1.eq.trim(tensor_ham(i))) then
!              idx_type(1)=1
!          end if
!          if (contr_info%label_op2.eq.trim(tensor_ham(i))) then
!              idx_type(2)=1
!          end if
!          if (contr_info%label_res.eq.trim(tensor_ham(i))) then
!              idx_type(3)=1
!          end if
!      end do
      itf1%idx_convention=0
      itf2%idx_convention=0
      itf3%idx_convention=0

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
      select case(itf1%idx_convention)
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

       itf1%ncops=len(trim(adjustl(t1_array(1)))//
     &          trim(adjustl(t1_array(2)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(5)))//
     &          trim(adjustl(t1_array(6)))//trim(adjustl(t1_array(7))))
       itf1%naops=len(trim(adjustl(t1_array(1)))//
     &          trim(adjustl(t1_array(2)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(5)))//
     &          trim(adjustl(t1_array(6)))//trim(adjustl(t1_array(7))))


      ! Operator 2
      ! c_array annhilations correspond to t2 creations and vice versa
      select case(itf2%idx_convention)
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
      itf2%ncops=len(trim(adjustl(e2_array(1)))//
     &           trim(adjustl(c_array(5)))//
     &           trim(adjustl(e2_array(2)))//trim(adjustl(c_array(6)))//
     &           trim(adjustl(e2_array(3)))//trim(adjustl(c_array(7))))

      itf2%a_ops=trim(adjustl(e2_array(5)))//trim(adjustl(c_array(1)))//
     &           trim(adjustl(e2_array(6)))//trim(adjustl(c_array(2)))//
     &           trim(adjustl(e2_array(7)))//trim(adjustl(c_array(3)))
      itf2%naops=len(trim(adjustl(e2_array(5)))//
     &           trim(adjustl(c_array(1)))//
     &           trim(adjustl(e2_array(6)))//trim(adjustl(c_array(2)))//
     &           trim(adjustl(e2_array(7)))//trim(adjustl(c_array(3))))

      ! Result
      select case(itf3%idx_convention)
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
     &          trim(adjustl(e2_array(1)))//
     &          trim(adjustl(e1_array(2)))//trim(adjustl(e2_array(2)))//
     &          trim(adjustl(e1_array(3)))//trim(adjustl(e2_array(3)))
      itf3%ncops=len(trim(adjustl(e1_array(1)))//
     &          trim(adjustl(e2_array(1)))//
     &          trim(adjustl(e1_array(2)))//trim(adjustl(e2_array(2)))//
     &          trim(adjustl(e1_array(3)))//trim(adjustl(e2_array(3))))

      itf3%a_ops=trim(adjustl(e1_array(5)))//
     &          trim(adjustl(e2_array(5)))//
     &          trim(adjustl(e1_array(6)))//trim(adjustl(e2_array(6)))//
     &          trim(adjustl(e1_array(7)))//trim(adjustl(e2_array(7)))
      itf3%naops=len(trim(adjustl(e1_array(5)))//
     &          trim(adjustl(e2_array(5)))//
     &          trim(adjustl(e1_array(6)))//trim(adjustl(e2_array(6)))//
     &          trim(adjustl(e1_array(7)))//trim(adjustl(e2_array(7))))

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

*----------------------------------------------------------------------*
      subroutine tensor_idx(contr_info,itf1,itf2,itf3)
*----------------------------------------------------------------------*
!     
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_tensor.h'
      
      type(binary_contr), intent(in) ::
     &     contr_info   ! Inofrmation about binary contraction
!      type(itf_tensor), intent(inout)  ::
!     &     itf1,       ! Operator 1 index
!     &     itf2,       ! Operator 2 index
!     &     itf3        ! Result index
      character(len=index_len) ::
     &     itf1, itf2, itf3
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
     &     t1_array,    ! Operator 1 array
     &     c_array,     ! Contraction index array
     &     e1_array,    ! External index of operator 1 array
     &     e2_array     ! External index of operator 2 array
      integer, parameter ::
     &     t_amp = 0,       ! [apij] (aacc)
     &     ham   = 1        ! [abip]
      character(len=20), dimension(4) ::
     &     tensor_ham=(/ 'H', 'INT_D', 'INT_HH', 'INT_PP' /)  ! Tensor names to use ham index convention


      ! Get occuation info
      do i = 1, 1   ! Just one block for now
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
       itf1=trim(adjustl(t1_array(1)))//
     &          trim(adjustl(t1_array(2)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(t1_array(5)))//
     &          trim(adjustl(t1_array(6)))//trim(adjustl(t1_array(7)))

      ! Operator 2
      ! c_array annhilations correspond to t2 creations and vice versa
       itf2=trim(adjustl(e2_array(1)))//
     &          trim(adjustl(c_array(5)))//
     &          trim(adjustl(e2_array(2)))//trim(adjustl(c_array(6)))//
     &          trim(adjustl(e2_array(3)))//trim(adjustl(c_array(7)))//
     &          trim(adjustl(e2_array(5)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(e2_array(6)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(e2_array(7)))//trim(adjustl(c_array(3)))
      
      ! Result
       itf3=trim(adjustl(e1_array(1)))//
     &          trim(adjustl(e2_array(1)))//
     &          trim(adjustl(e1_array(2)))//trim(adjustl(e2_array(2)))//
     &          trim(adjustl(e1_array(3)))//trim(adjustl(e2_array(3)))//
     &          trim(adjustl(e1_array(5)))//trim(adjustl(e2_array(5)))//
     &          trim(adjustl(e1_array(6)))//trim(adjustl(e2_array(6)))//
     &          trim(adjustl(e1_array(7)))//trim(adjustl(e2_array(7)))

      return
      end
