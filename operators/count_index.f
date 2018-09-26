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
      subroutine bcontr_to_itf(contr_info,itflog)
*----------------------------------------------------------------------*
!     Take GeCco binary contraction and produce ITF algo code.
!     Includes antisymmetry of residual equations and spin summation.
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      type(binary_contr), intent(in) ::
     &     contr_info      ! Inofrmation about binary contraction
      integer, intent(in) ::
     &     itflog          ! Output file

      type(itf_contr) ::
     &     itf_item        ! ITF contraction object; holds all info about the ITF algo line
      integer ::
     &    perm_array(4),   ! Info of antisymmetrisation permutation factors
     &    i                ! Loop index
      logical ::
     &    antisymm=.false. ! Include extra permuation factors for (1-P_ab)(1-P_ij)

      ! Initalise permutation factors to 0 == no permutation
      perm_array=0
      ! Determine if result needs antisymmetrising
      !call permute_tensors(contr_info,perm_array,itflog)
      !do i=1, size(perm_array)
      !   if (perm_array(i)==4) then
      !      ! Add extra permtation factos
      !      antisymm=.true.
      !   end if
      !end do
      ! Skip antisymmetrising each contraction and antisymm the complete
      ! tensor at the end

      if (perm_array(1)==0) then
         ! No permutations
         call itf_contr_init(contr_info,itf_item,perm_array(1),
     &                       antisymm,itflog)
         call assign_spin(itf_item)
      else
         do i=1, size(perm_array)
            ! Loop over permuation cases and send seperatley to
            ! assign_spin
            call itf_contr_init(contr_info,itf_item,perm_array(i),
     &                          antisymm,itflog)
            call assign_spin(itf_item)
            if (perm_array(i+1)==0) exit
         end do
      end if
      
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
      subroutine t_index(idx,tidx)
*----------------------------------------------------------------------*
!     Transpose ITF index string, abij => abji
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=index_len), intent(in) ::
     &     idx       ! ITF index string
      character(len=index_len), intent(inout) ::
     &     tidx      ! Transpose of ITF index string

      character(len=1) ::
     &     tmp

      tidx=idx
      tmp=idx(1:1)
      tidx(1:1)=idx(2:2)
      tidx(2:2)=tmp
      
      return
      end

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
     &     ires, it1, it2          ! Name of tensors + index
      character(len=index_len) ::
     &     tidx1, tidx2
      character(len=5) ::
     &     s_int                 ! Intermdiate tensor number
      character(len=70) ::
     &     itf_line,             ! Line of ITF code
     &     d_line,             ! Line of ITF code
     &     l_line,             ! Line of ITF code
     &     a_line,             ! Line of ITF code
     &     s_line,             ! Line of ITF code
     &     st1, st2           ! Name of spin summed tensors + index
      character(len=25) ::
     &     sfact,             ! String representation of factor
     &     sfact_star,        ! String representation of factor formatted for output
     &     s_idx              ! Spin index used to group lines
      character(len=maxlen_bc_label) ::
     &     rename_tensor
      integer ::
     &     i

      sfact='                         '
      sfact_star='                         '
      ! Remove these...
      nres=item%label_res
      nt1=item%label_t1
      nt2=item%label_t2

      ! Change names of specific tensors
      nres=rename_tensor(nres)
      
      ! Spin summ
      if (s1) then
         ! Pure spin
         tidx1=''
         call t_index(item%idx1,tidx1)
         st1='('//trim(adjustl(nt1))//'['//trim(item%idx1)//']'//
     &       ' - '//trim(adjustl(nt1))//'['//trim(tidx1)//']'//')'
      else
         st1=trim(adjustl(nt1))//'['//trim(item%idx1)//']'
      end if

      if (s2) then
         ! Pure spin
         tidx2=''
         call t_index(item%idx2,tidx2)
         st2='('//trim(adjustl(nt2))//'['//trim(item%idx2)//']'//
     &       ' - '//trim(adjustl(nt2))//'['//trim(tidx2)//']'//')'
      else
         st2=trim(adjustl(nt2))//'['//trim(item%idx2)//']'
      end if

      ! Convert factor to string, ignore if 1.0 or -1.0
      if (item%fact.ne.1.0) then
          if (item%fact.ne.-1.0) then
              write(sfact,*) item%fact
              do i=1, len(sfact)
                 ! Start after decimal point
                 if (sfact(i:i)=='0' .or. sfact(i:i)=='-') then
                    sfact(i:i)=' '
                 end if
              end do
              sfact_star=trim(adjustl(sfact))//'*'
          end if
      end if

      ! Write spin index to string
      ! Probably not needed...
      !write(s_idx,*) item%spin_idx
      s_idx=''

      if (item%fact.lt.0.0) then
          itf_line=trim(adjustl(s_idx))//'.'//trim(adjustl(nres))//
     &        '['//trim(item%idx3)//
     &        '] -= '//trim(sfact_star)//
     &        trim(adjustl(st1))//' '//
     &        trim(adjustl(st2))
      else
          itf_line=trim(adjustl(s_idx))//'.'//trim(adjustl(nres))//
     7        '['//trim(item%idx3)//
     &        '] += '//trim(sfact_star)//
     &        trim(adjustl(st1))//' '//
     &        trim(adjustl(st2))
      end if

      write(item%logfile,*) trim(itf_line)
      !write(item%logfile,*) "---------------------------------------"

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
      real(8) ::
     &     tmp_fact

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

      tmp_fact=item%fact

      ! Permute indicies to get antisymm tensors
      if (item%permute==0) then
         ! No permutations
         do i=1, size(t1_array)
            t1_array(i)=e1_array(i)
            t2_array(i)=e2_array(i)
         end do
      else if (item%permute==1) then
         ! No permutations, but get a factor
         do i=1, size(t1_array)
            t1_array(i)=e1_array(i)
            t2_array(i)=e2_array(i)
         end do
         if (item%perm4) then
            item%fact=tmp_fact*-0.25
         else
            item%fact=tmp_fact*-0.5
         end if
      else if (item%permute==2) then
         ! Permute creations
         do i=1, size(t1_array)/2
            t1_array(i)=e2_array(i)
            t1_array(i+4)=e1_array(i+4)
            t2_array(i)=e1_array(i)
            t2_array(i+4)=e2_array(i+4)
         end do
         if (item%perm4) then
            item%fact=tmp_fact*-0.25
         else
            item%fact=tmp_fact*-0.5
         end if
      else if (item%permute==3) then
         ! Permute annhilations
         do i=1, size(t1_array)/2
            t1_array(i)=e1_array(i)
            t1_array(i+4)=e2_array(i+4)
            t2_array(i)=e2_array(i)
            t2_array(i+4)=e1_array(i+4)
         end do
         if (item%perm4) then
            item%fact=tmp_fact*-0.25
         else
            item%fact=tmp_fact*-0.5
         end if
      else if (item%permute==4) then
         ! Permute creations and annhilations
         do i=1, size(t1_array)
            t1_array(i)=e2_array(i)
            t2_array(i)=e1_array(i)
         end do
         item%fact=tmp_fact*0.25
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
      item%idx1=trim(adjustl(t1_array(1)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(t1_array(2)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(c_array(3)))//
     &          trim(adjustl(t1_array(5)))//trim(adjustl(c_array(5)))//
     &          trim(adjustl(t1_array(6)))//trim(adjustl(c_array(6)))//
     &          trim(adjustl(t1_array(7)))//trim(adjustl(c_array(7)))
      case default
      ! [apij] (aacc), ie. T[abij]
      item%idx1=trim(adjustl(t1_array(1)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(t1_array(2)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(t1_array(3)))//trim(adjustl(c_array(3)))//
     &          trim(adjustl(t1_array(5)))//trim(adjustl(c_array(5)))//
     &          trim(adjustl(t1_array(6)))//trim(adjustl(c_array(6)))//
     &          trim(adjustl(t1_array(7)))//trim(adjustl(c_array(7)))
      end select

      ! Operator 2
      ! c_array annhilations correspond to t2 creations and vice versa
      select case(idx_type(2))
      case(ham)
      item%idx2=trim(adjustl(t2_array(1)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(t2_array(5)))//trim(adjustl(c_array(5)))//
     &          trim(adjustl(t2_array(3)))//trim(adjustl(c_array(3)))//
     &          trim(adjustl(t2_array(7)))//trim(adjustl(c_array(7)))//
     &          trim(adjustl(t2_array(2)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(t2_array(6)))//trim(adjustl(c_array(6)))
      case default
      item%idx2=trim(adjustl(t2_array(1)))//trim(adjustl(c_array(5)))//
     &          trim(adjustl(t2_array(2)))//trim(adjustl(c_array(6)))//
     &          trim(adjustl(t2_array(3)))//trim(adjustl(c_array(7)))//
     &          trim(adjustl(t2_array(5)))//trim(adjustl(c_array(1)))//
     &          trim(adjustl(t2_array(6)))//trim(adjustl(c_array(2)))//
     &          trim(adjustl(t2_array(7)))//trim(adjustl(c_array(3)))
      end select

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

         write(lulog,*) "sum_c1: ", sum_c1
         write(lulog,*) "sum_c2: ", sum_c1
         write(lulog,*) "sum_a1: ", sum_c1
         write(lulog,*) "sum_c1: ", sum_c1

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
      else if (item%rank3==4 .and.
     &         item%rank1==6 .or. item%rank2==6) then
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
      item%spin_idx=1

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
        write(item%logfile,*) "Error, didn't print out spin case"
      end if

      if (eloop) then
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
     &     r1,r2

       s1=.false.
       s2=.false.

      ! Determine rank; may have swapped t1 and t2 to move larger tensor
      ! to t1
      if (item%rank1==2 .and. item%rank2==4) then
         r1=item%rank2
         r2=item%rank1
      else if (item%rank1==0 .and. item%rank2==4) then
         r1=item%rank2
         r2=item%rank1
      else
         r1=item%rank1
         r2=item%rank2
      end if

       ! Pick out specific spin cases here
       if (sum(s1a)==sum(s1b) .and.
     &     sum(s2a)==sum(s2b)) then
         if (modulo(sum(s1a)+sum(s1b),2)==0 .and.
     &       modulo(sum(s2a)+sum(s2b),2)==0) then

            ! Doesn't work for rank 6 tensors yet...
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

!            DEBUG
!            if (logi) then
!               write(item%logfile,*) "s1b: ", s1b
!               write(item%logfile,*) "s1a: ", s1a
!               write(item%logfile,*)
!               write(item%logfile,*) "s2b: ", s2b
!               write(item%logfile,*) "s2a: ", s2a
!               write(item%logfile,*)
!               write(item%logfile,*)
!            else
!               write(item%logfile,*) "s1b: ", s1a
!               write(item%logfile,*) "s1a: ", s1b
!               write(item%logfile,*)
!               write(item%logfile,*) "s2b: ", s2b
!               write(item%logfile,*) "s2a: ", s2a
!               write(item%logfile,*)
!               write(item%logfile,*)
!            end if

            if (eloop==.false.) write(item%logfile,*) 'BEGIN'

            if (item%rank1==2 .and. item%rank2==4 .or.
     &          item%rank1==0 .and. item%rank2==4) then
               ! t1 and t2 were swapped in summation
               call print_itf_line(item,s2,s1)
            else 
               call print_itf_line(item,s1,s2)
            end if

            eloop=.true.
            item%spin_idx=item%spin_idx+1
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
      subroutine itf_contr_init(contr_info,itf_item,perm,antis,lulog)
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
     &     antis,
     &     lulog
      ! This should be in a 'constructor'

      ! Assign output file
      itf_item%logfile=lulog

      ! Assign permutation number
      itf_item%permute=perm
      itf_item%perm4=antis

      ! Assign labels
      itf_item%label_t1=contr_info%label_op1
      itf_item%label_t2=contr_info%label_op2
      itf_item%label_res=contr_info%label_res

      ! Assign factor
      itf_item%fact=contr_info%fact

      ! Check if an intermediate
      call check_inter(itf_item%label_t1,itf_item%inter(1))
      call check_inter(itf_item%label_t2,itf_item%inter(2))
      call check_inter(itf_item%label_res,itf_item%inter(3))

      ! Assign index string
      call assign_index(contr_info,itf_item)

      ! Assign rank
      call itf_rank(itf_item%idx1,itf_item%rank1)
      call itf_rank(itf_item%idx2,itf_item%rank2)
      call itf_rank(itf_item%idx3,itf_item%rank3)

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
!    Check if tensor is an intermediate
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
