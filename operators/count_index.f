*----------------------------------------------------------------------*
      subroutine count_index(lulog,idx,iocc,irstr,ngas,nspin,nops)
*----------------------------------------------------------------------*
!     Count index of operator
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     lulog, iocc(ngastp,2), ngas, nspin,
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
      subroutine remove_whitespace(string)
*----------------------------------------------------------------------*
!     Remove whitespace inbetween ITF index string
*----------------------------------------------------------------------*

      implicit none

      character(len=4), intent(inout) ::
     &     string   ! IFT index string with whitespaces
      character(len=4) ::
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
      subroutine print_itf_line(t1, t2, idx1, idx2,
     &                          idx3, inter, lulog)
*----------------------------------------------------------------------*
!     Print line of ITF code
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      character(len=maxlen_bc_label) ::
     &     t1, t2      ! Name of tensors involved in the contraction
      character(len=4), intent(in) ::
     &     idx1, idx2, idx3   ! Index strings
      integer, intent(in) ::
     &     inter,                ! Number of current intermediate
     &     lulog                 ! File to write to
      character(len=5) ::
     &     s_int,                ! Intermdiate tensor number
     &     p_int                 ! Previous intermdiate tensor number
      character(len=50) ::
     &     itf_line              ! Line of ITF code

      write(s_int,'(i0)') inter

      write(lulog,*) 'TENSOR:'
      if (t1.eq.'_STIN0001') then
        ! Previous intermediate appears on rhs
        t1='I'
        write(p_int,'(i0)') inter-1
        itf_line='I'//trim(s_int)//'['//trim(idx3)//']+='//
     &      trim(t1)//trim(p_int)//'['//trim(idx1)//']'//
     &      trim(t2)//'['//trim(idx2)//']'
        write(lulog,*) trim(itf_line)
      else
        itf_line='I'//trim(s_int)//'['//trim(idx3)//']+='//
     &      trim(t1)//'['//trim(idx1)//']'//trim(t2)//'['//
     &      trim(idx2)//']'
        write(lulog,*) trim(itf_line)
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine assign_index(arr,nops)
*----------------------------------------------------------------------*
!     Assign index letters to arrays
*----------------------------------------------------------------------*

      implicit none
      
      character, intent(inout)  ::
     &     arr(2,2)        ! Array to assign letters
      integer, intent(in) ::
     &     nops(4,2)      ! Matrix of index info
      integer ::
     &     i,j
      character, dimension(4) ::
     &     hol=(/ 'i','j','k','l' /),
     &     par=(/ 'a','b','c','d' /),
     &     val=(/ 'u','v','w','x' /)
      

      ! Annhilators first
      do i=1, nops(1,2)
          arr(1,i)=hol(i)
      end do
      do i=1, nops(2,2)
        arr(1,i+nops(1,2))=par(i)
      end do
      do i=1, nops(3,2)
        arr(1,i+nops(1,2)+nops(2,2))=val(i)
      end do
      ! Creations second
      do i=1, nops(1,1)
        arr(2,i)=hol(i+nops(1,2))
      end do
      do i=1, nops(2,1)
        arr(2,i+nops(1,1))=par(i+nops(2,2))
      end do
      do i=1, nops(3,1)
        arr(2,i+nops(1,1)+nops(2,1))=val(i+nops(3,2))
      end do

      return
      end


*----------------------------------------------------------------------*
      subroutine index_array(lulog,tensor,i_array,j_array,k_array)
*----------------------------------------------------------------------*
!     Determine tensor index for ITF
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lulog
      type(binary_contr), intent(in) ::
     &     tensor
      character, intent(inout)  ::
     &     i_array(2,2),
     &     j_array(2,2),
     &     k_array(2,2)
      integer ::
     &     nops(4,2),     ! Matrix of index info
     &     eops(4,2),     ! Matrix of external index info
     &     eops1(4,2),    ! Matrix of external index info, tensor 1
     &     shift_l,       ! Letter shift when assging index of second tensor
     &     counter1, counter2,
     &     i,j
      character, dimension(4) ::
     &     hol=(/ 'i','j','k','l' /),
     &     par=(/ 'a','b','c','d' /),
     &     val=(/ 'u','v','w','x' /)
      character ::
     &     tmp1(2,2),
     &     tmp2(2,2)


      ! Count index of first tensor
      do i = 1, tensor%nj_op1
        call count_index(lulog,i,
     &       tensor%occ_op1(1:,1:,i),
     &       tensor%rst_op1(1:,1:,1:,1:,1:,i),tensor%ngas,
     &       tensor%nspin,nops)
      end do
       
      call assign_index(i_array,nops)

      ! Count contraction index
      do i = 1, 1 !Just get one block
        call count_index(lulog,i,
     &       tensor%occ_cnt(1:,1:,i),
     &       tensor%rst_cnt(1:,1:,1:,1:,1:,i),
     &       tensor%ngas,tensor%nspin,nops)
       end do
      
      call assign_index(j_array,nops)

      ! Count external index of second tensor
      do i = 1, 1 ! Just the first block
        call count_index(lulog,i,
     &     tensor%occ_ex2(1:,1:,i),
     &     tensor%rst_ex2(1:,1:,1:,1:,1:,i),
     &     tensor%ngas,tensor%nspin,eops)
      end do
      
      ! Count external index of first tensor
      do i = 1, 1 ! Just the first block
        call count_index(lulog,i,
     &     tensor%occ_ex1(1:,1:,i),
     &     tensor%rst_ex1(1:,1:,1:,1:,1:,i),
     &     tensor%ngas,tensor%nspin,eops1)
      end do
      
      ! Creations first
      ! Each index must be shifted so as a) not to overwrite previous index
      !                                  b) not to use same index twice
      ! Using 'transpose' of E (so no. of ann + cre !> 3
      ! shift_l shifts position in letter array, so external indicies are unique
      shift_l=0
      do i=1, eops(1,1)
        shift_l=i+nops(1,2)+nops(1,1)+eops1(1,2)+eops(1,1)
        j_array(1,i+nops(1,2)+nops(2,2)+
     &  nops(3,2))=hol(shift_l)
      end do
      do i=1, eops(2,1)
        shift_l=i+nops(2,2)+nops(2,1)+eops1(2,2)+eops1(2,1)
        j_array(1,i+nops(1,2)+nops(2,2)+nops(3,2)+
     &  eops(1,1))=par(shift_l)
      end do
      do i=1, eops(3,1)
        shift_l=i+nops(3,2)+nops(3,1)+eops1(3,2)+eops1(3,1)
        j_array(1,i+nops(1,2)+nops(2,2)+nops(3,2)+eops(1,1)+
     &  eops(2,1))=val(shift_l)
      end do
      ! Annhilators second
      do i=1, eops(1,2)
        shift_l=i+nops(1,2)+nops(1,1)+eops(1,1)+eops1(1,2)+eops1(1,1)
        j_array(2,i+nops(1,1)+nops(2,1)+
     &  nops(3,1))=hol(shift_l)
      end do
      do i=1, eops(2,2)
        shift_l=i+nops(2,2)+nops(2,1)+eops(2,1)+eops1(2,2)+eops1(2,1)
        j_array(2,i+nops(1,1)+nops(2,1)+nops(3,1)+
     &  eops(1,2))=par(shift_l)
      end do
      do i=1, eops(3,2)
        shift_l=i+nops(3,2)+nops(3,1)+eops(3,1)+eops1(3,2)+eops1(3,2)
        j_array(2,i+nops(1,1)+nops(2,1)+nops(3,1)+eops(1,2)+
     &  eops(2,2))=val(shift_l)
      end do


      ! Compare index of tensors, different is used to label result tensor
      do i=1, 2
        do j=1, 2
          tmp1(i,j)='>'
          tmp2(i,j)='>'
        end do
      end do

      
      do i=1, 2
        do j=1, 2
          if (i_array(i,j).ne.j_array(1,1)) then
            if (i_array(i,j).ne.j_array(1,2)) then
              if (i_array(i,j).ne.j_array(2,1)) then
                if (i_array(i,j).ne.j_array(2,2)) then
                  tmp1(i,j)=i_array(i,j)
                end if
              end if
            end if
          end if
          if (j_array(i,j).ne.i_array(1,1)) then
            if (j_array(i,j).ne.i_array(1,2)) then
              if (j_array(i,j).ne.i_array(2,1)) then
                if (j_array(i,j).ne.i_array(2,2)) then
                  tmp2(i,j)=j_array(i,j)
                end if
              end if
            end if
          end if
        end do
      end do

      counter1=1
      counter2=1
      do i=1, 2
        if(tmp1(1,i).ne.'>') then
          k_array(1,counter1)=tmp1(1,i)
          counter1=counter1+1
        endif
        if(tmp1(2,i).ne.'>') then
          k_array(2,counter2)=tmp1(2,i)
          counter2=counter2+1
        endif
      end do
      do i=1, 2
        if(tmp2(1,i).ne.'>') then
          k_array(1,counter1)=tmp2(1,i)
          counter1=counter1+1
        endif
        if(tmp2(2,i).ne.'>') then
          k_array(2,counter2)=tmp2(2,i)
          counter2=counter2+1
        endif
      end do

      return
      end
