*------------------------------------------------------------------------*
      subroutine prod_pack_op(mel_in,mel_in2,mel_out,ndim)
*------------------------------------------------------------------------*
*     Multiply "packed" operators
*
*     mel_in -> The first "packed" operator
*     mel_in2-> The second "packed" operator
*     mel_out -> return the "packed" product
*
*     --------
*     Structure of a packed ME-list: the matrix is stored
*     columnwise in a sequence of records:
*
*     a11  a12  a13
*     a21  a22  a23
*     a31  a32  a33      is stored as
*
*     ndim = 3
*     rec:     1    2    3    4    5    6    7    8    8
*     elem:   a11, a21, a31, a12, a22, a32, a13, a23, a33
*     
*     yuri, set 2015
*------------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'

      type(me_list), intent(inout), target ::
     &     mel_in, mel_in2, mel_out
      integer, intent(in) ::
     &     ndim

      type(me_list), pointer ::
     &     mel
      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      logical ::
     &     closeit
      integer ::
     &     ic, il, irec, ioff, irec_ini, i, j, k, imat
      real(8), allocatable ::
     &     matrix(:,:), matrix2(:,:), matrix_out(:,:)

      allocate(matrix(ndim,ndim),
     &     matrix2(ndim,ndim),
     &     matrix_out(ndim,ndim))

      do imat=1,2
       if (imat.eq.1) then
        mel => mel_in
       else
        mel => mel_in2
       end if

       ffop => mel%fhand
       op => mel%op
       if (.not.associated(ffop))
     &      call quit(1,'prod_packed_op','No file assigned to list: '//
     &      trim(mel%label))
       if (ffop%unit.le.0) then
        call file_open(ffop)
        closeit = .true.
       else
        closeit = .false.
       end if
       if(ndim*ndim.ne.
     &      (ffop%active_records(2)-ffop%active_records(1)+1))
     &      call quit(1,'prod_packed_op',
     &      'inconsistent number of records in input list!'//
     &      ' It must be equal the square of ndim')
       
!     Get matrix
       irec_ini = ffop%current_record
       il = 1
       ic = 1
       do irec = ffop%active_records(1),ffop%active_records(2)
        call switch_mel_record(mel,irec)
        ioff = ffop%length_of_record*(ffop%current_record-1)
        if (imat.eq.1) then
         call get_vec(ffop,matrix(il:il,ic:ic),ioff+1,ioff+1)
        else
         call get_vec(ffop,matrix2(il:il,ic:ic),ioff+1,ioff+1)
        end if
        if(mod(irec - ffop%active_records(1) + 1,ndim).eq.0)then
         il = 1
         ic = ic+1
        else
         il = il+1
        end if
       end do
       call switch_mel_record(mel,irec_ini)
       if (closeit)
     &      call file_close_keep(ffop)
      end do

!     Multiply
      do i=1,ndim
       do j=1,ndim
        matrix_out(i,j) = 0.0d0
        do k=1,ndim
         matrix_out(i,j) = matrix_out(i,j) +
     &        matrix(i,k)*matrix2(k,j)
        end do
       end do
      end do

c     dbg
c      print*,"------------------------"
c      print*,"debug for prod_packed_op"
c      print*,"Matrix 1:"
c      do i=1,ndim
c       print*,matrix(i,:)
c      end do
c      print*,"Matrix 2:"
c      do i=1,ndim
c       print*,matrix2(i,:)
c      end do
c      print*,"Product:"      
c      do i=1,ndim
c       print*,matrix_out(i,:)
c      end do
c     dbg
      
! Put result
      ffop => mel_out%fhand
      op => mel_out%op
      if (.not.associated(ffop))
     &     call quit(1,'inv_packed_op','No file assigned to list: '//
     &     trim(mel_out%label))
      if (ffop%unit.le.0) then
       call file_open(ffop)
       closeit = .true.
      else
       closeit = .false.
      end if
      if(ndim*ndim.ne.
     &     (ffop%active_records(2)-ffop%active_records(1)+1))
     &     call quit(1,'inv_packed_op',
     &     'inconsistent number of records in input list!'//
     &     ' It must be equal the square of ndim')

      irec_ini = ffop%current_record
      il = 1
      ic = 1
      do irec = ffop%active_records(1),ffop%active_records(2)
       call switch_mel_record(mel_out,irec)
       ioff = ffop%length_of_record*(ffop%current_record-1)
       call put_vec(ffop,matrix_out(il:il,ic:ic),ioff+1,ioff+1)
       if(mod(irec - ffop%active_records(1) + 1,ndim).eq.0)then
        il = 1
        ic = ic+1
       else
        il = il+1
       end if
      end do
      call switch_mel_record(mel_out,irec_ini)
      if (closeit)
     &     call file_close_keep(ffop)
      
      end
