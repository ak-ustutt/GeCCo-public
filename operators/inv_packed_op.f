*------------------------------------------------------------------------*
      subroutine inv_pack_op(mel_in,mel_out,ndim)
*------------------------------------------------------------------------*
*     Diagonalise a "packed" operator
*
*     mel_in -> contains a "packed" operator of dimension ndim
*     mel_out -> return the "packed" inverse
*
*     See file prod_packed_op.f for a description on how the elements
*     are packed.
*     
*     yuri, set 2015
*------------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'

      type(me_list), intent(inout) ::
     &     mel_in, mel_out
      integer, intent(in) ::
     &     ndim

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      logical ::
     &     closeit
      integer ::
     &     ic, il, irec, ioff, irec_ini, ierr, i, j, k
      real(8), allocatable ::
     &     matrix_ini(:,:), matrix_inv(:,:),matrix(:,:),
     &     eigenvec(:,:),
     &     eigenval_r(:), eigenval_i(:),
     &     xscr(:,:)

      real(8), parameter ::
     &     zero_thr = 1d-10

      allocate(matrix(ndim,ndim),
     &     matrix_ini(ndim,ndim),
     &     matrix_inv(ndim,ndim),
     &     eigenvec(ndim,ndim),
     &     eigenval_r(ndim),
     &     eigenval_i(ndim),
     &     xscr(ndim,ndim))

      ffop => mel_in%fhand
      op => mel_in%op
      if (.not.associated(ffop))
     &     call quit(1,'inv_packed_op','No file assigned to list: '//
     &     trim(mel_in%label))
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
     &     ' It must be equal the square ndim')

      ! Get matrix
      irec_ini = ffop%current_record
      il = 1
      ic = 1
      do irec = ffop%active_records(1),ffop%active_records(2)
       call switch_mel_record(mel_in,irec)
       ioff = ffop%length_of_record*(ffop%current_record-1)
       call get_vec(ffop,matrix_ini(il:il,ic:ic),ioff+1,ioff+1)
       if(mod(irec - ffop%active_records(1) + 1,ndim).eq.0)then
        il = 1
        ic = ic+1
       else
        il = il+1
       end if
      end do
      call switch_mel_record(mel_in,irec_ini)
      if (closeit)
     &     call file_close_keep(ffop)

      matrix = matrix_ini

      ! Solve the EVP
      call eigen_asym(ndim,matrix,eigenval_r,eigenval_i,
     &     eigenvec,xscr,ierr)

      if(ierr.ne.0)
     &     call quit(1,'inv_packed_op','Error in eigen_asym!')
      if(any([(abs(eigenval_i(ic)).gt.zero_thr,ic=1,ndim)]))
     &     call quit(1,'inv_packed_op','Non real eigenvalue!')
      if(any([(abs(eigenval_r(ic)).lt.zero_thr,ic=1,ndim)]))
     &     call quit(1,'inv_packed_op','Zero eigenvalue!')

! Get the inverse: U D^-1 U^T
      do i=1,ndim
       do j=1,ndim
        matrix_inv(i,j) = 0.0d0
        do k=1,ndim
         matrix_inv(i,j) = matrix_inv(i,j) +
     &        eigenvec(i,k)*eigenvec(j,k)/eigenval_r(k)
        end do
       end do
      end do

cdbg
c      call check_inv(ndim, matrix_ini, matrix_inv)

! Put inverse
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
       call put_vec(ffop,matrix_inv(il:il,ic:ic),ioff+1,ioff+1)
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
