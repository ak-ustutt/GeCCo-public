*------------------------------------------------------------------------*
      subroutine diag_packed_op(mel_in,mel_evec,mel_eval,ndim,mel_S,
     &     verbose)
*------------------------------------------------------------------------*
*     Diagonalise an the "packed" operator
*
*     mel_in -> contains the "packed" operator: each entries of
*               the operator is in a different record, in a total of
*               ndim x ndim records
*     mel_evec -> resulting eigenvectors, in ndim records
*     mel_eval -> resulting eigenvalues, in ndim records
*     mel_S    -> overlap matrix, in ndim x ndim records
*     
*     yuri, oct 2014
*------------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'

      interface
       subroutine check_diag(ndim,nvec,mat,eigenvec,eigenval,s_mat)
       implicit none
       integer, intent(in) ::
     &      ndim, nvec
       real(8), intent(in) ::
     &      mat(ndim,ndim)
       real(8), intent(in) ::
     &      eigenvec(ndim,nvec), eigenval(nvec)
       real(8), intent(in), optional ::
     &      S_mat(ndim,ndim)
       end subroutine check_diag
      end interface

      type(me_list), intent(inout) ::
     &     mel_in, mel_evec, mel_eval
      integer, intent(in) ::
     &     ndim
      type(me_list), intent(inout), optional ::
     &     mel_S
      logical, intent(in), optional ::
     &     verbose

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      logical ::
     &     closeit
      integer ::
     &     ic, il, irec, ioff, irec_ini, ierr, i, j, k, l
      real(8), allocatable ::
     &     matrix(:,:), matrix_ini(:,:),
     &     S_matrix(:,:), S_matrix_ini(:,:), A_matrix(:,:),
     &     eigenvec(:,:), eigenvec_ini(:,:),
     &     eigenval_r(:), eigenval_i(:),
     &     xscr(:,:)

      real(8), parameter ::
     &     zero_thr = 1d-10
      character(50) ::
     &     out_format
      character(3) ::
     &     indicator

c dbg
c      real(8) :: K_delta
c dbg

      allocate(matrix(ndim,ndim),
     &     matrix_ini(ndim,ndim),
     &     S_matrix(ndim,ndim),
     &     S_matrix_ini(ndim,ndim),
     &     A_matrix(ndim,ndim),
     &     eigenvec(ndim,ndim),
     &     eigenvec_ini(ndim,ndim),
     &     eigenval_r(ndim),
     &     eigenval_i(ndim),
     &     xscr(ndim,ndim))

      if (present(mel_S)) then
       ffop => mel_S%fhand
       op => mel_S%op
       if (.not.associated(ffop))
     &      call quit(1,'diag_packed_op','No file assigned to list: '//
     &      trim(mel_S%label))
       if (ffop%unit.le.0) then
        call file_open(ffop)
        closeit = .true.
       else
        closeit = .false.
       end if
       if(ndim*ndim.ne.
     &      (ffop%active_records(2)-ffop%active_records(1)+1))
     &      call quit(1,'diag_packed_op',
     &      'inconsistent number of records in the S-matrix list!'//
     &      ' It must be equal the square of the number or states')

       ! Get S matrix
       irec_ini = ffop%current_record
       il = 1
       ic = 1
       do irec = ffop%active_records(1),ffop%active_records(2)
        call switch_mel_record(mel_S,irec)
        ioff = ffop%length_of_record*(ffop%current_record-1)
        call get_vec(ffop,S_matrix(il:il,ic:ic),ioff+1,ioff+1)
        if(mod(irec - ffop%active_records(1) + 1,ndim).eq.0)then
         il = 1
         ic = ic+1
        else
         il = il+1
        end if
       end do
       call switch_mel_record(mel_S,irec_ini)
       if (closeit)
     &      call file_close_keep(ffop)

       S_matrix_ini = S_matrix
       ! Solve the S-EVP
       call eigen_asym(ndim,S_matrix,eigenval_r,eigenval_i,A_matrix,
     &      xscr, ierr)
c     dbg
c       call check_diag(ndim,ndim,S_matrix_ini,A_matrix,eigenval_r)
c     dbgend
       if(ierr.ne.0)
     &      call quit(1,'diag_packed_op',
     &      'Error in eigen_asym - Smatrix!')
       if(any([(abs(eigenval_i(ic)).gt.zero_thr,ic=1,ndim)]))
     &      call quit(1,'diag_packed_op',
     &      'Non real eigenvalue in S matrix!')
       if(any([(eigenval_r(ic).lt.zero_thr,ic=1,ndim)]))
     &      call quit(1,'diag_packed_op',
     &      'Non positive definite S matrix!')

       do ic=1,ndim
        A_matrix(:,ic) = A_matrix(:,ic)/sqrt(eigenval_r(ic))
       end do
      end if

c dbg
c      print*, "Check: It must be a Kronecker delta:"
c      do i=1,ndim
c       do j=1,ndim
c        K_delta = 0.0d0
c        do k=1,ndim
c         do l=1,ndim
c          K_delta = K_delta +
c     &         A_matrix(k,i)*S_matrix_ini(k,l)*A_matrix(l,j)
c         end do
c        end do
c        print*,"Kronecker delta (",i,j,") :", K_delta
c       end do
c      end do
c dbg

      ffop => mel_in%fhand
      op => mel_in%op
      if (.not.associated(ffop))
     &     call quit(1,'diag_packed_op','No file assigned to list: '//
     &     trim(mel_in%label))
      if (ffop%unit.le.0) then
       call file_open(ffop)
       closeit = .true.
      else
       closeit = .false.
      end if
      if(ndim*ndim.ne.
     &     (ffop%active_records(2)-ffop%active_records(1)+1))
     &     call quit(1,'diag_packed_op',
     &     'inconsistent number of records in input list!'//
     &     ' It must be equal the square of the number or states')

      ! Get matrix
      irec_ini = ffop%current_record
      il = 1
      ic = 1
      do irec = ffop%active_records(1),ffop%active_records(2)
       call switch_mel_record(mel_in,irec)
       ioff = ffop%length_of_record*(ffop%current_record-1)
       call get_vec(ffop,matrix(il:il,ic:ic),ioff+1,ioff+1)
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

      matrix_ini = matrix

!     matrix = A_matrix^+ x matrix x A_matrix
      if (present(mel_S)) then
       do i=1,ndim
        do j=1,ndim
         matrix(i,j) = 0.0d0
         do k=1,ndim
          do l=1,ndim
           matrix(i,j) = matrix(i,j) +
     &          A_matrix(k,i)*matrix_ini(k,l)*A_matrix(l,j)
          end do
         end do
        end do
       end do
      end if

      ! Solve the EVP
      call eigen_asym(ndim,matrix,eigenval_r,eigenval_i,eigenvec,xscr,
     &     ierr)

      if(ierr.ne.0)
     &     call quit(1,'diag_packed_op','Error in eigen_asym!')
      if(any([(abs(eigenval_i(ic)).gt.zero_thr,ic=1,ndim)]))
     &     call quit(1,'diag_packed_op','Non real eigenvalue!')

! A_matrix x eigenvec
      if (present(mel_S)) then
       eigenvec_ini = eigenvec
       do i=1,ndim
        do j=1,ndim
         eigenvec(i,j) = 0.0d0
         do k=1,ndim
          eigenvec(i,j) = eigenvec(i,j) +
     &         A_matrix(i,k)*eigenvec_ini(k,j)
         end do
        end do
       end do
      end if
c dbg
c      if (present(mel_S)) then
c       call check_diag(ndim,ndim,matrix_ini,eigenvec,eigenval_r,
c     &      S_matrix_ini)
c      else
c       call check_diag(ndim,ndim,matrix_ini,eigenvec,eigenval_r)
c      end if
c dbgend

      ! A nice output, specific for the multistate problem
      indicator = ">>>"
      if (present(verbose)) then
       if(.not.verbose) indicator = "   "
      end if

      write(out_format,fmt='(A,i0,A)')
     &     '("'//indicator//'",x,f24.12,',ndim,'(x,f12.6))'
      if (present(mel_S)) then
       write(lulog,'("'//indicator//
     &      ' Multistate energies and eigenvectors ",'//
     &      '"(with overlap):")')
      else
       write(lulog,'("'//indicator//
     &      ' Multistate energies and eigenvectors:")')
      endif
      do ic = 1,ndim
       write(lulog,out_format)
     &      eigenval_r(ic),eigenvec(:,ic)
      end do

      ! Put eigenvectors
      ffop => mel_evec%fhand
      op => mel_evec%op
      if (.not.associated(ffop))
     &     call quit(1,'diag_packed_op','No file assigned to list: '//
     &     trim(mel_evec%label))
      if (ffop%unit.le.0) then
       call file_open(ffop)
       closeit = .true.
      else
       closeit = .false.
      end if
      if(ndim.ne.
     &     (ffop%active_records(2)-ffop%active_records(1)+1))
     &     call quit(1,'diag_packed_op',
     &     'inconsistent number of records in eigenvector list!'//
     &     ' It must be equal the number or states')

      irec_ini = ffop%current_record
      do irec = ffop%active_records(1),ffop%active_records(2)
       ic = irec - ffop%active_records(1) + 1
       call switch_mel_record(mel_evec,irec)
       ioff = ffop%length_of_record*(ffop%current_record-1)
       call put_vec(ffop,eigenvec(:,ic),ioff+1,ioff+ndim)
      end do
      call switch_mel_record(mel_evec,irec_ini)
      if (closeit)
     &     call file_close_keep(ffop)

      ! Put eigenvalues
      ffop => mel_eval%fhand
      op => mel_eval%op
      if (.not.associated(ffop))
     &     call quit(1,'diag_packed_op','No file assigned to list: '//
     &     trim(mel_eval%label))
      if (ffop%unit.le.0) then
       call file_open(ffop)
       closeit = .true.
      else
       closeit = .false.
      end if
      if(ndim.ne.
     &     (ffop%active_records(2)-ffop%active_records(1)+1))
     &     call quit(1,'diag_packed_op',
     &     'inconsistent number of records in eigenvector list!'//
     &     ' It must be equal the number or states')

      irec_ini = ffop%current_record
      do irec = ffop%active_records(1),ffop%active_records(2)
       ic = irec - ffop%active_records(1) + 1
       call switch_mel_record(mel_eval,irec)
       ioff = ffop%length_of_record*(ffop%current_record-1)
       call put_vec(ffop,eigenval_r(ic:ic),ioff+1,ioff+1)
      end do
      call switch_mel_record(mel_eval,irec_ini)
      if (closeit)
     &     call file_close_keep(ffop)

      end
