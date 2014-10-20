*------------------------------------------------------------------------*
      subroutine diag_packed_op(mel_in,mel_evec,mel_eval,ndim,orb_info)
*------------------------------------------------------------------------*
*     Diagonalise an the "packed" operator
*
*     mel_in -> contains the "packed" operator: each entries of
*               the operator is in a different record, in a total of
*               ndim x ndim records
*     mel_evec -> resulting eigenvectors, in ndim records
*     mel_eval -> resulting eigenvalues, in ndim records
*     
*     yuri, oct 2014
*------------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'

      type(me_list), intent(inout) ::
     &     mel_in, mel_evec, mel_eval
      integer, intent(in) ::
     &     ndim
      type(orbinf) ::
     &     orb_info

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      logical ::
     &     closeit
      integer ::
     &     ic, il, irec, ioff, irec_ini, ierr
      real(8), allocatable ::
     &     matrix(:,:), matrix_ini(:,:), eigenvec(:,:),
     &     eigenval_r(:), eigenval_i(:),
     &     xscr(:,:)

      real(8), parameter ::
     &     zero_thr = 1d-10
      character(50) ::
     &     out_format

      allocate(matrix(ndim,ndim),
     &     matrix_ini(ndim,ndim),
     &     eigenvec(ndim,ndim),
     &     eigenval_r(ndim),
     &     eigenval_i(ndim),
     &     xscr(ndim,ndim))

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

      ! Solve the EVP
      call eigen_asym(ndim,matrix,eigenval_r,eigenval_i,eigenvec,xscr,
     &     ierr)

      if(ierr.ne.0)
     &     call quit(1,'diag_packed_op','Error in eigen_asym!')
      if(any([(abs(eigenval_i(ic)).gt.zero_thr,ic=1,ndim)]))
     &     call quit(1,'diag_packed_op','Non real eigenvalue!')
c dbg
c      call check_diag(ndim,ndim,matrix_ini,eigenvec,eigenval_r)
c dbgend

      ! A nice output, specific for the multistate problem
      write(out_format,fmt='(A,i0,A)')
     &     '(">>>",x,f24.12,',ndim,'(x,f12.6))'
      write(lulog,'(">>> Multistate energies and eigenvectors:")')
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
