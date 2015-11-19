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
     &     ic, il, irec, ioff, irec_ini
      real(8), allocatable ::
     &     M_ini(:,:), M_inv(:,:),M(:,:)

      real(8) ::
     &     det, A, B, C, D, E, F, G, H, I

      real(8), parameter ::
     &     zero_thr = 1d-10

      allocate(M(ndim,ndim),
     &     M_ini(ndim,ndim),
     &     M_inv(ndim,ndim))

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

!--------------------
!     Get matrix
      irec_ini = ffop%current_record
      il = 1
      ic = 1
      do irec = ffop%active_records(1),ffop%active_records(2)
       call switch_mel_record(mel_in,irec)
       ioff = ffop%length_of_record*(ffop%current_record-1)
       call get_vec(ffop,M_ini(il:il,ic:ic),ioff+1,ioff+1)
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


! ndim = 2, 3 can be directly
      if (ndim.EQ.2) then
       det = M_ini(1,1)*M_ini(2,2) - M_ini(1,2)*M_ini(2,1)

       M_inv(1,1) =  M_ini(2,2)/det
       M_inv(2,1) = -M_ini(1,2)/det
       M_inv(1,2) = -M_ini(2,1)/det
       M_inv(2,2) =  M_ini(1,1)/det

      elseif (ndim.EQ.3) then

       A = M_ini(2,2)*M_ini(3,3) - M_ini(2,3)*M_ini(3,2)
       B = M_ini(2,1)*M_ini(3,3) - M_ini(2,3)*M_ini(3,1)
       C = M_ini(2,1)*M_ini(3,2) - M_ini(2,2)*M_ini(3,1)

       D = M_ini(1,2)*M_ini(3,3) - M_ini(1,3)*M_ini(3,2)
       E = M_ini(1,1)*M_ini(3,3) - M_ini(1,3)*M_ini(3,1)
       F = M_ini(1,1)*M_ini(3,2) - M_ini(1,2)*M_ini(3,1)

       G = M_ini(1,2)*M_ini(2,3) - M_ini(1,3)*M_ini(2,2)
       H = M_ini(1,1)*M_ini(2,3) - M_ini(1,3)*M_ini(2,1)
       I = M_ini(1,1)*M_ini(2,2) - M_ini(1,2)*M_ini(2,1)

       B = -B
       D = -D
       F = -F
       H = -H

       det = M_ini(1,1)*A + M_ini(1,2)*B + M_ini(1,3)*C

       M_inv(1,1) = A/det
       M_inv(2,1) = B/det
       M_inv(3,1) = C/det

       M_inv(1,2) = D/det
       M_inv(2,2) = E/det
       M_inv(3,2) = F/det

       M_inv(1,3) = G/det
       M_inv(2,3) = H/det
       M_inv(3,3) = I/det

      else

       M = M_ini
       call gaussj(M, ndim, ndim)
       M_inv = M

      end if

c     dbg
c      call check_inv(ndim, M_ini, M_inv)

!--------------------
!     Put inverse
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
       call put_vec(ffop,M_inv(il:il,ic:ic),ioff+1,ioff+1)
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
