*----------------------------------------------------------------------*
!>    returns the nroot eigenvectors of a davidson subspace
!!
!!    @param[in] dvdsbsp davidson subspace
!!    @param[out] vecs(ndim:nroot) eigenvectors
!!    @param[out] eigr(nroot) real part of eigenvalues
!!    @param[out] eigi(nroot) imaginary part of eigenvalues
!!    @param[in] nroot number or requested eigenvectors
!!    @param[in] ndim dimension of the eigenvectors has to be identical to the dimension of the subspace
*----------------------------------------------------------------------*
      subroutine dvdsbsp_get_eigenvec(dvdsbsp, vecs, eigr, eigi, nroot,
     &     ndim)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'

      integer, parameter::
     &     ntest = 000
      character(len=*),parameter::
     &     i_am="dvdsbsp_get_eigenvec"

      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp

      integer,intent(in)::
     &     nroot, ndim

      real(8),intent(inout)::
     &     eigr(nroot), eigi(nroot), vecs(ndim,nroot)

      real(8),allocatable::
     &     xmat(:),           !submatrix only
     &     xscr(:),             !scratch
     &     totvecs(:,:),        ! all eigenvectors
     &     toteigr(:),toteigi(:) ! all eigenvalues


      integer::
     &     maxdim,
     &     ii

      integer::
     &     ierr

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) "eigenvectors for ",nroot,"roots requested"
         write(lulog,*)  ndim," dimensions available"
      end if
      if(dvdsbsp%ncursub.ne.ndim) call quit(1,i_am,
     &     "requested dimension inconsistent with subspace.")

      allocate(xmat(ndim*ndim),xscr(ndim*ndim),totvecs(ndim,ndim),
     &     toteigr(ndim),toteigi(ndim) )

      maxdim=dvdsbsp%nmaxsub
      if (.not.dvdsbsp%with_metric)then
!      if(.True.)then
         if (ndim .eq. maxdim)then
            call eigen_asym(ndim, dvdsbsp%vMv_mat, toteigr, toteigi,
     &           totvecs, xscr, ierr)
         else
            call extract_submatrix(dvdsbsp%vMv_mat,  maxdim,  xmat,ndim) !vMv_mat -> xmat
            call eigen_asym(ndim, xmat, toteigr, toteigi,
     &           totvecs, xscr, ierr)
         end if
      else !with metric
         if (ndim .eq. maxdim)then
            call eigen_asym_met(ndim, dvdsbsp%vMv_mat, dvdsbsp%vSv_mat,
     &           toteigr, toteigi,
     &           totvecs, xscr, ierr)
         else
            call extract_submatrix(dvdsbsp%vMv_mat,  maxdim,  xmat,ndim) !vMv_mat -> xmat
            call extract_submatrix(dvdsbsp%vSv_mat,  maxdim,  xscr,ndim) !vSv_mat -> xscr
            call eigen_asym_met(ndim, xmat, xscr, toteigr, toteigi,
     &           totvecs, xscr, ierr) !xscr =vSv_matrix and scratch

         end if
      end if

      eigi=toteigi(1:nroot)
      eigr=toteigr(1:nroot)
      do ii=1,nroot
         if (eigi(ii).gt.0) call quit(2,i_am,
     &        "imaginary eigenvalue detected")
      end do

      if(ntest.ge.100)then
         write (lulog,*) "eigenvalues"
         do ii=1,ndim
            write (lulog,*)  toteigr(ii), toteigi(ii)
         end do
      end if
      do ii=1,nroot
        vecs(1:ndim,ii)=totvecs(1:ndim,ii)
      end do

      deallocate(xmat,xscr,totvecs,
     &     toteigr,toteigi)
      return
      contains

*----------------------------------------------------------------------*
!>   extracts th upper left  submatrix  of mat_in to mat_out
!!
!!   @param mat_in a nmat_in x nmat_in matrix
!!   @param nmat_in dimension of mat_in
!!   @param mat_out extracted nmat_out x nmat_out matrix
!!   @param nmat_out dimension of mat_out
*----------------------------------------------------------------------*
      subroutine extract_submatrix(mat_in,  nmat_in, mat_out, nmat_out)
*----------------------------------------------------------------------*
      implicit none
      character(len=*),parameter::
     &     i_am="extract_submatrix"
      integer,parameter::
     &     ntest=00

      integer,intent(in)::
     &     nmat_in,                !> actual dimension of matrix
     &     nmat_out                 !> dimension of filled part

      real(8),dimension(:),intent(in)::
     &     mat_in
      real(8),dimension(:),intent(out)::
     &     mat_out

      integer::
     &     kdx, idx, jdx

      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,i_am)

      do idx=1,nmat_out
         do jdx=1,nmat_out
            mat_out(nmat_out*(idx-1) + jdx)
     &           =mat_in(nmat_in*(idx-1) + jdx)
         end do
      end do
      return
      end subroutine



      end subroutine
