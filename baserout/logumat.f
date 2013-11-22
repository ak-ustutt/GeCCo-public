*----------------------------------------------------------------------*
      subroutine logumat(ndim,xlogx,xmat,xscr1,xscr2,xscr3)
*----------------------------------------------------------------------*
*     calculate the logarithm of a unitary matrix
*
*     the algorithm will use the eispack-routine rg() to calculate the
*     eigenvalues of the matrix which are decomposed into modulus and 
*     angle.
*     the modulus should be one always, else the routine exits.
*
*     andreas, aug 2004
*     
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter :: ntest = 00, maxn = 100

      integer, intent(in) ::
     &     ndim
      real(8), intent(inout) ::
     &     xlogx(ndim,ndim), xmat(ndim,ndim),
     &     xscr1(ndim,ndim), xscr2(ndim,ndim), xscr3(ndim,ndim)

      integer ::
     &     ii, ierr
      real(8) ::
     &     ang, xmod     ,ang1,ang2

* O(N) scratch
      real(8) ::
     &     eigr(ndim), eigi(ndim), scr(ndim)
      integer ::
     &     iscr(ndim)


      if (ntest.gt.0) then
        write(lulog,*) ' ==================== '
        write(lulog,*) '  LOGUMAT at work !!  '
        write(lulog,*) ' ==================== '

      end if

      if (ntest.ge.100) then
        write(lulog,*) ' xmat on entry:'
        call wrtmat2(xmat,ndim,ndim,ndim,ndim)
      end if

      xscr2(1:ndim,1:ndim) = xmat(1:ndim,1:ndim)

      ! get eigenvalues and -vectors ...
      call rg(ndim,ndim,xscr2,
     &        eigr,eigi,1,xscr1,
     &        iscr,scr,ierr)
      ! and normalize vectors (not done by rg)
      call nrmvec(ndim,xscr1,eigi)

      if (ierr.ne.0) then
        write(lulog,*) 'error code from rg: ',ierr
        stop 'logumat (1)'
      end if

      if (ntest.ge.10) write(lulog,*) ' eigenvalues of matrix:'

*----------------------------------------------------------------------*
*     the eigenvalues are v = exp(a+ib) so the logarithm log(v) yields
*     a and b. as the matrix is unitary, a is always 0 and we are left
*     with b, which is the angle in the complex plane.
*     the angles will be collected in eigr(), later referred to as 
*     matrix D
*----------------------------------------------------------------------*
      ierr = 0
      if (ntest.ge.10) write(lulog,'(x,a,g10.3)') 
     &                      "threshold = ",1d4*epsilon(1d0)
      do ii = 1, ndim
        xmod = eigr(ii)*eigr(ii) + eigi(ii)*eigi(ii)
        if (abs(xmod-1d0).gt.1d4*epsilon(1d0)) ierr = ierr+1
        ang1 = atan2(eigi(ii),eigr(ii))
c        ang2 = acos(eigr(ii))*sign(1d0,eigi(ii))
        if (ntest.ge.10)
     &       write(lulog,'(i4,2(2x,e20.10),3(2x,f15.10),g10.3)')
     &       ii,eigr(ii),eigi(ii),xmod,ang1,ang2,xmod-1d0
        eigr(ii) = ang1
      end do

      if (ierr.gt.0) then
        write(lulog,*) 'error: detected eigenvalues with |v| != 1'
        stop 'logumat (2)'
      end if

      ! sort components of transformation matrix into 
      ! real and imaginary part: U = A + iB
      !  A on xscr1
      !  B on xscr2

      ii = 0
      do while(ii.lt.ndim)
        ii = ii+1
        if (eigi(ii).eq.0d0) then ! real eigenvalue
          xscr2(1:ndim,ii) = 0d0
        else ! complex pair
          xscr2(1:ndim,ii) = xscr1(1:ndim,ii+1)  ! imag. part
          xscr2(1:ndim,ii+1) = -xscr2(1:ndim,ii) ! and cmplx. conj.
          xscr1(1:ndim,ii+1) = xscr1(1:ndim,ii)
          ii = ii+1                              ! add. increment
        end if
      end do

      if (ntest.ge.100) then
        write(lulog,*) ' eigenvectors (Re):'
        call wrtmat2(xscr1,ndim,ndim,ndim,ndim)
        write(lulog,*) ' eigenvectors (Im):'
        call wrtmat2(xscr1,ndim,ndim,ndim,ndim)
      end if

*----------------------------------------------------------------------*
*
*     now we have to calculate U iD U^+
*     
*     the real part is
*                A iD (iB)^+ + (iB) iD A^+ = A D B^T - B D A^T
*
*     the imaginary part is  
*                A iD A^+ + iB iD (iB)^+ = i (A D A^T + B D B^T)
*     
*     as iD has either zero or pairwise conjugate entries, the imaginary
*     part vanishes (note that we started from a real unitary matrix)
*
*----------------------------------------------------------------------*
      
      ! A on xscr1
      ! B on xscr2

      ! A D --> xscr3
      do ii = 1, ndim
        xscr3(1:ndim,ii) = xscr1(1:ndim,ii)*eigr(ii)
      end do

      ! AD B^T --> xlogx
c      call matml7(xlogx,xscr3,xscr2,
c     &            ndim,ndim,
c     &            ndim,ndim,
c     &            ndim,ndim,
c     &            0d0,1d0, 2 )
      call dgemm('n','t',ndim,ndim,ndim,
     &           1d0,xscr3,ndim,
     &               xscr2,ndim,
     &           0d0,xlogx,ndim)

      ! B D --> xscr3
      do ii = 1, ndim
        xscr3(1:ndim,ii) = xscr2(1:ndim,ii)*eigr(ii)
      end do

      !-BD A^T --> xlogx
c      call matml7(xlogx,xscr3,xscr1,
c     &            ndim,ndim,
c     &            ndim,ndim,
c     &            ndim,ndim,
c     &            1d0,-1d0, 2 )
      call dgemm('n','t',ndim,ndim,ndim,
     &          -1d0,xscr3,ndim,
     &               xscr1,ndim,
     &           0d0,xlogx,ndim)

      if (ntest.ge.100) then
        write(lulog,*) ' result on xlogx:'
        call wrtmat2(xlogx,ndim,ndim,ndim,ndim)
      end if

      return

      end
