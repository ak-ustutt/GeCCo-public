*----------------------------------------------------------------------*
      subroutine optc_diis_extr(xmat,xvec,ndim,ndel,xscr,kpiv)
*----------------------------------------------------------------------*
*
*     xmat contains the DIIS B-matrix; obtain the new extrapolation
*     weights on xvec
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ndim
      integer, intent(out) ::
     &     ndel
      integer, intent(inout) ::
     &     kpiv(*)
      real(8), intent(in) ::
     &     xmat(*)
      real(8), intent(out) ::
     &     xvec(*)
      real(8), intent(inout) ::
     &     xscr(*)

      logical ::
     &     again
      integer ::
     &     ncdim, ii, jj, iioff, iioff2
      real(8) ::
     &     cond, xcorsum

      if (ntest.gt.0) then
        write(lulog,*) '===================='
        write(lulog,*) ' DIIS extrapolation' 
        write(lulog,*) '===================='
      end if

      again = .true.
      ncdim = ndim
      do while(again)
        ndel = ndim - ncdim
        do ii = 1, ncdim
          iioff  = ii*(ii-1)/2
          iioff2 = (ii+ndel)*(ii+ndel-1)/2 + ndel
          do jj = 1, ii
            xscr(iioff+jj) = xmat(iioff2+jj)
          end do
        end do
        iioff = ncdim*(ncdim+1)/2
        do jj = 1, ncdim
          xscr(iioff+jj) = -1d0
        end do
        xscr(iioff+ncdim+1) = 0d0

        if (ntest.ge.20) then
          write(lulog,*) 'Augmented matrix passed to dspco:'
          call prtrlt(xscr,ncdim+1)
        end if
        
* solve DIIS problem to obtain new weights
        ! first factorize the DIIS matrix and get its condition number
        call dspco(xscr,ncdim+1,kpiv,cond,xvec)
        if (ntest.ge.15) write(lulog,*)'ncdim,condition number: ',
     &       ncdim,cond
        if (cond.eq.0d0) call quit(1,'optc_diis_extr',
     &       'condition number = 0. maxsub > # of parameters?')
        if (ntest.ge.30) then
          write(lulog,*) 'Factorized matrix from dspco:'
          call prtrlt(xscr,ncdim+1)
          write(lulog,*) 'Pivot array:'
          call iwrtma(kpiv,1,ncdim+1,1,ncdim+1)
        end if

        ! ... we set up the RHS vector and solve the DIIS equations
        xvec(1:ncdim) =  0d0
        xvec(ncdim+1) = -1d0
        call dspsl(xscr,ncdim+1,kpiv,xvec)
        if (ntest.ge.10) then
          write(lulog,*) 'result of dspsl:'
          write(lulog,*) 'w:', xvec(1:ncdim)
          write(lulog,*) 'l:', xvec(ncdim+1)
        end if

* analyze solution and request deletion of vectors, if necessary
        xcorsum = 0d0
c        do ii = 1, ncdim-1
c          xcorsum = xcorsum + abs(xvec(ii))
c        end do
c        if (xcorsum/abs(xvec(ncdim)).gt.1.2d0) then
c          ncdim = ncdim -1
c        else
          again = .false.
c        end if

      end do

      if (ncdim.lt.ndim) then
        ndel = ndim-ncdim
        do ii = ncdim, 1, -1
          xvec(ii+ndel) = xvec(ii)
        end do
        xvec(1:ndel) = 0d0
      end if
        
      if (ntest.ge.10) then
        write(lulog,*) 'final result for extrapolation weights:'
        call wrtmat2(xvec,1,ndim,1,ndim)
      end if

      return
      end
