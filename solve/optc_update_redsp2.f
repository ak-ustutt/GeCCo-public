*----------------------------------------------------------------------*
      subroutine optc_update_redsp2(xmat,ndim1,ndim2,mxdim,shift,
     &     iord1,ff_sbsp1,iord2,ff_sbsp2,
     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
*----------------------------------------------------------------------*
*
*     update reduced-space matrix xmat(i,j) = <v(i)|w(j)>, where
*     xmat has the leading dimension mxdim
*     
*     for incore-version: most recent vector v expected on xbuf1
*                         most recent vector w expected on xbuf2
*
*     shift==.true.  xmat(i,j) := xmat(i+1,j+1), xmat(1,j) is discarded
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 100

      type(filinf) ::
     &     ff_sbsp1, ff_sbsp2
      logical, intent(in) ::
     &     shift
      integer, intent(in) ::
     &     ndim1, ndim2, mxdim, nincore, nwfpar, lenbuf,
     &     iord1(*), iord2(*) 
      real(8), intent(inout) ::
     &     xmat(mxdim,*), xbuf1(*), xbuf2(*), xbuf3(*)
      
      integer ::
     &     iioff, iioff2, ii, jj, irec, jrec1, jrec2

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     ddot, da_ddot

      if (ntest.gt.0) then
        write(luout,*) '==============================='
        write(luout,*) ' welcome to optc_update_redsp2'
        write(luout,*) '==============================='
      end if

      if (ntest.ge.20) then
        write(luout,*) 'subspace matrix on input:'
        call wrtmat2(xmat,ndim1,ndim1,mxdim,mxdim)
      end if

      ! subspace is full: delete contribution from first vector
      ! i.e. shift all matrix-elements down by one
      if (shift) then
        do ii = 1, ndim1 - 1
          do jj = 1, ndim1 - 1
            xmat(jj,ii) = xmat(jj+1,ii+1)
          end do
        end do
        if (ntest.ge.20) then
          write(luout,*) 'subspace matrix after shift:'
          call wrtmat2(xmat,ndim1,ndim1,mxdim,mxdim)
        end if
      end if

      ! a) update diagonal element
      jrec1 = ioptc_get_sbsp_rec(ndim1,iord1,ndim1,mxdim)
      jrec2 = ioptc_get_sbsp_rec(ndim1,iord2,ndim2,mxdim)
      if (nincore.ge.2) then
        xmat(ndim1,ndim1) = ddot(nwfpar,xbuf1,1,xbuf2,1)
      else
        xmat(ndim1,ndim1) = da_ddot(ff_sbsp1,jrec1,ff_sbsp2,jrec2,
     &       nwfpar,xbuf1,xbuf2,lenbuf)
      end if

      ! b) update last column
      do irec = 1, ndim1

        jj = iord2(irec)

        if (jj.eq.ndim1) cycle

        if (nincore.eq.3) then
          call vec_from_da(ff_sbsp2,irec,xbuf3,nwfpar)
          xmat(ndim1,jj) = ddot(nwfpar,xbuf1,1,xbuf3,1)
        else if (nincore.eq.2) then
          call vec_from_da(ff_sbsp2,irec,xbuf2,nwfpar)
          xmat(ndim1,jj) = ddot(nwfpar,xbuf1,1,xbuf2,1)
        !else if (nincore.eq.1) then

          ! routine needed that makes inner product betw. vector in core
          ! and vector on disc
        else 
          xmat(ndim1,jj) =
     &         da_ddot(ff_sbsp1,jrec1,ff_sbsp2,irec,
     &         nwfpar,xbuf1,xbuf2,lenbuf)
        end if

      end do

      ! c) update last row
c      if (nincore.gt.1.and.nincore.lt.3) the
      if (nincore.gt.2.and.nincore.lt.3) then
        ! reload w(n)
        call vec_from_da(ff_sbsp2,jrec2,xbuf2,nwfpar)
      end if

      do irec = 1, ndim1

        jj = iord1(irec)

        if (jj.eq.ndim1) cycle

        if (nincore.eq.3) then
          call vec_from_da(ff_sbsp1,irec,xbuf3,nwfpar)
          xmat(jj,ndim1) = ddot(nwfpar,xbuf2,1,xbuf3,1)
        else if (nincore.eq.2) then
          call vec_from_da(ff_sbsp1,irec,xbuf1,nwfpar)
          xmat(jj,ndim1) = ddot(nwfpar,xbuf1,1,xbuf2,1)
        !else if (nincore.eq.1) then

          ! routine needed that makes inner product betw. vector in core
          ! and vector on disc
        else 
          xmat(jj,ndim1) =
     &         da_ddot(ff_sbsp2,jrec2,ff_sbsp1,irec,
     &         nwfpar,xbuf1,xbuf2,lenbuf)
        end if

      end do

      if (ntest.ge.20) then
        write(luout,*) 'updated subspace matrix:'
        call wrtmat2(xmat,ndim1,ndim1,mxdim,mxdim)
      end if

      return

      end
