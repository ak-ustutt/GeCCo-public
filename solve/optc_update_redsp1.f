*----------------------------------------------------------------------*
      subroutine optc_update_redsp1(xmat,ndim,mxdim,shift,init,
     &     iord,ff_sbsp,
     &     nincore,nwfpar,
     &     lenbuf,xbuf1,xbuf2,xbuf3)
*----------------------------------------------------------------------*
*
*     update reduced-space matrix xmat(i,j) = <v(i)|v(j)>, where
*     xmat is triangular packed
*     
*     for incore-version: most recent vector expected on xbuf1
*                         (remains unchanged)
*
*     shift==.true.  xmat(i,j) := xmat(i+1,j+1), xmat(1,j) is discarded
*     init==.true.   init new row/col to 0d0
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00
 
      type(filinf), intent(in) ::
     &     ff_sbsp
      logical, intent(in) ::
     &     shift, init
      integer, intent(in) ::
     &     ndim, mxdim, nincore, nwfpar, lenbuf,
     &     iord(*)
      real(8), intent(inout), target ::
     &     xmat(*), xbuf1(*), xbuf2(*), xbuf3(*)
      
      integer ::
     &     iioff, iioff2, ii, jj, irec, jrec, idx, iblk, len

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     ddot, da_ddot

      if (ntest.gt.0) then
        write(luout,*) '==============================='
        write(luout,*) ' welcome to optc_update_redsp1'
        write(luout,*) '==============================='
        write(luout,*) 'nwfpar = ',nwfpar
      end if

      if (ntest.ge.20) then
        write(luout,*) 'subspace matrix on input:'
        call prtrlt(xmat,ndim)
      end if

      ! subspace is full: delete contribution from first vector
      ! i.e. shift all matrix-elements down by one
      if (shift) then
        do ii = 1, ndim-1
          iioff  = ii*(ii-1)/2
          iioff2 = (ii+1)*ii/2 + 1
          do jj = 1, ii
            xmat(iioff+jj) = xmat(iioff2+jj)
          end do
        end do
        if (ntest.ge.20) then
          write(luout,*) 'subspace matrix after shift:'
          call prtrlt(xmat,ndim)
        end if
      end if

      ! update last column
      iioff = ndim*(ndim-1)/2

      ! init to zero, if requested:
      if (init) then
        do ii = 1, ndim
          xmat(iioff+ii) = 0d0
        end do
      end if

      if (ndim.eq.0)
     &     call quit(1,'optc_update_redsp1','ndim.eq.0')
      jrec = ioptc_get_sbsp_rec(ndim,iord,ndim,mxdim)
      do irec = 1, ndim

        ii = iord(irec)

        if (nincore.ge.2) then
          if (irec.eq.jrec) then
            xmat(iioff+ii) = xmat(iioff+ii)
     &           +ddot(nwfpar,xbuf1,1,xbuf1,1)
          else
            call vec_from_da(ff_sbsp,irec,xbuf2,nwfpar)
            xmat(iioff+ii) = xmat(iioff+ii)
     &           +ddot(nwfpar,xbuf1,1,xbuf2,1)
          end if

        !else if (nincore.eq.1) then

          ! routine needed that makes inner product betw. vector in core
          ! and vector on disc

        else
          xmat(iioff+ii) = xmat(iioff+ii) +
     &         da_ddot(ff_sbsp,jrec,1,ff_sbsp,irec,1,
     &         nwfpar,xbuf1,xbuf2,lenbuf)

        end if

      end do

      if (ntest.ge.20) then
        write(luout,*) 'updated subspace matrix:'
        call prtrlt(xmat,ndim)
      end if

      return

      end
