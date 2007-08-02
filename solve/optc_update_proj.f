*----------------------------------------------------------------------*
      subroutine optc_update_proj
     &     (xvec,nvec,ndim,mxdim,shift,
     &     iord,ff_sbsp,ffvec,
     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2)
*----------------------------------------------------------------------*
*
*     update reduced-space vectors xvec(i,j) = <v(i)|x_j>
*
*     i = 1, ndim;  j = 1, nvec; mxdim: leading dim for xvec
*     
*     for incore-version: most recent vector v expected on xbuf1
*
*     shift==.true.  xvec(i) := xvec(i+1), xvec(1) is discarded
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00

      type(filinf) ::
     &     ff_sbsp, ffvec
      logical, intent(in) ::
     &     shift
      integer, intent(in) ::
     &     ndim, mxdim, nincore, nwfpar, lenbuf,
     &     iord(*), nvec
      real(8), intent(inout) ::
     &     xvec(mxdim,nvec), xbuf1(*), xbuf2(*)
      
      integer ::
     &     iioff, iioff2, ii, jj, jrec

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     ddot, da_ddot

      if (ntest.gt.0) then
        call write_title(luout,wst_dbg_subr,
     &       'welcome to optc_update_proj')
      end if

      if (ntest.ge.20) then
        write(luout,*) 'subspace vectors on input:'
        call wrtmat2(xvec,ndim,nvec,mxdim,nvec)
      end if

      ! subspace is full: delete contribution from first vector
      ! i.e. shift all elements down by one
      if (shift) then
        do jj = 1, nvec
          do ii = 1, ndim - 1
            xvec(ii,jj) = xvec(ii+1,jj)
          end do
        end do
        if (ntest.ge.20) then
          write(luout,*) 'subspace vector after shift:'
          call wrtmat2(xvec,ndim,nvec,mxdim,nvec)
        end if
      end if

      ! add new projection (= last element in subspace)
      jrec = ioptc_get_sbsp_rec(ndim,iord,ndim,mxdim)
      if (nincore.ge.2) then
        do jj = 1, nvec
          call vec_from_da(ffvec,jj,xbuf2,nwfpar)
          xvec(ndim,jj) = ddot(nwfpar,xbuf1,1,xbuf2,1)
        end do
      else
        do jj = 1, nvec
          xvec(ndim,jj) = da_ddot(ff_sbsp,jrec,ffvec,jj,
     &         nwfpar,xbuf1,xbuf2,lenbuf)
        end do
      end if

      if (ntest.ge.20) then
        write(luout,*) 'updated subspace vector:'
        call wrtmat2(xvec,ndim,nvec,mxdim,nvec)
      end if

      return

      end
