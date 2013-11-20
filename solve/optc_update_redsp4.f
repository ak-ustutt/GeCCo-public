*----------------------------------------------------------------------*
      subroutine optc_update_redsp4(
     &     xmat,smat,xvec,ndim,nrhs,mxdim,nadd,ndel,init,
     &     iordv,ff_vsbsp,iordw,ff_wsbsp,iordx,ff_xsbsp,
     &     ff_rhs,
     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
*----------------------------------------------------------------------*
*
*     update reduced-space matrices xmat(i,j) = <v(i)|w(j)>,
*                               and smat(i,j) = <v(i)|x(j)>, 
*     and the reduced-sp. vector  xvec(i,k) = <v(i)|rhs(k)> (if nrhs>0)
*
*     ndim is the current subspace dimenstion
*     nadd is the number of entries that need be updated (always from
*           (ndim-nadd+1) to (ndim)
*     ndel is the number of entries to be deleted first (matrix is
*          shifted accordingly, always entries 1 to ndel are deleted)
*     init if (.true.) init new entries to 0d0 and perform shift of 
*          matrix if requested by ndel
*
*     xmat and xvec have the leading dimension mxdim
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00

      logical ::
     &     init
      type(filinf) ::
     &     ff_vsbsp, ff_wsbsp, ff_xsbsp, ff_rhs
      integer, intent(in) ::
     &     ndim, nrhs, nadd, ndel, mxdim, nincore, nwfpar, lenbuf,
     &     iordv(*), iordw(*), iordx(*)
      real(8), intent(inout) ::
     &     xmat(mxdim,*), smat(mxdim,*),
     &     xvec(mxdim,*), xbuf1(*), xbuf2(*), xbuf3(*)
      
      integer ::
     &     ii, jj, irec, jrec, jrec_last, rhsrec, rhsrec_last

      real(8), external ::
     &     ddot, da_ddot

      if (ntest.gt.0) then
        write(lulog,*) '==============================='
        write(lulog,*) ' welcome to optc_update_redsp4'
        write(lulog,*) '==============================='
        write(lulog,*) ' ndim, nadd, ndel: ',ndim, nadd, ndel
        write(lulog,*) ' init = ',init
      end if

      if (ntest.ge.20) then
        write(lulog,*) 'subspace matrix 1 on input:'
        call wrtmat2(xmat,ndim,ndim,mxdim,mxdim)
        write(lulog,*) 'subspace matrix 2 on input:'
        call wrtmat2(smat,ndim,ndim,mxdim,mxdim)
        if (nrhs.gt.0) then
          write(lulog,*) 'subspace rhs on input:'
          call wrtmat2(xvec,ndim,nrhs,mxdim,nrhs)
        end if
      end if

      ! subspace is full: delete contribution from first ndel vectors
      ! i.e. shift all matrix-elements down by one
      if (ndel.gt.0.and.init) then
        do jj = 1, ndim - ndel
          do ii = 1, ndim - ndel
            xmat(ii,jj) = xmat(ii+ndel,jj+ndel)
          end do
        end do
        do jj = 1, ndim - ndel
          do ii = 1, ndim - ndel
            smat(ii,jj) = smat(ii+ndel,jj+ndel)
          end do
        end do
        if (ntest.ge.20) then
          write(lulog,*) 'subspace matrix 1 after shift:'
          call wrtmat2(xmat,ndim,ndim,mxdim,mxdim)
          write(lulog,*) 'subspace matrix 2 after shift:'
          call wrtmat2(smat,ndim,ndim,mxdim,mxdim)
        end if
      end if

      if (nrhs.gt.0.and.ndel.gt.0.and.init) then
        do jj = 1, nrhs
          do ii = 1, ndim - ndel
            xmat(ii,jj) = xmat(ii+ndel,jj)
          end do
        end do
        if (ntest.ge.20) then
          write(lulog,*) 'subspace rhs after shift:'
          call wrtmat2(xmat,ndim,nrhs,mxdim,nrhs)
        end if
      end if

      if (init) then
        if (nrhs.gt.0) then
          xvec(ndim-nadd+1:ndim,1:nrhs) = 0d0
        end if
        do ii = 1, ndim
          xmat(ii,ndim-nadd+1:ndim) = 0d0
          xmat(ndim-nadd+1:ndim,ii) = 0d0
        end do
        do ii = 1, ndim
          smat(ii,ndim-nadd+1:ndim) = 0d0
          smat(ndim-nadd+1:ndim,ii) = 0d0
        end do
      end if

      ! loop over records of <v|-file
      jrec_last = -1
      rhsrec_last = -1
      do irec = 1, ndim

        ii = iordv(irec)
        
        ! if we can hold the vector in core, load it
        if (nincore.ge.2)
     &       call vec_from_da(ff_vsbsp,irec,xbuf1,nwfpar)
        
        ! loop over records of |w>-file
        do jrec = 1, ndim
          
          jj = iordw(jrec)

          ! is this a block of xmat to be updated?
          if (ii.gt.ndim-nadd.or.jj.gt.ndim-nadd) then

            ! incore version:
            if (nincore.ge.2) then
c              ! avoid loading if record is in core
c              if (jrec_last.ne.jrec) then
                jrec_last = jrec
                call vec_from_da(ff_wsbsp,jrec,xbuf2,nwfpar)
c              end if
c dbg
              if (ii.eq.jj) print *,'cartesian scalar product: ',
     &          ddot(nwfpar,xbuf1,1,xbuf1,1)
c dbgend
              xmat(ii,jj) = xmat(ii,jj)
     &                    + ddot(nwfpar,xbuf1,1,xbuf2,1)
            else
              xmat(ii,jj) = xmat(ii,jj)
     &                    + da_ddot(ff_vsbsp,irec,ff_wsbsp,jrec,
     &                              nwfpar,xbuf1,xbuf2,lenbuf)
            end if
          end if

        end do

        ! loop over records of |x>-file
        do jrec = 1, ndim
          
          jj = iordx(jrec)

          ! is this a block of xmat to be updated?
          if (ii.gt.ndim-nadd.or.jj.gt.ndim-nadd) then

            ! incore version:
            if (nincore.ge.2) then
c              ! avoid loading if record is in core
c              if (jrec_last.ne.jrec) then
                jrec_last = jrec
                call vec_from_da(ff_xsbsp,jrec,xbuf2,nwfpar)
c              end if
              smat(ii,jj) = smat(ii,jj)
     &              + ddot(nwfpar,xbuf1,1,xbuf2,1)
            else
              smat(ii,jj) = smat(ii,jj)
     &            + da_ddot(ff_vsbsp,irec,ff_xsbsp,jrec,
     &                      nwfpar,xbuf1,xbuf2,lenbuf)
            end if
          end if

        end do

        if (ii.gt.ndim-nadd) then
          do rhsrec = 1, nrhs

            ! we can use xbuf3 then ...
            if (nincore.eq.3) then
              ! avoid loading if record is in core
              if (rhsrec_last.ne.rhsrec) then
                rhsrec_last = rhsrec
                call vec_from_da(ff_rhs,rhsrec,xbuf3,nwfpar)                
              end if
                xvec(ii,rhsrec) = xvec(ii,rhsrec)
     &                  + ddot(nwfpar,xbuf1,1,xbuf3,1)
            else if (nincore.eq.2) then
              jrec_last = -1 ! signal that xbuf2 is destroyed
              call vec_from_da(ff_rhs,rhsrec,xbuf2,nwfpar)
              xvec(ii,rhsrec) = xvec(ii,rhsrec)
     &                + ddot(nwfpar,xbuf1,1,xbuf2,1)
            else
              xvec(ii,rhsrec) = xvec(ii,rhsrec)
     &            + da_ddot(ff_vsbsp,irec,ff_rhs,rhsrec,
     &                      nwfpar,xbuf1,xbuf2,lenbuf)
            end if
          end do
        end if

      end do

      if (ntest.ge.10) then
        write(lulog,*) 'updated subspace matrix 1:'
        call wrtmat2(xmat,ndim,ndim,mxdim,mxdim)
        write(lulog,*) 'updated subspace matrix 2:'
        call wrtmat2(smat,ndim,ndim,mxdim,mxdim)
        if (nrhs.gt.0) then
          write(lulog,*) 'updated subspace rhs:'
          call wrtmat2(xvec,ndim,nrhs,mxdim,nrhs)
        end if
      end if

      return

      end
