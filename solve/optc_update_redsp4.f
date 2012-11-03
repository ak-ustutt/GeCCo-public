*----------------------------------------------------------------------*
      subroutine optc_update_redsp4(
     &     xmat,smat,xvec,ndim,nrhs,mxdim,nadd,ndel,init,
     &     iordv,ff_vsbsp,iordw,ff_wsbsp,iordx,ff_xsbsp,
     &     ff_rhs,
     &     nsec,nwfpsec,idstsec,signsec,
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
     &     iordv(*), iordw(*), iordx(*), nsec, nwfpsec(*), idstsec(*)
      real(8), intent(in) ::
     &     signsec(*)
      real(8), intent(inout) ::
     &     xmat(mxdim,*), smat(mxdim,*),
     &     xvec(mxdim,*), xbuf1(*), xbuf2(*), xbuf3(*)
      
      integer ::
     &     ii, jj, irec, jrec, jrec_last, rhsrec, rhsrec_last, isec

      real(8), external ::
     &     ddot, da_ddot

      if (ntest.gt.0) then
        write(luout,*) '==============================='
        write(luout,*) ' welcome to optc_update_redsp4'
        write(luout,*) '==============================='
        write(luout,*) ' ndim, nadd, ndel: ',ndim, nadd, ndel
        write(luout,*) ' init = ',init
      end if

      if (ntest.ge.20) then
        write(luout,*) 'subspace matrix 1 on input:'
        call wrtmat2(xmat,ndim,ndim,mxdim,mxdim)
        write(luout,*) 'subspace matrix 2 on input:'
        call wrtmat2(smat,ndim,ndim,mxdim,mxdim)
        if (nrhs.gt.0) then
          write(luout,*) 'subspace rhs on input:'
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
          write(luout,*) 'subspace matrix 1 after shift:'
          call wrtmat2(xmat,ndim,ndim,mxdim,mxdim)
          write(luout,*) 'subspace matrix 2 after shift:'
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
          write(luout,*) 'subspace rhs after shift:'
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
              do isec = 1, nsec
                xmat(ii,jj) = xmat(ii,jj) + signsec(isec)
     &               * ddot(nwfpsec(isec),xbuf1(idstsec(isec)),1,
     &                      xbuf2(idstsec(isec)),1)
              end do
            else
              do isec = 1, nsec
c dbg  see evpc_core for explanation
                  if (idstsec(isec).ne.1)
     &             call quit(1,'**08**','bug in the code!')
c dbgend
                xmat(ii,jj) = xmat(ii,jj) + signsec(isec)
     &              * da_ddot(ff_vsbsp,irec,idstsec(isec),ff_wsbsp,jrec,
     &                   idstsec(isec),nwfpsec(isec),xbuf1,xbuf2,lenbuf)
              end do
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
              do isec = 1, nsec
                smat(ii,jj) = smat(ii,jj) + signsec(isec)
     &                * ddot(nwfpsec(isec),xbuf1(idstsec(isec)),1,
     &                       xbuf2(idstsec(isec)),1)
              end do
            else
              do isec = 1, nsec
c dbg  see evpc_core for explanation
                  if (idstsec(isec).ne.1)
     &             call quit(1,'**09**','bug in the code!')
c dbgend
                smat(ii,jj) = smat(ii,jj) + signsec(isec)
     &              * da_ddot(ff_vsbsp,irec,idstsec(isec),ff_xsbsp,jrec,
     &                   idstsec(isec),nwfpsec(isec),xbuf1,xbuf2,lenbuf)
              end do
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
              do isec = 1, nsec
                xvec(ii,rhsrec) = xvec(ii,rhsrec) + signsec(isec)
     &                  * ddot(nwfpsec(isec),xbuf1(idstsec(isec)),1,
     &                         xbuf3(idstsec(isec)),1)
              end do
            else if (nincore.eq.2) then
              jrec_last = -1 ! signal that xbuf2 is destroyed
              call vec_from_da(ff_rhs,rhsrec,xbuf2,nwfpar)
              do isec = 1, nsec
                xvec(ii,rhsrec) = xvec(ii,rhsrec) + signsec(isec)
     &                  * ddot(nwfpsec(isec),xbuf1(idstsec(isec)),1,
     &                         xbuf2(idstsec(isec)),1)
              end do
            else
              do isec = 1, nsec
c dbg  see evpc_core for explanation
                  if (idstsec(isec).ne.1)
     &             call quit(1,'**10**','bug in the code!')
c dbgend
                xvec(ii,rhsrec) = xvec(ii,rhsrec) + signsec(isec)
     &              * da_ddot(ff_vsbsp,irec,idstsec(isec),ff_rhs,rhsrec,
     &                   idstsec(isec),nwfpsec(isec),xbuf1,xbuf2,lenbuf)
              end do
            end if
          end do
        end if

      end do

      if (ntest.ge.10) then
        write(luout,*) 'updated subspace matrix 1:'
        call wrtmat2(xmat,ndim,ndim,mxdim,mxdim)
        write(luout,*) 'updated subspace matrix 2:'
        call wrtmat2(smat,ndim,ndim,mxdim,mxdim)
        if (nrhs.gt.0) then
          write(luout,*) 'updated subspace rhs:'
          call wrtmat2(xvec,ndim,nrhs,mxdim,nrhs)
        end if
      end if

      return

      end
