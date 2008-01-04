*----------------------------------------------------------------------*
      subroutine optc_expand_vec(xvec,ndim,xnrm,getnrm,
     &     ffamp,irecamp,xfac,ff_sbsp,iord,
     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2)
*----------------------------------------------------------------------*
*     expand subspace vector to full basis
*
*     incore: xbuf1 contains ffamp on entry and on exit
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00

      type(filinf), intent(in) ::
     &     ffamp, ff_sbsp
      integer, intent(in) ::
     &     ndim, nincore, nwfpar, lenbuf,
     &     iord(*), irecamp
      real(8), intent(in) ::
     &     xvec(*), xfac
      logical, intent(in) ::
     &     getnrm
      real(8), intent(out) ::
     &     xnrm
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*)

      integer ::
     &     irec, ii
      real(8) ::
     &     xveco(ndim), xscr

      real(8), external ::
     &     ddot, da_ddot

      if (ntest.gt.0) then
        write(luout,*) '==========================='
        write(luout,*) ' info from optc_expand_vec' 
        write(luout,*) '==========================='
        write(luout,*) ' ndim : ', ndim
        write(luout,*) ' ndim : ', nwfpar
        write(luout,*) ' xvec : ', xvec(1:ndim)
        write(luout,*) ' xfac : ', xfac
      end if

      if (nincore.ge.2) then

        if (xfac.eq.0d0) then
          xbuf1(1:nwfpar) = 0d0
        else if (xfac.ne.1d0) then
          xbuf1(1:nwfpar) = xfac*xbuf1(1:nwfpar)
        end if

        if (ntest.gt.30) then
          xscr = sqrt(ddot(nwfpar,xbuf1,1,xbuf1,1))
          write(luout,*) ' |initial| = ',xscr
        end if

        do irec = 1, ndim
          ii = iord(irec)
          if (ntest.ge.20) then
            write(luout,*) 'record: ', irec,ii,xvec(ii)
          end if
          if (xvec(ii).eq.0d0) cycle
          call vec_from_da(ff_sbsp,irec,xbuf2,nwfpar)
          if (ntest.ge.30) then
            xscr = sqrt(ddot(nwfpar,xbuf2,1,xbuf2,1))
            write(luout,*) ' |v| = ',xscr
          end if
          xbuf1(1:nwfpar) = xbuf1(1:nwfpar)+xvec(ii)*xbuf2(1:nwfpar)
        end do

        if (getnrm.or.ntest.gt.30)
     &       xscr = sqrt(ddot(nwfpar,xbuf1,1,xbuf1,1))
        if (getnrm) xnrm = xscr

        if (ntest.gt.30) then
          write(luout,*) ' |final| = ',xscr
        end if
        call vec_to_da(ffamp,irecamp,xbuf1,nwfpar)

!      else if (nincore.eq.1) then

        ! needed: add vector on disc to vector in core

      else
      
        do irec = 1, ndim
          xveco(irec) = xvec(iord(irec))
        end do

        call da_matvec(ffamp,irecamp,xfac,
     &       ff_sbsp,1,xveco,ndim,0d0,
     &       nwfpar,xbuf1,xbuf2,lenbuf)
        
        if (getnrm)
     &       xnrm = da_ddot(ffamp,irecamp,ffamp,irecamp,
     &       nwfpar,xbuf1,xbuf2,lenbuf)

      end if

      return

      end
