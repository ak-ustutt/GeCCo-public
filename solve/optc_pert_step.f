
*----------------------------------------------------------------------*
      subroutine optc_pert_step(ffvec,ffgrd,ffdia,trrad,
     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3,ffscr)
*----------------------------------------------------------------------*
*
*     make |vec> = |vec> - (D^-1+damp)|grd>
*
*     where damp is chosen to have <grd|(D^-1+damp)^2|grd> <= trrad
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00

      type(filinf), intent(in) ::
     &     ffvec, ffgrd, ffdia, ffscr
      integer, intent(in) ::
     &     nincore, nwfpar, lenbuf
      real(8), intent(in) ::
     &     trrad
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)

      logical ::
     &     converged
      integer ::
     &     ixd_iter
      real(8) ::
     &     xdamp, xdamp_min, xnrm, xval

      real(8), external ::
     &     ddot, da_ddot, fndmnx, da_fndmnx

      if (ntest.gt.0) then
        write(luout,*) '============================'
        write(luout,*) ' hi, this is optc_pert_step'
        write(luout,*) '============================'
      end if

      if (nincore.ge.2) then
        call vec_from_da(ffdia,1,xbuf2,nwfpar)
        xdamp_min = max(0d0,-fndmnx(xbuf2,nwfpar,-1))
      else
        xdamp_min = max(0d0,-da_fndmnx(ffdia,1,-1,
     &       nwfpar,xbuf1,lenbuf))
      end if

      xdamp = 0d0
      if (xdamp_min.gt.0d0) xdamp = xdamp_min + 0.1d0
c     if (nincore.ge.1)
      if (nincore.ge.2)
     &     call vec_from_da(ffgrd,1,xbuf1,nwfpar)
      
      converged = .false.
      ixd_iter = 0

      do while (.not.converged)

        if (nincore.eq.3) then
        
          xbuf3(1:nwfpar) = xbuf1(1:nwfpar)/(xbuf2(1:nwfpar) + xdamp)
          xnrm = sqrt(ddot(nwfpar,xbuf3,1,xbuf3,1))

        else if (nincore.eq.2) then

          xbuf2(1:nwfpar) = xbuf1(1:nwfpar)/(xbuf2(1:nwfpar) + xdamp)
          xnrm = sqrt(ddot(nwfpar,xbuf2,1,xbuf2,1))

c      else if (nincore.eq.1) then
        ! vector in-core = vector incore / vector o/o/core 

        else

          call da_diavec(ffscr,1,0d0,
     &         ffgrd,1,1d0,
     &         ffdia,1,xdamp,-1d0,
     &         nwfpar,xbuf1,xbuf2,lenbuf)
          xnrm = sqrt(da_ddot(ffscr,1,ffscr,1,
     &         nwfpar,xbuf1,xbuf1,lenbuf))
          
        end if

        xval = xnrm - trrad
        call optc_xdamp_ctl(xdamp,xdamp_min,xval,converged,ixd_iter)

      end do

      if (nincore.ge.2) then
        call vec_from_da(ffvec,1,xbuf1,nwfpar)

        if (nincore.eq.2) then
          xbuf1(1:nwfpar) = xbuf1(1:nwfpar) - xbuf2(1:nwfpar)
        else
          xbuf1(1:nwfpar) = xbuf1(1:nwfpar) - xbuf3(1:nwfpar)
        end if

        call vec_to_da(ffvec,1,xbuf1,nwfpar)
        
c      else if (nincore.eq.1) then
        ! add incore vector to o/o/core vector

      else
        call da_vecsum(ffvec,1,ffvec,1,1d0,ffscr,1,-1d0,
     &       nwfpar,xbuf1,xbuf2,lenbuf)

      end if

      return
      end
