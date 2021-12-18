*----------------------------------------------------------------------*
      subroutine optc_xdamp_ctl(xdamp,xmindamp,xval,converged,iter)
*----------------------------------------------------------------------*
*
*     control search after optimal xdamp when no analytical derivative
*     is available
*
*     not re-entrant !
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 000, mxd_iter = 100
      real(8), parameter ::
     &     thrdamp = 1d-4, xinc = 1d-5, xfailfac = 1d5,
     &     xmxdxdamp = 0.1d0

      logical, intent(out) ::
     &     converged
      integer, intent(inout) ::
     &     iter
      real(8), intent(inout) ::
     &     xdamp
      real(8), intent(in) ::
     &     xval, xmindamp
      
      logical, save ::
     &     inuse = .false.

      real(8), save ::
     &     xlow, xhigh, xvlow, xvhigh, xdamp_min,
     &     xd(mxd_iter), xv(mxd_iter)
      integer, save ::
     &     isubcnt

      integer ::
     &     ii
      real(8) ::
     &     xdiff, xgr, dxdamp, dxdsign

      if (inuse.and.iter.eq.0) then
        write(lulog,*) 'fatal: initialization call to optc_xdamp_ctl '//
     &       'while used by other instance!'
        call quit(1,'optc_xdamp_ctl','fatal -- re-entrant call')
      end if

      if (iter.eq.0) then
        xdamp_min = xmindamp
        xlow = xdamp_min
        xhigh = 1d13
        xvlow = 100d0
        xvhigh = -100d0
        isubcnt = 0
      end if
      
      converged = (abs(xdamp).lt.1d-6.and.xval.le.0d0).or.
     &     (abs(xval).lt.max(abs(xdamp)*1d-6,1d-6))

c        if (ixd_iter.gt.0)
      write(lulog,'(x,">>",2i4,2e12.3,l)')
     &       iter,isubcnt,xdamp,xval,converged

      if (converged) then
        inuse = .false.
        return
      end if

      if (xval.gt.0d0) then
        xlow = xdamp
        xvlow = xval
        dxdsign = +1d0          ! expect step in this direction
      else
        xhigh = xdamp
        xvhigh = xval
        dxdsign = -1d0
      end if

      if (xdamp.lt.0d0) then    ! that should never ever happen!!
        write(lulog,*) 'What did you do? xdamp = ',xdamp
        call quit(1,'optc_xdamp_ctl','xdamn')
      end if

      iter = iter+1
      if (iter.ge.mxd_iter) then
        call quit(0,'optc_xdamp_ctl','unsolvable problems in damping')
      end if
      xd(iter) = xdamp
      xv(iter) = xval

      if (isubcnt.eq.0) then
        ! restart at a very new point
        ! make a small increment to get the gradient              
        if (xv(iter).gt.0d0) then
          dxdamp = + max(xinc,abs(xdamp*xinc))
        else
          dxdamp = - max(xinc,abs(xdamp*xinc))
        end if
        if (ntest.ge.150)
     &       write(lulog,*)
     &       'xdamp> initialized num. gradient with xinc = ',
     &       dxdamp
             ! next time we see us in step two
        isubcnt = 2
      else if (isubcnt.eq.1) then
        ! this is step one of a numerical Newton optimization
        ! make a small increment to get the gradient              
        if (xv(iter).gt.0d0) then
          dxdamp = + max(xinc,abs(xdamp*xinc))
        else
          dxdamp = - max(xinc,abs(xdamp*xinc))
        end if
        if (ntest.ge.150)
     &       write(lulog,*)
     &       'xdamp> initialized num. gradient with xinc = ',
     &       dxdamp
        ! next time we see us in step two
        isubcnt = 2
      else if (isubcnt.eq.2) then
        ! get gradient from finite difference
        if (ntest.ge.150) then
          write(lulog,*) 'xdamp> linear model '
          write(lulog,*) ' points used:'
          do ii = 0, 1
            write(lulog,'(3x,i2,2(2x,e25.8))')
     &           ii, xd(iter-ii), xv(iter-ii)
          end do
        end if
        xdiff = xd(iter)-xd(iter-1)
        xgr = (xv(iter)-xv(iter-1))/xdiff
        if (ntest.ge.150)
     &             write(lulog,*) 'xdamp> current num. gradient: ',xgr    
        ! gradient has to be negative
        if (xgr.lt.0d0) then
          if (ntest.ge.150)
     &             write(lulog,*) 'xdamp> accepting step'  
          ! make a Newton step
          dxdamp = - xv(iter)/xgr
          isubcnt = 1
          if (ntest.eq.150)
     &         write(lulog,*) 'xdamp> step = ',dxdamp
          if (dxdamp.lt.0.5d0) isubcnt = 2
        else
          if (ntest.ge.150)
     &             write(lulog,*) 'xdamp> search at a new place'  
          ! retry with a larger xdamp
          dxdamp = +10d0       !+ xfailfac*xinc
          iter = iter-1
          isubcnt = 0
        end if
      else
        write(lulog,*) 'unknown isubcnt = ',isubcnt
        call quit(1,'optc_xdamp_ctl','unexpected event')
      end if

      ! unless isubcnt was reset to 0
      if (isubcnt.ne.0) then
        if (xdamp+dxdamp.ge.xhigh) then
          if (xval.eq.xvhigh)
     &         call quit(1,'optc_xdamp_ctl','unexpected event')
          xdamp = xdamp-xval*(xdamp-xhigh)/(xval-xvhigh)
        else if (xdamp+dxdamp.le.xlow) then
          if (xlow.eq.xdamp_min) xvlow=-xval
          if (xval.eq.xvlow)
     &         call quit(1,'optc_xdamp_ctl','unexpected event')
          xdamp = xdamp-xval*(xdamp-xlow)/(xval-xvlow)
        else
          xdamp = xdamp+dxdamp
        end if
      else
        xdamp = xdamp+dxdamp
      end if
      
      return
      end
