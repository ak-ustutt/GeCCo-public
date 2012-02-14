*----------------------------------------------------------------------*
      subroutine spinsym_traf(mode,ndim,mat,map,nsing,sing,trip,half)
*----------------------------------------------------------------------*
*     does a unitary (half-)transformation using a spin-symmetrization
*     map as input for the unitary matrix U.
*     the result is split into a "singlet" and a "triplet" block;
*     in case of the back-transformation (mode=2), these two blocks
*     are rejoined.
*
*     matthias, feb 2010
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00
      integer, intent(in) ::
     &     mode, ndim, nsing, map(ndim)
      real(8), intent(inout), target ::
     &     mat(ndim,ndim), sing(nsing,nsing),
     &     trip(ndim-nsing,ndim-nsing)
      logical, intent(in) ::
     &     half
      real(8) ::
     &     sig, fac
      real(8), pointer ::
     &     scr(:,:)

      integer ::
     &     icol, iline, ntrip, ising, itrip
      fac = 1d0/sqrt(2d0)

      ntrip = ndim - nsing
      if (mode.lt.1.or.mode.gt.2.or.ntrip.lt.0
     &    .or.(half.and.mode.eq.1))
     &     call quit(1,'spinsym_traf','are you kidding me?')

      if (ntest.ge.100) then
        write(luout,'(x,a)') '-------------------------------'
        write(luout,'(x,a,i1,a)') 'spinsym_traf at work (mode = ',
     &                            mode,')'
        write(luout,'(x,a)') '-------------------------------'
      end if

      if (mode.eq.1) then
        if (ntest.ge.100) then
          write(luout,*) 'input matrix:'
          call wrtmat2(mat,ndim,ndim,ndim,ndim)
        end if
        allocate(scr(ndim,ndim))
        scr = 0d0
        sing = 0d0
        trip = 0d0

        ! first mat*U
        ising = 0
        itrip = 0
        do icol = 1, ndim
          if (abs(map(icol)).gt.ndim) call quit(1,'spinsym_traf',
     &           'map must stay within matrix dimension!')
          if (abs(map(icol)).eq.icol) then
            if (map(icol).gt.0) then
              ising = ising + 1
              scr(1:ndim,ising) = mat(1:ndim,icol)
            else
              itrip = itrip + 1
              scr(1:ndim,nsing+itrip) = mat(1:ndim,icol)
            end if
          else if (abs(map(icol)).gt.icol) then
            sig = dble(sign(1,map(icol)))
            ising = ising + 1
            scr(1:ndim,ising)
     &       = fac*mat(1:ndim,icol)
     &       + sig*fac*mat(1:ndim,abs(map(icol)))
          else
            sig = - dble(sign(1,map(icol)))
            itrip = itrip + 1
            scr(1:ndim,nsing+itrip)
     &       = fac*mat(1:ndim,abs(map(icol)))
     &       + sig*fac*mat(1:ndim,icol)
          end if
        end do
        if (ising.ne.nsing.or.ising+itrip.ne.ndim) call quit(1,
     &      'spinsym_traf','error in singlet triplet splitting')

        ! now U^+*(mat*U) --> (sing|trip)
        ising = 0
        itrip = 0
        do iline = 1, ndim
          if (abs(map(iline)).eq.iline) then
            if (map(iline).gt.0) then
              ising = ising + 1
              sing(ising,1:nsing) = scr(iline,1:nsing)
            else
              itrip = itrip + 1
              trip(itrip,1:ntrip) = scr(iline,nsing+1:ndim)
            end if
          else if (abs(map(iline)).gt.iline) then
            sig = dble(sign(1,map(iline)))
            ising = ising + 1
            sing(ising,1:nsing)
     &       = fac*scr(iline,1:nsing)
     &       + sig*fac*scr(abs(map(iline)),1:nsing)
          else
            sig = - dble(sign(1,map(iline)))
            itrip = itrip + 1
            trip(itrip,1:ntrip)
     &       = fac*scr(abs(map(iline)),nsing+1:ndim)
     &       + sig*fac*scr(iline,nsing+1:ndim)
          end if
        end do
        deallocate(scr)
        if (ntest.ge.100) then
          write(luout,*) 'output "singlet" block:'
          call wrtmat2(sing,nsing,nsing,nsing,nsing)
          write(luout,*) 'output "triplet" block:'
          call wrtmat2(trip,ntrip,ntrip,ntrip,ntrip)
        end if

      else if (mode.eq.2) then
        if (ntest.ge.100) then
          write(luout,*) 'input "singlet" block:'
          call wrtmat2(sing,nsing,nsing,nsing,nsing)
          write(luout,*) 'input "triplet" block:'
          call wrtmat2(trip,ntrip,ntrip,ntrip,ntrip)
        end if

        if (half) then
          scr => mat
        else
          allocate(scr(ndim,ndim))
        end if

        ! first U*(sing|trip)
        scr = 0d0
        ising = 0
        itrip = 0
        do iline = 1, ndim
          if (abs(map(iline)).gt.ndim) call quit(1,'spinsym_traf',
     &           'map must stay within matrix dimension!')
          if (abs(map(iline)).eq.iline) then
            if (map(iline).gt.0) then
              ising = ising + 1
              scr(iline,1:nsing) = sing(ising,1:nsing)
            else
              itrip = itrip + 1
              scr(iline,nsing+1:ndim) = trip(itrip,1:ntrip)
            end if
          else if (abs(map(iline)).gt.iline) then
            sig = dble(sign(1,map(iline)))
            ising = ising + 1
            scr(iline,1:nsing) = fac*sing(ising,1:nsing)
            scr(abs(map(iline)),1:nsing) = sig*fac*sing(ising,1:nsing)
          else
            sig = - dble(sign(1,map(iline)))
            itrip = itrip + 1
            scr(abs(map(iline)),nsing+1:ndim) = fac*trip(itrip,1:ntrip)
            scr(iline,nsing+1:ndim) = sig*fac*trip(itrip,1:ntrip)
          end if
        end do
        if (ising.ne.nsing.or.ising+itrip.ne.ndim) call quit(1,
     &      'spinsym_traf','error in singlet triplet splitting')

        ! if full trafo is required, do [U*(sing|trip)]*U^+
        if (.not.half) then
          mat = 0d0
          ising = 0
          itrip = 0
          do icol = 1, ndim
            if (abs(map(icol)).eq.icol) then
              if (map(icol).gt.0) then
                ising = ising + 1
                mat(1:ndim,icol) = scr(1:ndim,ising)
              else
                itrip = itrip + 1
                mat(1:ndim,icol) = scr(1:ndim,nsing+itrip)
              end if
            else if (abs(map(icol)).gt.icol) then
              sig = dble(sign(1,map(icol)))
              ising = ising + 1
              mat(1:ndim,icol) = fac*scr(1:ndim,ising)
              mat(1:ndim,abs(map(icol))) = sig*fac*scr(1:ndim,ising)
            else
              sig = - dble(sign(1,map(icol)))
              itrip = itrip + 1
              mat(1:ndim,abs(map(icol))) = mat(1:ndim,abs(map(icol)))
     &                           + fac*scr(1:ndim,nsing+itrip)
              mat(1:ndim,icol) = mat(1:ndim,icol)
     &                           + sig*fac*scr(1:ndim,nsing+itrip)
            end if
          end do
          deallocate(scr)
        end if
c dbg   can be used to set 1
c        mat(1:ndim,1:ndim) = 0d0
c        do iline = 1, ndim
c          mat(iline,iline) = 1d0
c        end do
c dbgend

        if (ntest.ge.100) then
          write(luout,*) 'output matrix:'
          call wrtmat2(mat,ndim,ndim,ndim,ndim)
        end if
      end if

      return
      end
