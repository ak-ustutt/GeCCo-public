*----------------------------------------------------------------------*
      subroutine normalize_guess(
     &     ff_v,irecv,ff_mvp,irecmvp,ff_met,irecmet,
     &     use_s,iopt,nopt,opti_info)
*----------------------------------------------------------------------*
*
*     Normalizes a trial vector and scales the coreesponding
*     matrix vector product (and metric vector product) accordingly.
*
*     matthias, mar 2013
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_file_array.h'
      include 'def_optimize_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     irecv, irecmvp, irecmet, iopt, nopt
      logical ::
     &     use_s(nopt)
      type(filinf) ::
     &     ff_v
      type(file_array) ::
     &     ff_mvp(nopt), ff_met(nopt)
      type(optimize_info), intent(in) ::
     &     opti_info
      
      integer ::
     &     ifree, nsec, isec, stsec, ndsec, nwfpar(nopt), ioff, jopt
      real(8) ::
     &     xnrm
      real(8), pointer ::
     &     xbuf1(:), xbuf2(:), signsec(:)
      integer, pointer ::
     &     nwfpsec(:), idstsec(:)

      real(8), external ::
     &     ddot, da_ddot

      if (ntest.ge.100) then
        write(luout,*) '-------------------------'
        write(luout,*) ' normalize_guess at work '
        write(luout,*) '-------------------------'
      end if

      nwfpar(1:nopt) = opti_info%nwfpar(1:nopt)

      ifree = mem_setmark('normalize_guess')
      ifree = mem_alloc_real(xbuf1,maxval(nwfpar(1:nopt)),'xbuf1')
      ifree = mem_alloc_real(xbuf2,nwfpar(iopt),'xbuf2')
      if (use_s(iopt)) then
        nsec = opti_info%nsec(iopt)
        ioff = sum(opti_info%nsec(1:iopt))-nsec
        nwfpsec => opti_info%nwfpsec(ioff+1:ioff+nsec)
        idstsec => opti_info%idstsec(ioff+1:ioff+nsec)
        signsec => opti_info%signsec(ioff+1:ioff+nsec)
        call vec_from_da(ff_v,irecv,xbuf1,nwfpar(iopt))
        call vec_from_da(ff_met(iopt)%fhand,irecmet,xbuf2,nwfpar(iopt))
        xnrm = 0d0
        do isec = 1, nsec
          xnrm = xnrm + signsec(isec)
     &         * ddot(nwfpsec(isec),xbuf1(idstsec(isec)),1,
     &                xbuf2(idstsec(isec)),1)
        end do
      else
        xnrm = da_ddot(ff_v,irecv,1,
     &                 ff_v,irecv,1,
     &                 nwfpar(iopt),
     &                 xbuf1,xbuf2,
     &                 nwfpar(iopt))
      end if
      xnrm = sqrt(xnrm)

      if (ntest.ge.100)
     &     write(luout,*) 'Normalizing root with norm ',xnrm

      call da_sccpvec(ff_v,irecv,
     &                ff_v,irecv,
     &                1d0/xnrm,nwfpar(iopt),
     &                xbuf1,nwfpar(iopt))
      do jopt = 1, nopt
        call da_sccpvec(ff_mvp(jopt)%fhand,irecmvp,
     &                  ff_mvp(jopt)%fhand,irecmvp,
     &                  1d0/xnrm,nwfpar(jopt),
     &                  xbuf1,nwfpar(jopt))
        if (use_s(jopt))
     &     call da_sccpvec(ff_met(jopt)%fhand,irecmet,
     &                  ff_met(jopt)%fhand,irecmet,
     &                  1d0/xnrm,nwfpar(jopt),
     &                  xbuf1,nwfpar(jopt))
      end do
      ifree = mem_flushmark()

      return

      end
