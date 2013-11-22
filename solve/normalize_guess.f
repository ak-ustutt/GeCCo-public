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
     &     ifree, nwfpar(nopt), ioff, jopt
      real(8) ::
     &     xnrm
      real(8), pointer ::
     &     xbuf1(:), xbuf2(:)

      real(8), external ::
     &     ddot, da_ddot

      if (ntest.ge.100) then
        write(lulog,*) '-------------------------'
        write(lulog,*) ' normalize_guess at work '
        write(lulog,*) '-------------------------'
      end if

      nwfpar(1:nopt) = opti_info%nwfpar(1:nopt)

      ifree = mem_setmark('normalize_guess')
      ifree = mem_alloc_real(xbuf1,maxval(nwfpar(1:nopt)),'xbuf1')
      ifree = mem_alloc_real(xbuf2,nwfpar(iopt),'xbuf2')
      if (use_s(iopt)) then
        call vec_from_da(ff_v,irecv,xbuf1,nwfpar(iopt))
        call vec_from_da(ff_met(iopt)%fhand,irecmet,xbuf2,nwfpar(iopt))
        xnrm = ddot(nwfpar(iopt),xbuf1,1,xbuf2,1)
      else
        xnrm = da_ddot(ff_v,irecv,
     &                 ff_v,irecv,
     &                 nwfpar(iopt),
     &                 xbuf1,xbuf2,
     &                 nwfpar(iopt))
      end if
      xnrm = sqrt(xnrm)

      if (ntest.ge.100)
     &     write(lulog,*) 'Normalizing root with norm ',xnrm

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
