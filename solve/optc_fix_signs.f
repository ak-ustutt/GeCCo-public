*----------------------------------------------------------------------*
      subroutine optc_fix_signs(xrsnrm,xvec,
     &     ff_sbsp,ndim_sbsp,iord_sbsp,
     &     ff_vec,irecvec,
     &     opti_info,iopt,
     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2)
*----------------------------------------------------------------------*
*
*  Matthias' sign fix for internal contraction stuff (as collected
*     from evpc_core
*     
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_optimize_info.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ndim_sbsp,iord_sbsp(ndim_sbsp),irecvec
      real(8), intent(in) ::
     &     xvec(ndim_sbsp)
      real(8), intent(inout) ::
     &     xrsnrm
      type(filinf) ::
     &     ff_vec,ff_sbsp
      integer, intent(in) ::
     &     nincore, nwfpar, lenbuf, iopt
      real(8), intent(inout) ::
     &     xbuf1(lenbuf), xbuf2(lenbuf)
      type(optimize_info) ::
     &     opti_info


      integer ::
     &     ioff, nsec, idx, irec, isec
      integer, pointer ::
     &     nwfpsec(:), idstsec(:)
      real(8), pointer ::
     &     signsec(:)

      real(8), external ::
     &     ddot
      


      nsec = opti_info%nsec(iopt)
c dbg
      print *,'nsec = ',nsec
c dbg
      if (nsec.le.1) return

      ioff = sum(opti_info%nsec(1:iopt))-nsec
      nwfpsec => opti_info%nwfpsec(ioff+1:ioff+nsec)
      idstsec => opti_info%idstsec(ioff+1:ioff+nsec)
      signsec => opti_info%signsec2(ioff+1:ioff+nsec)
      do irec = 1, ndim_sbsp
        idx = iord_sbsp(irec)
        if (xvec(idx).eq.0d0) cycle
        call vec_from_da(ff_sbsp,irec,xbuf2,
     &                           nwfpar)
        do isec = 1, nsec
          xbuf1(idstsec(isec):idstsec(isec)+nwfpsec(isec)-1) = 
     &         xbuf1(idstsec(isec):idstsec(isec)+nwfpsec(isec)-1)
     &         +signsec(isec)*xvec(idx)*
     &         xbuf2(idstsec(isec):idstsec(isec)+nwfpsec(isec)-1)
        end do
      end do
      if (opti_info%typ_prc(iopt).ne.optinf_prc_traf)
     &     xrsnrm = sqrt(ddot(nwfpar,xbuf1,1,xbuf1,1))
      call vec_to_da(ff_vec,irecvec,xbuf1,nwfpar)

      return
      end
