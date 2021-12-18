      logical function is_proper_restr_for_hpvx(rst,hpvx,
     &                                          hpvxgas,ngas)
      implicit none

      integer, parameter ::
     &     ntest = 00
      character(len=24), parameter ::
     &     i_am = 'is_proper_restr_for_hpvx'

      integer, intent(in) ::
     &     ngas, hpvx, hpvxgas(ngas),
     &     rst(2,ngas)

      integer ::
     &     igas
      logical ::
     &     ok
      
      ! simple version, maybe some more logic will be neccessary
      ok = .true.
      do igas = 1, ngas
        if (hpvxgas(igas).ne.hpvx) cycle
        ok = ok.and.rst(1,igas).ge.0
        ok = ok.and.rst(2,igas).ge.0
        ok = ok.and.rst(1,igas).le.rst(2,igas)
      end do

      is_proper_restr_for_hpvx = ok

      end
