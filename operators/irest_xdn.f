*----------------------------------------------------------------------*
      function irest_xdn(ixdn,irest,hpvxgas,ngas)
*----------------------------------------------------------------------*
*     return (dependent on ixdn) 1: e(x)citation, 2: (d)e-excitation, or
*     3: (n)on-excitation part of restriction
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer ::
     &     irest_xdn(2,ngas,2,2)

      integer, intent(in) ::
     &     ixdn, ngas, hpvxgas(ngas),
     &     irest(2,ngas,2,2)

      integer ::
     &     ica, ihpv, igas, iscr(2,ngas,2,2)

      iscr(1:2,1:ngas,1:2,1:2) = 0
      do igas = 1, ngas
        ihpv = hpvxgas(igas)
        if (ixdn.eq.3.and.ihpv.ne.ivale) cycle
        if (ixdn.ne.3.and.ihpv.eq.ivale) cycle
        ica = 1
        if (ixdn.eq.1.and.ihpv.eq.ihole .or.
     &      ixdn.eq.2.and.(ihpv.eq.ipart.or.ihpv.eq.iextr)) ica = 2
        iscr(1:2,igas,ica,1:2) = irest(1:2,igas,ica,1:2)
        if (ixdn.eq.3)
     &       iscr(1:2,igas,2,1:2) = irest(1:2,igas,2,1:2)
      end do

      irest_xdn = iscr

      return
      end
