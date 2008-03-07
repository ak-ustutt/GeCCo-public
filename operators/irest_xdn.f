*----------------------------------------------------------------------*
      function irest_xdn(ixdn,irest,hpvxgas,ngas,nspin)
*----------------------------------------------------------------------*
*     return (dependent on ixdn) 1: e(x)citation, 2: (d)e-excitation, or
*     3: (n)on-excitation part of restriction
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer ::
     &     nspin,
     &     irest_xdn(2,ngas,2,2,nspin)

      integer, intent(in) ::
     &     ixdn, ngas, hpvxgas(ngas),
     &     irest(2,ngas,2,2,nspin)

      integer ::
     &     ica, ihpv, igas, ispin, iscr(2,ngas,2,2,nspin)

      iscr(1:2,1:ngas,1:2,1:2,1:nspin) = 0
      do ispin = 1, nspin
        do igas = 1, ngas
          ihpv = hpvxgas(igas)
          if (ixdn.eq.3.and.ihpv.ne.ivale) cycle
          if (ixdn.ne.3.and.ihpv.eq.ivale) cycle
          ica = 1
          if (ixdn.eq.1.and.ihpv.eq.ihole .or.
     &      ixdn.eq.2.and.(ihpv.eq.ipart.or.ihpv.eq.iextr)) ica = 2
          iscr(1:2,igas,ica,1:2,ispin) = irest(1:2,igas,ica,1:2,ispin)
          if (ixdn.eq.3)
     &         iscr(1:2,igas,2,1:2,ispin) = irest(1:2,igas,2,1:2,ispin)
        end do
      end do

      irest_xdn = iscr

      return
      end
