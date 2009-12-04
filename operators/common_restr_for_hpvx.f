      subroutine common_restr_for_hpvx(rst_res,rst1,rst2,hpvx,
     &                             hpvxgas,ngas)
      ! 
      ! find max restriction that describes the common paths
      ! through an occupation graph
      !
      !  e.g.  0 0  2 2   means   +
      !                           |   
      !                           +-+-+
      !        0 1  2 2   means   +        +-+
      !                           |     +    |
      !                           +-+-+      +-+
      !
      !        the max. common occ. graph is   0 0  2 2
      !
      implicit none

      integer, parameter ::
     &     ntest = 00
      character(len=21), parameter ::
     &     i_am = 'common_restr_for_hpvx'

      integer, intent(in) ::
     &     ngas, hpvx, hpvxgas(ngas),
     &     rst1(2,ngas),rst2(2,ngas)
      integer, intent(out) ::
     &     rst_res(2,ngas)

      integer ::
     &     igas

      ! simple version, maybe some more logic will be neccessary
      do igas = 1, ngas
c dbg
        print *,'>> igas, hpvxgas, hpvx: ',igas,hpvxgas(igas),hpvx,
     &       hpvxgas(igas).ne.hpvx
c dbg
        if (hpvxgas(igas).ne.hpvx) cycle
        rst_res(1,igas) = max(rst1(1,igas),rst2(1,igas))
        rst_res(2,igas) = min(rst1(2,igas),rst2(2,igas))
      end do

      end
