*----------------------------------------------------------------------*
      subroutine make_refmat(dao,cmo,orb_info)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_orbinf.h'

      real(8), intent(out) ::
     &     dao(*)
      type(orbinf), intent(in), target ::
     &     orb_info
      real(8), intent(in) ::
     &     cmo(*)
      
      integer ::
     &     idxcmo, idxdao, igas, isym, ngas, nsym, norb, nao, nspin
      integer, pointer ::
     &     hpvxgas(:,:), mostnd(:,:,:), nbas(:)

      real(8) ::
     &     fac

      hpvxgas => orb_info%ihpvgas
      mostnd => orb_info%mostnd
      nbas   => orb_info%nbas
      
      ngas  = orb_info%ngas
      nspin = orb_info%nspin
      nsym  = orb_info%nsym

      idxcmo = 1
      do igas = 1, ngas
! COMMENTING OUT NEXT THREE LINES --- SHOULD NOT BE A PROBLEM @pradipta
!       if (hpvxgas(igas,1).eq.ivale)
!    &       call quit(1,'make_refmat',
!    &                   'not prepared for valence spaces')
        idxdao = 1
        do isym = 1, nsym
          norb = mostnd(2,isym,igas)-mostnd(1,isym,igas)+1
          nao  = nbas(isym)

          fac = 0d0
          if (hpvxgas(igas,1).eq.IHOLE)
     &         fac = 2d0        ! doubly occupied ?
          if (nspin.eq.2.and.hpvxgas(igas,2).ne.IHOLE)
     &         fac = 1d0        ! only singly occupied !

          ! do not put outside this loop, as we need to
          ! update idxcmo correctly!
          if (fac.gt.0d0.and.nao*norb.gt.0) then

            call dgemm('n','t',nao,nao,norb,
     &           fac,cmo(idxcmo),nao,
     &               cmo(idxcmo),nao,
     &           1d0,dao(idxdao),nao)
          end if
            
          idxdao = idxdao + nao*nao
          idxcmo = idxcmo + nao*norb
        end do
      end do

      return
      end
