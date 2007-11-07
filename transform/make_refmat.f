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
     &     idxcmo, idxdao, igas, isym, ngas, nsym, norb, nao
      integer, pointer ::
     &     hpvxgas(:), mostnd(:,:,:), nbas(:)


      hpvxgas => orb_info%ihpvgas
      mostnd => orb_info%mostnd
      nbas   => orb_info%nbas
      
      ngas = orb_info%ngas
      nsym = orb_info%nsym

      idxcmo = 1
      do igas = 1, ngas
        if (hpvxgas(igas).eq.ivale)
     &       call quit(1,'make_refmat',
     &                   'not prepared for valence spaces')
        idxdao = 1
        do isym = 1, nsym
          norb = mostnd(2,isym,igas)-mostnd(1,isym,igas)+1
          nao  = nbas(isym)

          ! do not put outside this loop, as we need to
          ! update idxcmo correctly!
          if (hpvxgas(igas).eq.ihole.and.nao*norb.gt.0) then

            call dgemm('n','t',nao,nao,norb,
     &           2d0,cmo(idxcmo),nao,
     &               cmo(idxcmo),nao,
     &           1d0,dao(idxdao),nao)
          end if
            
          idxdao = idxdao + nao*nao
          idxcmo = idxcmo + nao*norb
        end do
      end do

      return
      end
